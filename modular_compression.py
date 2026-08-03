#!/usr/bin/env python3
"""
modular_compression.py
======================

Measures the compression predicted by Theorem 6.6 of the manuscript
"Maximal Cliques as a Normal Form in the Algebra of Graphs".

For an input graph it computes:

    n, m          order and size
    mw            modular width  (largest number of children of a prime node)
    ||Phi||       size of the shared factored decomposition (straight-line
                  program), i.e. the total number of symbol occurrences in the
                  node formulas of the modular decomposition tree
    lambda        |mu(g)|, the number of maximal cliques, obtained WITHOUT
                  enumerating them, by evaluating Phi in the semiring (N,+,x)
    omega         clique number, by evaluating Phi in (N,max,+)
    ratio         lambda / ||Phi||   -- the compression factor

IMPORTANT.  The bound of Theorem 6.6 is a bound on the SHARED representation.
Each node's formula is emitted once and referred to by name; writing Phi out as
a literal tree, with each child's expression duplicated once per maximal clique
of the quotient that contains it, is exponentially larger in general.  This
script measures the shared size, which is what the theorem is about.

Usage
-----
    python3 modular_compression.py edgelist.txt
    python3 modular_compression.py --selftest     # correctness checks

The edge list is whitespace-separated, one edge per line; '#' and '%' start a
comment.  Vertices may be arbitrary tokens.
"""

import sys
import time
from itertools import combinations

# the decomposition recurses on the modular tree; large sparse graphs can be deep
sys.setrecursionlimit(200000)


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------
def read_edgelist(path):
    adj = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#%":
                continue
            p = line.split()
            if len(p) < 2:
                continue
            u, v = p[0], p[1]
            adj.setdefault(u, set())
            adj.setdefault(v, set())
            if u != v:
                adj[u].add(v)
                adj[v].add(u)
    return adj


# ---------------------------------------------------------------------------
# basic graph utilities on an explicit vertex subset
# ---------------------------------------------------------------------------
def components(adj, V):
    V = set(V)
    seen, out = set(), []
    for s in V:
        if s in seen:
            continue
        comp, stack = set(), [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.add(x)
            for y in adj[x] & V:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        out.append(comp)
    return out


def co_components(adj, V):
    """Connected components of the complement of g[V]."""
    V = set(V)
    seen, out = set(), []
    for s in V:
        if s in seen:
            continue
        comp, stack = set(), [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.add(x)
            for y in (V - adj[x] - {x}):
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        out.append(comp)
    return out


# ---------------------------------------------------------------------------
# maximal modules by partition refinement
# ---------------------------------------------------------------------------
def maximal_modules_not_containing(adj, V, v):
    """Maximal modules of g[V] that do not contain v (classical M(G,v)),
    computed by partition refinement."""
    V = set(V)
    rest = V - {v}
    P = [s for s in (adj[v] & rest, rest - adj[v]) if s]
    changed = True
    while changed:
        changed = False
        for x in V:
            newP = []
            for X in P:
                if x in X:
                    newP.append(X)
                    continue
                A = X & adj[x]
                B = X - adj[x]
                if A and B:
                    newP.append(A)
                    newP.append(B)
                    changed = True
                else:
                    newP.append(X)
            P = newP
    return P


def maximal_modular_partition(adj, V):
    """Partition of V into maximal strong modules, for g[V] with both g[V] and
    its complement connected (the prime case)."""
    V = set(V)
    v = next(iter(V))
    P1 = maximal_modules_not_containing(adj, V, v)
    # find the part containing v, using a pivot from outside it
    u = None
    for X in P1:
        if v not in X:
            u = next(iter(X))
            break
    if u is None:
        return [{x} for x in V]
    P2 = maximal_modules_not_containing(adj, V, u)
    Mv = None
    for X in P2:
        if v in X:
            Mv = X
            break
    if Mv is None:
        Mv = {v}
    parts = [Mv] + [X for X in P1 if not (X & Mv)]
    covered = set().union(*parts) if parts else set()
    for x in V - covered:          # safety net; should not trigger
        parts.append({x})
    return parts


# ---------------------------------------------------------------------------
# modular decomposition tree
# ---------------------------------------------------------------------------
class Node:
    __slots__ = ("kind", "children", "vertices", "quotient_mu", "formula_size")

    def __init__(self, kind, children, vertices):
        self.kind = kind              # 'leaf' | 'parallel' | 'series' | 'prime'
        self.children = children
        self.vertices = vertices
        self.quotient_mu = None       # list of index-sets, = mu(h_N)
        self.formula_size = 0


def bk_small(adj, V):
    """All maximal cliques of g[V]; V is small (a prime quotient)."""
    V = set(V)
    if not V:
        return [frozenset()]
    out = []

    def rec(R, P, X):
        if not P and not X:
            out.append(frozenset(R))
            return
        piv = max(P | X, key=lambda u: len(adj[u] & P))
        for w in list(P - adj[piv]):
            rec(R | {w}, P & adj[w], X & adj[w])
            P = P - {w}
            X = X | {w}

    rec(set(), set(V), set())
    return out


def modular_decomposition(adj, V=None):
    if V is None:
        V = set(adj)
    V = set(V)
    if len(V) == 1:
        return Node("leaf", [], V)
    comps = components(adj, V)
    if len(comps) > 1:
        ch = [modular_decomposition(adj, C) for C in comps]
        return Node("parallel", ch, V)
    cocomps = co_components(adj, V)
    if len(cocomps) > 1:
        ch = [modular_decomposition(adj, C) for C in cocomps]
        return Node("series", ch, V)
    parts = maximal_modular_partition(adj, V)
    if len(parts) <= 1:               # cannot decompose further
        return Node("prime", [modular_decomposition(adj, {x}) for x in V], V)
    ch = [modular_decomposition(adj, M) for M in parts]
    return Node("prime", ch, V)


def annotate(node, adj):
    """Attach mu(h_N) to each internal node and record its formula size."""
    for c in node.children:
        annotate(c, adj)
    p = len(node.children)
    if node.kind == "leaf":
        node.quotient_mu = None
        node.formula_size = 1
        return
    if node.kind == "parallel":
        node.quotient_mu = [frozenset([i]) for i in range(p)]
    elif node.kind == "series":
        node.quotient_mu = [frozenset(range(p))]
    else:
        # build the quotient graph on the children and take its maximal cliques
        reps = [next(iter(c.vertices)) for c in node.children]
        qadj = {i: set() for i in range(p)}
        for i, j in combinations(range(p), 2):
            if reps[j] in adj[reps[i]]:
                qadj[i].add(j)
                qadj[j].add(i)
        node.quotient_mu = bk_small(qadj, range(p))
    node.formula_size = sum(len(I) for I in node.quotient_mu)


def slp_size(node):
    """Total size of the shared (straight-line-program) representation."""
    return node.formula_size + sum(slp_size(c) for c in node.children)


def modular_width(node):
    w = 2
    if node.kind == "prime":
        w = max(w, len(node.children))
    for c in node.children:
        w = max(w, modular_width(c))
    return w


def count_cliques(node):
    """|mu(g)| by evaluating Phi in (N, +, x).  No enumeration."""
    if node.kind == "leaf":
        return 1
    vals = [count_cliques(c) for c in node.children]
    total = 0
    for I in node.quotient_mu:
        prod = 1
        for i in I:
            prod *= vals[i]
        total += prod
    return total


def clique_number(node):
    """omega(g) by evaluating Phi in (N, max, +)."""
    if node.kind == "leaf":
        return 1
    vals = [clique_number(c) for c in node.children]
    return max(sum(vals[i] for i in I) for I in node.quotient_mu)


# ---------------------------------------------------------------------------
# self-test
# ---------------------------------------------------------------------------
def _selftest():
    import random
    from itertools import chain

    def is_module(adj, V, M):
        M, V = set(M), set(V)
        for x in V - M:
            k = len(adj[x] & M)
            if k not in (0, len(M)):
                return False
        return True

    def brute_maximal_cliques(adj, V):
        return set(bk_small(adj, V))

    rnd = random.Random(5)
    bad_mod = bad_cnt = bad_om = 0
    for t in range(400):
        n = rnd.randint(1, 9)
        V = list(range(n))
        adj = {v: set() for v in V}
        for a, b in combinations(V, 2):
            if rnd.random() < rnd.random():
                adj[a].add(b)
                adj[b].add(a)
        T = modular_decomposition(adj)
        annotate(T, adj)
        # (a) every internal node's children are modules of the node's subgraph
        def check(node):
            ok = True
            if node.children:
                for c in node.children:
                    if not is_module(adj, node.vertices, c.vertices):
                        ok = False
                if set().union(*[c.vertices for c in node.children]) != node.vertices:
                    ok = False
            return ok and all(check(c) for c in node.children)
        if not check(T):
            bad_mod += 1
        # (b) counting and omega via the semiring evaluations
        true = brute_maximal_cliques(adj, V)
        if count_cliques(T) != len(true):
            bad_cnt += 1
        if clique_number(T) != max(len(c) for c in true):
            bad_om += 1
    print("modularity of children : 400 graphs, %d failures" % bad_mod)
    print("count via (N,+,x)      : 400 graphs, %d mismatches" % bad_cnt)
    print("omega  via (N,max,+)   : 400 graphs, %d mismatches" % bad_om)

    # cocktail-party graph: the extreme compression case
    for k in (3, 5, 8, 12):
        V = list(range(2 * k))
        adj = {v: set() for v in V}
        for a, b in combinations(V, 2):
            if a // 2 != b // 2:
                adj[a].add(b)
                adj[b].add(a)
        T = modular_decomposition(adj)
        annotate(T, adj)
        print("cocktail-party n=%3d : mw=%d  ||Phi||=%3d  lambda=%6d  ratio=%.1f"
              % (2 * k, modular_width(T), slp_size(T), count_cliques(T),
                 count_cliques(T) / slp_size(T)))


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    if sys.argv[1] == "--selftest":
        _selftest()
        return

    path = sys.argv[1]
    adj = read_edgelist(path)
    n = len(adj)
    m = sum(len(a) for a in adj.values()) // 2

    t0 = time.perf_counter()
    T = modular_decomposition(adj)
    annotate(T, adj)
    t_md = time.perf_counter() - t0

    w = modular_width(T)
    size = slp_size(T)
    t1 = time.perf_counter()
    lam = count_cliques(T)
    om = clique_number(T)
    t_eval = time.perf_counter() - t1

    print("file           : %s" % path)
    print("n              : %d" % n)
    print("m              : %d" % m)
    print("modular width  : %d" % w)
    print("||Phi|| (SLP)  : %d" % size)
    print("lambda         : %d      (by semiring evaluation, no enumeration)" % lam)
    print("omega          : %d" % om)
    print("compression    : lambda/||Phi|| = %.3f" % (lam / size if size else 0))
    print("time (decomp)  : %.3f s" % t_md)
    print("time (evaluate): %.4f s" % t_eval)
    print()
    print("LaTeX row:")
    print("%s & %d & %d & %d & %d & %d & %.2f \\\\" %
          (path.rsplit("/", 1)[-1].rsplit(".", 1)[0], n, m, w, size, lam,
           lam / size if size else 0))


if __name__ == "__main__":
    main()
