#!/usr/bin/env python3
"""
Reference implementation of the maximal-clique decomposition algorithm
(Algorithm 1 of "An Algorithm for the Maximal Clique Problem of Graphs
using Algebraic Techniques").

The point of this file is that it realises the complexity bound proved in
Theorem 8.3: monomials are stored as integer bit masks over the forward
neighbourhood N^r(v_i), so that intersection, union and inclusion of
monomials are single machine operations, and the maximality tests are the
ones used in the proof.  This is the only difference from a naive symbolic
implementation, and it is worth several orders of magnitude.

Usage
-----
    python3 mcd_reference.py edgelist.txt          # run and report statistics
    python3 mcd_reference.py edgelist.txt --check  # also verify against networkx

The edge list is a whitespace-separated file, one edge per line; lines
beginning with '#' or '%' are ignored.  Vertices may be arbitrary tokens.

Outputs the table columns used in Section 9 of the paper:
    |V|, |E|, d, omega, lambda, #omega, |F|, time.
"""

import sys
import time
from collections import defaultdict


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------
def read_edgelist(path):
    adj = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#%":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            u, v = parts[0], parts[1]
            if u != v:
                adj[u].add(v)
                adj[v].add(u)
    for v in list(adj):
        adj[v] = adj[v]
    return dict(adj)


# --------------------------------------------------------------------------
# degeneracy ordering (Matula-Beck, O(n + m))
# --------------------------------------------------------------------------
def degeneracy_ordering(adj):
    """Return (order, d).  In `order`, every vertex has at most d neighbours
    that follow it."""
    nodes = list(adj)
    deg = {v: len(adj[v]) for v in nodes}
    maxdeg = max(deg.values()) if deg else 0
    bucket = [set() for _ in range(maxdeg + 1)]
    for v in nodes:
        bucket[deg[v]].add(v)

    order, removed, d, i = [], set(), 0, 0
    for _ in range(len(nodes)):
        i = 0
        while i < len(bucket) and not bucket[i]:
            i += 1
        d = max(d, i)
        v = bucket[i].pop()
        order.append(v)
        removed.add(v)
        for u in adj[v]:
            if u not in removed:
                bucket[deg[u]].discard(u)
                deg[u] -= 1
                bucket[deg[u]].add(u)
    return order, d


# --------------------------------------------------------------------------
# Algorithm 1
# --------------------------------------------------------------------------
def maximal_clique_decomposition(adj, want_F=True):
    """Return (cliques, d, size_of_F).

    `cliques` is the list of maximal cliques, i.e. the monomials of mu(g).
    `size_of_F` is |F| of Proposition 7.4, the number of monomials before
    the Phase 2 reduction; it is returned for the experimental table.
    """
    order, d = degeneracy_ordering(adj)
    pos = {v: i for i, v in enumerate(order)}
    n = len(order)

    # forward neighbourhoods
    fwd = [[u for u in adj[v] if pos[u] > pos[v]] for v in order]

    cliques = []
    size_F = 0

    for i in range(n):
        vi = order[i]
        R = fwd[i]                       # local universe, |R| <= d
        if not R:
            # h_i = epsilon; the only candidate clique is {v_i} itself
            size_F += 1
            if not _dominated(adj, pos, vi, (), i):
                cliques.append((vi,))
            continue

        idx = {u: b for b, u in enumerate(R)}
        Rset = set(R)

        # local adjacency masks inside g[N^r(v_i)]
        nb = [0] * len(R)
        for u in R:
            mask = 0
            for w in adj[u]:
                if w in Rset:
                    mask |= 1 << idx[w]
            nb[idx[u]] = mask

        # ---- Phase 1: build mu(g[N^r(v_i)]) by adding local vertices ----
        H = [0]                          # list of monomials; 0 is epsilon
        for b in range(len(R)):
            bit = 1 << b
            Nb = nb[b]
            keep, cand = [], []
            for m in H:
                if m & ~Nb == 0:         # V_m subset of N(v_b): extend  (rule 1a)
                    keep.append(m | bit)
                else:                    # keep and generate              (rule 1b)
                    keep.append(m)
                    cand.append((m & Nb) | bit)
            if cand:
                prefix = (1 << (b + 1)) - 1   # vertices of R present so far
                seen = set()
                for c in cand:
                    if c in seen:
                        continue
                    seen.add(c)
                    # c is maximal in the current local graph iff the common
                    # neighbourhood of V_c, taken inside the prefix, lies in V_c
                    common = prefix
                    cc = c
                    while cc:
                        low = cc & -cc
                        common &= nb[low.bit_length() - 1]
                        cc ^= low
                    if common & ~c == 0:
                        keep.append(c)
            H = keep

        # ---- Phase 2: keep only the globally maximal ones ----
        size_F += len(H)
        for m in H:
            verts = tuple(R[b] for b in range(len(R)) if m >> b & 1)
            if not _dominated(adj, pos, vi, verts, i):
                cliques.append((vi,) + verts)

    return cliques, d, (size_F if want_F else None)


def _dominated(adj, pos, vi, verts, i):
    """True iff some u in N^l(v_i) is adjacent to every vertex of `verts`
    (Theorem 7.2(3)).  Then {v_i} u verts is not a maximal clique."""
    cand = None
    for w in (vi,) + verts:
        s = {x for x in adj[w] if pos[x] < i}
        cand = s if cand is None else (cand & s)
        if not cand:
            return False
    return bool(cand)


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    path = sys.argv[1]
    check = "--check" in sys.argv

    adj = read_edgelist(path)
    n = len(adj)
    m = sum(len(a) for a in adj.values()) // 2

    t0 = time.perf_counter()
    cliques, d, size_F = maximal_clique_decomposition(adj)
    elapsed = time.perf_counter() - t0

    omega = max((len(c) for c in cliques), default=0)
    n_omega = sum(1 for c in cliques if len(c) == omega)
    biggest = min((sorted(c) for c in cliques if len(c) == omega), default=[])

    print("file            : %s" % path)
    print("|V|             : %d" % n)
    print("|E|             : %d" % m)
    print("degeneracy d    : %d" % d)
    print("omega           : %d" % omega)
    print("lambda = |mu(g)|: %d" % len(cliques))
    print("#omega          : %d" % n_omega)
    print("|F|             : %d   (|F|/lambda = %.3f)"
          % (size_F, size_F / len(cliques) if cliques else 0))
    print("a maximum clique: %s" % biggest)
    print("time (s)        : %.4f" % elapsed)

    if check:
        try:
            import networkx as nx
        except ImportError:
            print("\n[--check skipped: networkx not installed]")
            return
        G = nx.Graph()
        G.add_nodes_from(adj)
        for u in adj:
            for v in adj[u]:
                G.add_edge(u, v)
        t1 = time.perf_counter()
        ref = {frozenset(c) for c in nx.find_cliques(G)}
        t_ref = time.perf_counter() - t1
        ours = {frozenset(c) for c in cliques}
        print("\nnetworkx find_cliques : %d cliques in %.4f s" % (len(ref), t_ref))
        print("agreement             : %s" % ("OK" if ours == ref else "MISMATCH"))


if __name__ == "__main__":
    main()
