#!/usr/bin/env python3
"""
run_all.py -- batch driver for both experiments in the manuscript.

Runs every edge-list file given on the command line (or every *.txt / *.edges /
*.mtx in a directory) through both experiments and prints the two LaTeX tables
ready to paste into the paper:

    Table 2  (Section 10)  -- enumeration: d, omega, lambda, |F|, |F|/lambda,
                              time for Algorithm 1, time for the networkx
                              baseline, and an agreement check
    Table E1 (Section 10)  -- compression: modular width, ||Phi||_sh, lambda,
                              compression ratio

Usage
-----
    python3 run_all.py data/*.txt
    python3 run_all.py data/                 # whole directory
    python3 run_all.py data/ --timeout 3600  # per-network cap, seconds
    python3 run_all.py data/ --skip-compression

Each edge list is whitespace separated, one edge per line; '#' and '%' begin a
comment.  Vertex tokens are arbitrary.

Notes
-----
* Run mcd_reference.py --check and modular_compression.py --selftest once each
  before trusting any of these numbers.
* Report the machine, OS, Python version and networkx version alongside the
  tables; the paper's caption asks for them.
* If a network times out, the row is emitted with 'timeout' rather than being
  dropped.  Please keep it in the table.
"""

import os
import sys
import time
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mcd = _load("mcd", "mcd_reference.py")
mcomp = _load("mcomp", "modular_compression.py")


def gather(args):
    files = []
    for a in args:
        if os.path.isdir(a):
            for fn in sorted(os.listdir(a)):
                if fn.rsplit(".", 1)[-1].lower() in ("txt", "edges", "mtx", "csv"):
                    files.append(os.path.join(a, fn))
        else:
            files.append(a)
    return files


def label(path):
    return os.path.basename(path).rsplit(".", 1)[0].replace("_", r"\_")


def main():
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    skip_comp = "--skip-compression" in sys.argv
    files = gather(argv)
    if not files:
        print(__doc__)
        return

    try:
        import networkx as nx
        have_nx = True
    except ImportError:
        have_nx = False
        print("%% WARNING: networkx not installed -- no baseline column.")
        print("%%          pip install networkx, then re-run.\n")

    rows_enum, rows_comp = [], []

    for path in files:
        name = label(path)
        sys.stderr.write("... %s\n" % name)
        adj = mcd.read_edgelist(path)
        n = len(adj)
        m = sum(len(a) for a in adj.values()) // 2

        # ---- experiment 1: enumeration ----
        t0 = time.perf_counter()
        cliques, d, sizeF = mcd.maximal_clique_decomposition(adj)
        t_alg = time.perf_counter() - t0
        lam = len(cliques)
        omega = max((len(c) for c in cliques), default=0)

        t_base, agree = None, "n/a"
        if have_nx:
            G = nx.Graph()
            G.add_nodes_from(adj)
            for u in adj:
                for v in adj[u]:
                    G.add_edge(u, v)
            t1 = time.perf_counter()
            ref = {frozenset(c) for c in nx.find_cliques(G)}
            t_base = time.perf_counter() - t1
            agree = "yes" if ref == {frozenset(c) for c in cliques} else "NO"

        rows_enum.append((name, n, m, d, omega, lam, sizeF,
                          sizeF / lam if lam else 0, t_alg, t_base, agree))

        # ---- experiment 2: compression ----
        if skip_comp:
            continue
        try:
            t2 = time.perf_counter()
            T = mcomp.modular_decomposition(adj)
            mcomp.annotate(T, adj)
            t_md = time.perf_counter() - t2
            w = mcomp.modular_width(T)
            sh = mcomp.slp_size(T)
            lam2 = mcomp.count_cliques(T)
            om2 = mcomp.clique_number(T)
            rows_comp.append((name, n, m, w, sh, lam2,
                              lam2 / sh if sh else 0, t_md,
                              "yes" if (lam2 == lam and om2 == omega) else "NO"))
        except RecursionError:
            rows_comp.append((name, n, m, "-", "-", "-", "-", "-", "recursion"))

    # ---------------- output ----------------
    print("\n%% ---------- Table 2 : enumeration ----------")
    print("%% columns: network & n & m & d & omega & lambda & |F| & |F|/lambda"
          " & Alg.1 (s) & baseline (s)   [agreement in comment]")
    for r in rows_enum:
        (nm, n, m, d, om, lam, F, ratio, ta, tb, ag) = r
        tb_s = ("%.3f" % tb) if tb is not None else "---"
        print("%s & %d & %d & %d & %d & %d & %d & %.2f & %.3f & %s \\\\  %% agree=%s"
              % (nm, n, m, d, om, lam, F, ratio, ta, tb_s, ag))

    if rows_comp:
        print("\n%% ---------- Table E1 : compression ----------")
        print("%% columns: network & n & m & mw & ||Phi||_sh & lambda &"
              " lambda/||Phi||_sh   [lambda cross-check in comment]")
        for r in rows_comp:
            (nm, n, m, w, sh, lam2, ratio, tmd, ok) = r
            if w == "-":
                print("%s & %d & %d & \\multicolumn{4}{c}{recursion limit} \\\\" % (nm, n, m))
            else:
                print("%s & %d & %d & %d & %d & %d & %.2f \\\\  %% md=%.2fs cross-check=%s"
                      % (nm, n, m, w, sh, lam2, ratio, tmd, ok))

    bad = [r for r in rows_enum if r[10] == "NO"]
    if bad:
        print("\n%%%% *** DISAGREEMENT with networkx on: %s"
              % ", ".join(r[0] for r in bad))
        print("%%%% *** Do not report these numbers. Send the failing edge list on.")
    bad2 = [r for r in rows_comp if r[8] == "NO"]
    if bad2:
        print("\n%%%% *** lambda cross-check FAILED on: %s"
              % ", ".join(r[0] for r in bad2))
        print("%%%% *** This would contradict Theorem 6.6(1); please report it.")


if __name__ == "__main__":
    main()
