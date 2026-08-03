# Experiment package — for Dr. Murali Krishna Enduri

Everything needed to produce the two tables in Section 10 of
*Maximal Cliques as a Normal Form in the Algebra of Graphs*.

Pure Python 3, no dependencies except `networkx` for the baseline column.

```
pip install networkx
```

## Files

| file | what it does |
|---|---|
| `mcd_reference.py` | Algorithm 1, bit-mask implementation. Produces Table 2. |
| `modular_compression.py` | Modular decomposition, `mw(g)`, `‖Φ‖_sh`, and λ, ω by semiring evaluation. Produces Table E1. |
| `run_all.py` | Runs both over a directory and prints both LaTeX tables. |
| `demo/` | Two tiny graphs for a smoke test. |

## Step 1 — validate before trusting anything

```
python3 mcd_reference.py --check        # (on any edge list)
python3 modular_compression.py --selftest
```

Expected:

```
randomized validation: 800 trials, 0 mismatches
modularity of children : 400 graphs, 0 failures
count via (N,+,x)      : 400 graphs, 0 mismatches
omega  via (N,max,+)   : 400 graphs, 0 mismatches
cocktail-party n= 24 : mw=2  ||Phi||= 60  lambda=  4096  ratio=68.3
```

If any line reports a failure, **stop and tell me** — it means a theorem in the
paper is wrong, not just the code.

## Step 2 — smoke test

```
python3 run_all.py demo/
```

Should print two LaTeX tables with `agree=yes` and `cross-check=yes` on both rows.

## Step 3 — the real runs

Put the eleven networks from the DMAA submission in a directory as edge lists
(whitespace-separated pairs, one per line, `#`/`%` comments ignored):

```
python3 run_all.py data/ > tables.tex
```

Paste the rows into Table 2 and Table E1. Then please record, for the captions:
machine and CPU, OS, Python version, `networkx` version.

Timing: `powerGrid` (n=4941) takes about a minute for the modular decomposition
on a modern laptop. Nothing should take hours. If something does, that is itself
worth reporting.

## What each table is for

**Table 2 — enumeration.** Columns `d`, `ω`, `λ`, `|F|`, `|F|/λ`, our time,
baseline time. Two things matter here:

- **The baseline column is not optional.** DMAA rejected the previous version
  partly because there was no comparison. Report it whatever it says. We expect
  to be *slower* than tuned Bron–Kerbosch by a small constant factor; saying so
  is much better received than omitting it.
- **`|F|/λ` tests Proposition 7.4**, which guarantees `|F| ≤ (d+1)λ`. In my
  trials it stayed near 1.

**Table E1 — compression.** This one matters more. Theorem 6.6 is the paper's
headline claim, and the DMAA referee's central objection was that the algebraic
viewpoint gave no new insight. A theorem with no measurement behind it is a weak
answer.

Columns `mw(g)`, `‖Φ‖_sh`, `λ`, `λ/‖Φ‖_sh`. The last is the compression factor.

**Please report it honestly whichever way it goes.** On uniform random sparse
graphs I get `mw ≈ n` (essentially prime, no modular structure) and a ratio
*below* 1 — no compression at all. Real networks contain many non-trivial modules
and may behave quite differently, but we do not know until you run it. If the
compression does not appear, that is a finding we need before a referee tells us,
and it should go in the paper with the observed `mw(g)` explaining why.

`run_all.py` also cross-checks that λ from the semiring evaluation matches λ from
actual enumeration. That is a direct empirical test of Theorem 6.6(1) on real
data and worth a sentence in the text.

## If something disagrees

`run_all.py` prints a loud comment block if either check fails. In that case send
back the failing edge list rather than the table — a disagreement means either a
bug or a false theorem, and we need to know which before submitting.

## One caution

Please do not fill any cell from the old DMAA table. Those timings came from a
prototype that manipulated symbolic expressions as strings; they were dominated
by interpreter overhead and contradicted our own complexity theorem. Every number
in both tables should come from a fresh run.
