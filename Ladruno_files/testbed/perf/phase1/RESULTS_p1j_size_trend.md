# ADR-75 P1j — the factorization share RISES with model size, and the exponents say why

**Closes P1h's one open follow-up.** P1h measured the phase split at a single size (11.5k DOF)
and *predicted* the trend: factorization is superlinear in N while triangular solve is ~linear in
nnz(L), so the factorization share — and with it the value of every factorization-reuse lever —
should rise with N. This measures it over an order of magnitude.

- **Model:** the ADR-75 Lane-B cube — `LadrunoBrick` + `LadrunoJ2`, `Newton` + `LoadControl`,
  **15 steps**, n = 15/20/25/30/35 ⇒ **11,520 / 26,460 / 50,700 / 86,490 / 136,080 DOF**.
- **Binary:** this worktree at `6e5e1ba0b` (post-P1i, so `threads`/`nElem`/`nNode` are honest).
- **Method:** `p1j_size_trend.py`, **median of 2**, COARSE profiling, `MKL_NUM_THREADS=4`.
  Solves per run are identical between the two storage modes at every size (44/45/47/48/50), so
  no row is confounded by a convergence difference. Tip displacement agrees between modes to
  <1e-9 at every size.
- **Reproduce:** `OPS_PYD=<dist\bin> MKL_NUM_THREADS=4 python3.12 -S p1j_size_trend.py`, then
  `p1h_parse_split.py` over each h5, then `p1j_rollup.py`.

---

## 1. The trend

`-matrixType 2` @4T (the shipped configuration):

| n | nDOF | step s | solver s | %step | symb ms | factor ms | trisolve ms | **fac/(fac+tri)** | fac/tri |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 15 | 11,520 | 13.28 | 2.95 | 22.2 | 160 | 1,991 | 798 | **71.4%** | 2.5 |
| 20 | 26,460 | 37.77 | 10.46 | 27.7 | 114 | 8,267 | 2,080 | **79.9%** | 4.0 |
| 25 | 50,700 | 79.75 | 29.78 | 37.3 | 250 | 25,333 | 4,196 | **85.8%** | 6.0 |
| 30 | 86,490 | 170.97 | 77.86 | 45.5 | 401 | 70,255 | 7,201 | **90.7%** | 9.8 |
| 35 | 136,080 | 299.51 | 180.13 | 60.1 | 712 | 168,187 | 11,230 | **93.7%** | 15.0 |

`-matrixType 0` @4T:

| n | nDOF | step s | solver s | %step | symb ms | factor ms | trisolve ms | **fac/(fac+tri)** | fac/tri |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 15 | 11,520 | 15.61 | 4.59 | 29.4 | 79 | 3,809 | 702 | **84.4%** | 5.4 |
| 20 | 26,460 | 45.82 | 15.66 | 34.2 | 158 | 13,393 | 2,109 | **86.4%** | 6.3 |
| 25 | 50,700 | 102.66 | 50.21 | 48.9 | 281 | 45,531 | 4,398 | **91.2%** | 10.4 |
| 30 | 86,490 | 209.73 | 130.89 | 62.4 | 530 | 123,221 | 7,135 | **94.5%** | 17.3 |
| 35 | 136,080 | 290.25 | 204.63 | 70.5 | 673 | 195,974 | 7,988 | **96.1%** | 24.5 |

## 2. The mechanism — measured, not argued

Log-log slope against nDOF over all five points:

| phase | `-matrixType 2` | `-matrixType 2` (warm n=35, see §5) | `-matrixType 0` |
|---|---|---|---|
| `soe.factor` | N^1.793 | **N^1.645** | **N^1.654** |
| `soe.trisolve` | N^1.071 | **N^0.935** | **N^1.015** |
| `soe.symbolic` | N^0.654 | — | N^0.899 |
| solver total | N^1.662 | — | N^1.592 |
| step | N^1.262 | — | N^1.211 |
| **factor − trisolve gap** | **0.722** | **0.710** | **0.639** |

**Triangular solve is essentially linear (N^0.94–1.07); factorization is strongly superlinear
(N^1.65–1.79).** That gap *is* the trend, and a slope fitted over five points is a much stronger
claim than two endpoint percentages. Confirmed on both storage modes independently.

⚠ **The absolute exponents are soft; the GAP is not.** The `-matrixType 2` top point was measured
under a cold process heap (§5), which inflates both phases there — recomputing with a warm n=35
moves `soe.factor` from N^1.793 to N^1.645 and `soe.trisolve` from N^1.071 to N^0.935. But the
**gap moves only 0.722 → 0.710**, and the unsymmetric mode (order-invariant, §5) gives 0.639
independently. **Quote the gap (~0.64–0.72), or quote the exponents as ranges — never a single
exponent to three digits.**

`soe.symbolic` is **sublinear** and becomes irrelevant with size (1.2% → 0.24% of step).

## 3. What it does to the reuse-lever ceiling — P1h's "16%" was a floor

Factorization as a share of **step** (`-matrixType 2` @4T):

| nDOF | factor %step | trisolve %step | a perfect factorization-reuse lever buys |
|---:|---:|---:|---:|
| 11,520 | 15.0% | 6.0% | 1.18× |
| 50,700 | 31.8% | 5.3% | 1.47× |
| 136,080 | **56.2%** | **3.7%** | **2.28×** |

Two things move in opposite directions, and both favour the levers:

1. **The part reuse CAN skip grows** — 15.0% → 56.2% of step.
2. **The part it cannot touch SHRINKS** — triangular solve 6.0% → 3.7% of step.

So `-krylov`, the `factored` gate, ModifiedNewton and IMPL-EX are all worth materially more at
production scale than the 11.5k-DOF measurement suggested. P1h's explicit warning that its
16%-of-step figure was "a small-model floor, not a ceiling for production" is confirmed.

**Independent corroboration.** Extrapolating the fitted `step` and `solver` exponents to the
540,675 DOF of ADR-75b's G-L3 run predicts a solver share saturating at ~100%; G-L3 measured
**98.5%** on the cluster by a completely different route. Two unrelated measurements agree.

## 4. Method notes / caveats

- **The `fac/(fac+tri)` denominator excludes `soe.symbolic` deliberately.** Symbolic runs once per
  sparsity pattern while factor/trisolve run once per solve, so symbolic's share depends on the
  STEP COUNT, not on model size. Including it would let a step-count choice leak into a size
  trend. Steps are held constant at 15 across all sizes for the same reason.
- **Absolute small-N values carry several points of scatter.** P1h's n=15 `-matrixType 2` row
  gives `fac/(fac+tri)` = 78.1% where this sweep gives 71.4% (the `-matrixType 0` rows agree
  closely: 85.1% vs 84.4%). The difference is in `trisolve`, the smallest and noisiest number in
  the set. The **trend** spans 22 percentage points and is monotone on both storage modes, so it
  is far larger than this scatter — but do not quote any single row to three digits.
- **Median of 2, one box.** Threaded PARDISO is not byte-reproducible (P1f) and this box swings
  ±30% on background load. The headline metric is a *within-run ratio*, so box load largely
  cancels; absolute walls do not get that protection.
- **`-matrixType 2` ran FIRST at every size**, so its absolute `step`/`solver`/`factor` walls at
  n ≥ 30 carry the cold-heap penalty of §5. The ratio columns do not. **No sym-vs-unsym step
  comparison is made in this document** — the two modes were never measured under matched
  allocation conditions at large N.

## 5. The allocation-order artifact — and why the headline survives it

`formTangent` appeared to *invert* with size: symmetric assembly was faster at n=15 (8.4 s vs
9.0 s, consistent with P1h's `soe.addA` finding) but ~37% **slower** at n=35 (97.9 s vs 68.7 s).
Solve counts were identical (50/50), so it was not a convergence difference. Since
`-matrixType 2` always ran first within a size, the obvious suspect was allocation order, so it
was tested rather than assumed — one n=35 pair with the config order reversed:

| n=35 run | position | step s | formTangent s | solver s | **fac/(fac+tri)** |
|---|---|---:|---:|---:|---:|
| `-matrixType 2` | 1st | 310.31 / 288.70 | 97.86 / 91.35 | 187.44 / 172.81 | 93.68 / 93.81% |
| `-matrixType 2` | **2nd** | **189.55** | **64.65** | **108.89** | **93.50%** |
| `-matrixType 0` | 2nd | 282.97 / 297.53 | 68.69 / 69.46 | 198.03 / 211.24 | 96.03 / 96.13% |
| `-matrixType 0` | **1st** | **294.47** | **69.24** | **208.66** | **96.02%** |

**It is an artifact, and a big one — 1.58× on the whole step.** Symmetric run second is 189.6 s
against 299.5 s run first, with *both* assembly (−34%) and solve (−40%) collapsing. Unsymmetric is
indifferent to position (294.5 s vs ~290 s).

The asymmetry makes sense: whichever config must **grow** the process heap pays first-touch page
faults, and unsymmetric CSR is ~2× the entries — so symmetric-after-unsymmetric is the only
ordering that fits entirely inside memory the process has already touched. Unsymmetric never
benefits, because nothing before it allocated more.

**Two consequences.**

1. **The apparent assembly inversion is not real.** Measured warm, symmetric assembly is *faster*
   at n=35 too (64.65 s vs 69.24 s) — the same direction as n=15 and as P1h's `soe.addA` result.
   Had this gone unchecked, the published claim would have been the exact opposite of the truth.
2. **The headline is untouched.** `fac/(fac+tri)` reads **93.68 / 93.81 / 93.50%** for symmetric
   across two orderings and a 1.58× wall swing, and **96.03 / 96.13 / 96.02%** for unsymmetric —
   stable to ±0.15 pp. `factor %step` moves only 56.4 / 55.9 → 53.5. This is the empirical
   justification for choosing a within-run ratio as the metric instead of an absolute wall: it was
   picked to be robust to box load, and it turns out to be robust to a far larger effect that was
   not anticipated at all.

**Rule for future perf scripts:** either run one configuration per process, or randomize/rotate the
order and report the position. A multi-config loop in one process silently confounds storage mode
with allocation order at large N.
