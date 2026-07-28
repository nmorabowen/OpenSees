# ADR-75 P1h — the PARDISO solver-phase split, measured

**Closes the measurement half of the instrumentation item.**
[#667](https://github.com/nmorabowen/OpenSees/pull/667) added the brackets and CI proves they
*exist*; nothing published what the split *is*. This does.

- **Model:** Lane B — `laneB_model.py`, 15×15×15 `LadrunoBrick` + `LadrunoJ2` (3375 elements,
  4096 nodes, **11,520 DOF**), plastic, `Newton` + `LoadControl`, 15 steps, **44 solves**.
- **Binary:** this worktree's `dist/bin/opensees.pyd` at `e3e575549` + the P1h harness
  (all 7 brackets present; `tests/test_pardiso_solver.py` **8/8 pass** on it).
- **Method:** `./sweep_p1h.sh 3` — interleaved rounds, **median of 3**, COARSE profiling for the
  headline (see §5 for why not deep), MKL/OMP threads pinned per row.
- **Reproduce:** `./sweep_p1h.sh 3`, then `p1h_parse_split.py <h5> --json <json>`.

---

## 1. The split

Every `soe.*` figure is summed over **every site** in the tree, not read off `linearSolve`
(the `LEDGER_quirks` multi-site rule). Times in ms.

| config | step | solver | %step | symbolic | **factor** | trisolve | cgs | **%fac** | %tri |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| UmfPack @1T | 24199 | 15708 | 64.9 | 72.8 | **13692.9** | 1663.4 | — | **88.0** | 11.5 |
| Pardiso `-matrixType 0` @1T | 17843 | 8882 | 49.8 | 208.9 | **7761.0** | 912.0 | — | **87.4** | 10.3 |
| Pardiso `-matrixType 0` @4T | 13645 | 4298 | 31.6 | 163.2 | **3529.5** | 617.7 | — | **81.8** | 14.1 |
| Pardiso `-matrixType 2` @4T | 11451 | 2484 | 21.8 | 170.9 | **1764.5** | 495.8 | — | **73.4** | 20.1 |
| Pardiso `-matrixType 0 -krylov 6` @4T | 10840 | 1511 | 13.7 | 184.6 | **75.1** | 14.6 | 1154.9 | **5.0** | 1.0 |

`%fac` / `%tri` are shares **of the solver**, which is the number that bounds a reuse lever.
`solver = symbolic + factor + trisolve + cgs`.

## 2. Headline — factorization dominates, but its share is NOT a constant

**Factorization is ~88% of solver time single-threaded, and UmfPack and PARDISO agree on that to
within 0.6 pp** (88.0 vs 87.4). PARDISO's win is therefore a *faster factorization*, not a
different phase mix — the phases are in the same proportion, PARDISO just runs them ~1.76× faster.

The non-obvious result is what happens as the solver gets **better**:

| lever | factor | trisolve | effect on the mix |
|---|---:|---:|---|
| 1T → 4T (`-matrixType 0`) | 7761.0 → 3529.5 = **2.20×** | 912.0 → 617.7 = **1.48×** | %fac 87.4 → 81.8 |
| `-matrixType 0` → `2` (@4T) | 3529.5 → 1764.5 = **2.00×** | 617.7 → 495.8 = **1.25×** | %fac 81.8 → 73.4 |

**Both levers speed up factorization roughly twice as much as they speed up the triangular solve**,
so every improvement to the solver makes the *non*-factorization remainder a larger share. This is
the direct measurement of what P1e could only infer from a size trend ("threading erodes `-krylov`
— factorization parallelizes, triangular solves don't"): it is now a measured 2.20× vs 1.48×
rather than an inference from 1.22→1.30→1.51×.

**Consequence for the reuse levers.** At the production-recommended config (`-matrixType 2` @4T)
the solver is 21.8% of step and factorization is 73.4% of that = **16.0% of step**. That is the
hard ceiling on *any* factorization-reuse lever at this size — a perfect one that made
factorization free would buy 1.19× on the step. Triangular solve is 20.1% of solver = **4.4% of
step**, and no reuse lever touches it. The levers are worth having, but the headroom shrinks
precisely on the configurations you would actually ship.

**It also explains the 1→8 thread ceiling** (1.58×, ADR-75 headline): Amdahl *inside* the solver.
The part that threads is the part that is being made small.

## 3. `-krylov` — and the concrete proof that `soe.factor` understates

`-krylov 6` @4T: **`soe.factor` drops from 44 calls to 1**, and `soe.cgs` picks up the other 43.
Solver time 4298 → 1511 ms (**2.84×** less solver work); step 13645 → 10840 (**1.26×**), which
matches P1e's 1.22–1.30× at this size.

Read the `%fac` column naively and you conclude factorization **collapsed to 5.0%** of the solver.
It did not. **77.4% is in `soe.cgs`, and any CGS give-up refactorizes *inside* that phase-23 call**
(Intel's automatic fallback), so its cost bills to `soe.cgs`, never to `soe.factor`. On this run
CGS won all 43 (consistent with P1e's finding that the fallback never fired in 340 solves), so the
`soe.cgs` figure here is genuinely iteration and not hidden factorization — but **the parser cannot
tell you that and neither can the profile.** Use `-stats` / `iparm[19]` for the win rate.

This is why `soe.cgs` is a separate bracket and must never be summed into `soe.factor`.

## 4. The assembly side — and a loop closed

At `-matrixType 2` @4T the step decomposes as **formTangent 62.5% / solver 21.8% /
formUnbalance 9.6% / update 5.6%**.

ADR-75b's L3-0 gate measured the element fraction as **74.9%** under `-matrixType 2` @4T — but it
derived that as the *complement* of a solve it could not resolve, because PARDISO had no brackets.
The element-side phases here sum to **77.7%**, and the solve is now resolved rather than inferred.
**L3-0's number holds.**

⚠ **But do not carry that 77.7% forward as "assembly dominates".** It is a property of an 11.5k-DOF
model, and it **inverts with size**: G-L3 measured the largest element loop at **0.95% against a
98.5% solve at 540,675 DOF** (ADR-75b §13, which is why Lane 3 is closed). At production scale
essentially the whole step is the phases this document splits — so the split below matters *more*
there, not less, and the §5 ceiling is a small-model floor.

`soe.addA` (DEEP-gated, n=1): **875.4 ms unsymmetric vs 511.8 ms symmetric**, against a ~7.5 s
`formTangent` — i.e. the scatter is **11.5% / 7.0% of assembly**. Half-storage nearly halves the
scatter, which is the same mechanism as its factor win (half the entries) and a second-order
benefit of `-matrixType 2` that the P1d study did not separate out.

`setSize` is **not** a cost at this size: `dc.s.fill` 16–20 ms unsym / 5.9 ms sym, `dc.s.verify`
0.3–0.5 ms, both **one call** ≈ 0.1% of step. Worth restating that these are **once per sparsity
pattern**, so a workload that re-derives the pattern (remeshing, ADR-51 element removal, contact
re-emission) pays them per `domainChanged`, not per solve — that is where they could matter.

## 5. Method notes / caveats

- **Coarse, not deep, for the headline.** The deep tax is measured, not argued: `-deep` costs
  **+4.5%** on the step (unsym 13645 → 14257) and **+2.8%** (sym 11451 → 11772). That tax lands on
  the assembly side, i.e. on the very fraction being reported, so the headline rows are coarse and
  only the `soe.addA` rows are deep.
- **Threaded PARDISO is not byte-reproducible** (P1f), and this box swings ±30% on background load.
  Median of 3, interleaved; the sweep **aborts** if any compiler is running.
- **The attribution gap is ~0 here** (−0.3 to −0.6 ms) because `LoadControl` solves once per
  iteration. Under `DisplacementControl`/`ArcLength` it is not — the integrator solves in
  `update()` and `newStep()` too, and solver total then *exceeds* `linearSolve`. See the
  `LEDGER_quirks` row; do not generalize this row's zero gap.
- **The 16%-of-step ceiling is a small-model FLOOR, not a ceiling for production.** ADR-75's
  standing correction applies and G-L3 quantifies it: the solve is 21.8% of step here at 11.5k DOF
  but **98.5% at 540,675 DOF**. Factorization is superlinear in N while triangular solve is
  ~linear in nnz(L), so **%fac should also rise with N** — both effects push the same way, making
  factorization-reuse levers *more* valuable at production scale, not less. **Measuring that trend
  (51k / 136k / 540k DOF) is the single obvious follow-up and is NOT done here** — every number in
  §1–§4 is one model at one size.
- **One box, one model class** (plastic `LadrunoBrick` + `LadrunoJ2`), `Newton` + `LoadControl`.

## 6. Bug found while measuring — profiler run attributes (FIXED, P1i)

The HDF5 run attributes were wrong, confirming the open item in ADR-75 §9. Measured on a
`-perStep` probe at `MKL_NUM_THREADS=4` on a 3375-element model:

| attribute | was | actual | verdict |
|---|---|---|---|
| `threads` | `1` | 4 | **bug — fixed** |
| `nElem` | `0` | 3375 | **bug — fixed** |
| `nNode` | `0` | 4096 | **bug — fixed** |
| `nnz` | `0` | — | left 0, deliberately (see below) |
| `nSteps` | `0` | 15 | **NOT a bug** — see below |

**`nSteps` is by design.** It reads 15 *with* `-perStep` and 0 without, because it is derived from
the per-step series, exactly like `dt_min`/`dt_max`. An earlier draft of this document listed it as
a defect; that was wrong. A coarse run has no series to count.

**`threads` was the real one.** `Profiler.cpp` set `m.threads = threads_.size()` — the count of
threads registered with the *profiler*, which is 1 for any run whose command layer is
single-threaded, no matter what MKL is doing. A wrong positive integer is worse than no answer,
because it reads as a measurement. Now resolved from `MKL_NUM_THREADS` → `OMP_NUM_THREADS` →
machine width, then overridden by `mkl_get_max_threads()` wherever MKL is compiled in (which is
authoritative: it sees a programmatic `mkl_set_num_threads()` and reports physical cores when
nothing is declared). The registered count remains a floor, so a genuinely threaded run can never
report fewer threads than it registered.

**`nElem`/`nNode` were promised and never delivered** — a comment in `Profiler.cpp` claimed the
command layer populated them; it only ever filled `nDOF`. Now filled at all four `buildMeta()` call
sites (`profiler report` and `profiler checkpoint`, each existing twice — once in the Python ladder
and once in the completely separate Tcl one).

**`nnz` stays 0 on purpose.** `LinearSOE` exposes no size-agnostic nnz accessor; only some concrete
SOEs carry the member. Filling it would mean adding a virtual to an upstream base class. A
normalizer that exists for some solvers and not others is worse than one that never does.

⚠ **The tables in §1–§4 predate the fix**, so they were produced by a binary that recorded
`threads=1`. Every thread count in this document is the pinned environment value, recorded by the
sweep, never read from the profile.
