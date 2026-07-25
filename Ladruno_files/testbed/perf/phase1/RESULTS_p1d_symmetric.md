# ADR-75 P1d — symmetric PARDISO (`-matrixType`): MEASURED

**Verdict: symmetric wins on both axes, and it is the memory lever BLR was supposed to be.**
**1.94–1.96× UmfPack** at 4 threads and **−42% peak memory**, exactly — not approximately.

> **Correction (2026-07-25, second run).** The first draft of this file headlined
> **1.35× vs unsymmetric PARDISO**. A second full sweep puts that at **1.24×**. The
> symmetric rows are reproducible; the *unsymmetric PARDISO* row is the volatile one (it alone
> improved 14.6% between runs while every other row moved ~7%). Use **~1.25× vs unsymmetric
> PARDISO** — reproduced 3 of 4 thread-count/run combinations — and treat 1.35× as the top of the
> range, not the figure. **Versus UmfPack the symmetric result is stable to ±1%** (1.94–1.96× @4T,
> 1.62–1.63× @1T), so that is the number to quote. The memory result is a deterministic counter,
> not a timing, and is unaffected.

Date: 2026-07-25 · worktree `adr-75-session-handoff-29effa` · same binary for every row.

## What was measured

`system Pardiso -matrixType 0|1|2` (aliases `-symmetric` = 2, `-spd` = 1), where the SOE stores the
**upper triangle only** and PARDISO factorizes with `mtype` 11 / 2 / −2 respectively.

Model: Lane B — 15³ `LadrunoBrick` + `LadrunoJ2` (associated flow ⇒ **genuinely symmetric
tangent**), 11,520 DOF, implicit static, 15 `LoadControl` steps, Newton, RCM.
Harness: `laneB_solver_bench.py` (median of 5 **interleaved** rounds) driven by `sweep_p1d.sh`;
memory by `p1d_memory.py`; Tcl path by `p1d_tcl_smoke.tcl`.

## 1. Time

Two independent full sweeps. **Run A** = the P1d code as first written. **Run B** = after the
adversarial-review fixes (which added an `addA` symmetry check and a perturbed-pivot read), so B also
answers "did the new diagnostics cost anything" — they did not.

| MKL threads | solver | median s (A) | median s (B) | vs UmfPack (A / B) | vs unsym PARDISO (A / B) |
|---|---|---|---|---|---|
| 1 | UmfPack | 18.6426 | 17.7098 | 1.00× | — |
| 1 | Pardiso (unsym) | 14.4015 | 13.5920 | 1.29× / 1.30× | 1.00× |
| 1 | Pardiso `-matrixType 2` | 11.4641 | 10.9578 | **1.63× / 1.62×** | 1.26× / 1.24× |
| 1 | Pardiso `-matrixType 1` | 10.3004 | 9.8248 | 1.81× / 1.80× | 1.40× / 1.38× |
| 4 | UmfPack | 17.6229 | 16.2942 | 1.00× | — |
| 4 | Pardiso (unsym) | 12.1545 | 10.3833 | 1.45× / 1.57× | 1.00× |
| 4 | Pardiso `-matrixType 2` | 8.9878 | 8.3946 | **1.96× / 1.94×** | 1.35× / 1.24× |
| 4 | Pardiso `-matrixType 1` | 9.0499 | 8.3632 | 1.95× / 1.95× | 1.34× / 1.24× |

Tip displacement rel err vs UmfPack is **0.0 in every row of both runs**, except a single 2.0e-16
(one ULP) on run A's unsymmetric 4-thread row.

**Which column to quote — this reversed on the second run.** Run A's draft said to read
"vs unsym PARDISO" because it is an interleaved within-run comparison. Two runs show that reasoning
was incomplete: the **unsymmetric PARDISO row is itself the least reproducible measurement in the
set** (10.38 vs 12.15 s @4T, a 14.6% swing, while every other row moved ~7%), so a ratio built on it
inherits that noise. The `-matrixType 2`-vs-UmfPack ratio is stable to **±1%** across both runs at
both thread counts. So:

- **Quote 1.94–1.96× vs UmfPack @4T** (and 1.62–1.63× @1T). Reproducible.
- **Quote ~1.25× vs unsymmetric PARDISO** — the value reproduced in 3 of 4 run/thread combinations.
  1.35× is the top of the range, not the figure.

Neither UmfPack anchor here is the locked P1 baseline (22.711 s); that came from a different session
on a differently loaded machine. Cross-session absolute times are ±20%.

## 2. Memory — the point of the exercise

`system Pardiso -stats` (new in P1d) dumps PARDISO's own counters: `iparm[14]` peak symbolic,
`iparm[15]` permanent, `iparm[16]` peak numeric, total peak = `max(iparm[14], iparm[15]+iparm[16])`.

| | nnz(A) stored | nnz in factors | **TOTAL PEAK** |
|---|---|---|---|
| unsymmetric (mtype 11) | 818,892 | 8,464,734 | **105.60 MB** |
| symmetric (mtype −2) | 415,206 (**−49.3%**) | 4,459,964 (**−47.3%**) | **61.48 MB (−41.8%)** |
| SPD (mtype 2) | 415,206 | 4,459,964 | **56.60 MB (−46.4%)** |

The stored-nnz figure is exact by construction: `(818892 − 11520)/2 + 11520 = 415206`.

> **Correction (adversarial review).** The first draft read these counters immediately after phase 22
> and reported 105.43 / 61.31 / 56.42 MB. Intel documents `iparm[16]` as the peak over numerical
> factorization **and solution**, so reading it before the first phase-33 solve gives a LOWER BOUND —
> wrong for a figure billed as "the number that decides whether a model fits". Now read after the
> first solve. The correction is ~0.17 MB and the **−41.8% ratio is unchanged**, because both
> configurations were mis-measured identically; the absolute numbers are what needed fixing.

Unlike the timings, these are deterministic counters — they reproduce exactly, so no range applies.

**Contrast with P2b's BLR result — this is the important part.** BLR shrank the *stored factors*
21.8% but moved *peak* memory only 4.6%, because peak is dominated by the active frontal/working
space, not by what is stored. Symmetric half-storage shrinks **both**: the fronts themselves are
half-size, so the fit/no-fit number drops ~42%. And it does so **exactly** — no accuracy tolerance,
no Newton degradation, legal on the byte-identical/oracle paths where BLR is forbidden.

Against the P1c capability wall (UmfPack OOM at 86,490 DOF; PARDISO solved 136,080 DOF), a 42% cut
in peak is the single largest memory lever ADR-75 has produced.

## 3. Correctness

- Tip displacement **bit-identical to UmfPack** (rel err 0.0) for every symmetric configuration at
  every thread count — one 2.0e-16 on the *unsymmetric* row, i.e. within a single ULP. Half-storage
  is exact on a symmetric tangent, as it must be.
- Tcl parity (`p1d_tcl_smoke.tcl`, elastic `stdBrick` cantilever, `OpenSees.exe`): `UmfPack`,
  `Pardiso`, `-matrixType 2`, `-matrixType 1`, `-symmetric`, `-spd` all return
  `ux = 1.035149231914e+00` — identical to 12 digits.
- `matrixType 0` is byte-identical to the pre-P1d code path (the half-storage filter is a
  `!halfStore ||` short-circuit; nothing else in the unsymmetric path changed).

## 4. Why the default stays UNSYMMETRIC

Despite winning on both axes here, `-matrixType 0` remains the default:

1. **Half-storage on an unsymmetric tangent silently solves the wrong system.** Only the `col >= row`
   half of each element matrix is read — no averaging, no detection. The fork has genuinely
   unsymmetric tangents: contact (`LadrunoContact`), non-associated flow, follower loads,
   `LadrunoUP`. Those must not be silently symmetrized by a default.
2. ADR-75 §5 quotes the Kratos precedent as *"default symmetric but expose the override"*. That
   holds for a code whose element library is symmetric by construction; this fork's is not.
3. The win is one explicit token away, and `-stats` now makes it self-evidently worth taking.

**Recommendation to model authors:** if your tangent is symmetric (elastic, J2/associated plasticity,
no contact), use `-matrixType 2`. Prefer **2 over 1**: at 4 threads they are within 0.4% on time,
`2` is only ~8% worse on peak memory, and `1` (SPD, Cholesky, no pivoting) **fails outright on an
indefinite tangent** — which any softening or buckling model eventually is. The `-4` zero-pivot error
now says so explicitly. `1`'s single-thread edge (1.38–1.40× vs 1.24–1.26×) is real but not worth the
fragility.

**It now warns if you get this wrong.** `addA` compares each discarded lower-triangle entry against
its mirror during the first few tangent assemblies and reports once if they differ. That converts the
worst failure mode here — a converged, plausible, wrong answer — into a loud message. It is a
first-pass detector, not a guarantee: a tangent that only becomes unsymmetric later (contact closing
at step 50, non-associated flow after first yield) is outside its window.

## 4b. What the adversarial review changed

A scoped adversarial review ran against the first implementation. It **cleared the two areas I was
least sure of** — it re-derived the graph semantics (`AnalysisModel::getDOFGraph` → `addEdgeFast`;
`Vertex::addEdge` drops self-edges; `ID::insert` keeps adjacency sorted and deduplicated, so the
off-diagonal count is always even) and then brute-forced the fill loop over 400 random graphs × 3
matTypes, finding `lastLoc == nnz` exactly every time. It also killed the `getTag() != a` worry:
`Graph::getVertexPtr` is a lookup **by tag**, so `vertexTag == a` inside that loop is a tautology.

What it found instead was concentrated in *undetected wrongness*, and all of it is now fixed:

| # | Finding | Status |
|---|---|---|
| 1 | Perturbed pivots (`iparm[13]`) never read — mtype ±2 runs eps 1e-8, five orders looser than the unsymmetric 1e-13, so a near-singular tangent returns `error == 0` and a solve to a *perturbed* matrix | warns per factorization + in `-stats` |
| 2 | The asymmetry hazard was documented but not detected, though the data is already in hand | detected in `addA` (§4) |
| 3 | The in-flight "row starts at the diagonal" check is false for `matType 0` | already presence-based, not position-based |
| 4 | A bounds *clamp* would turn overflow into silently dropped entries; `lastLoc != nnz` is exact and stronger | already `return -1` + the exact equality |
| 5 | The OOM path was dead code (throwing `new`), `Asize` clobbered the recovery, `delete[]` before allocation left danglers | `std::nothrow` + `return -1` |
| 6 | The Tcl parse loop silently swallowed unknown options and a valueless `-matrixType` — the same class of failure this feature exists to kill | both warn |
| 7 | The CMake comment's "mixing MKL layers" rationale was a non-sequitur; SP/MP already link the PARDISO objects. And `system Pardiso` there fell through to a **stale global `theSOE`** | rationale corrected; explicit refusal added |
| 8 | `-stats` read one phase too early | fixed (§2) |
| 9 | The `iparm[10]/[12]` comment misstated Intel's position — Intel *supports* scaling+matching for mtype −2 and recommends it for saddle-point systems | comment corrected; no code change |
| 10 | `n <= 0` guard, uninitialized `idum`/`ddum`, `setLinearSOE` not resetting per-SOE state, and a dead second symmetric implementation in the same directory | all fixed; the dead pair is now a `LEDGER_quirks` entry |

## 5. Reproduce

```bash
cd Ladruno_files/testbed/perf/phase1
./sweep_p1d.sh                                        # time, 1 + 4 threads
OPS_PYD=<worktree>/dist/bin python3.12 -S p1d_memory.py   # peak memory, -stats
<worktree>/dist/bin/OpenSees.exe p1d_tcl_smoke.tcl        # Tcl parity
```

## 6. Two traps banked while measuring this

- **`-matrixType` as a STRING silently costs you the solver.** openseespy
  `system('Pardiso','-matrixType','2')` fails `OPS_GetIntInput`; the factory returned 0, `theSOE`
  stayed null, and OpenSees fell back to the **default `ProfileSPDLinSOE`** with only a warning. On a
  3D solid mesh that is catastrophically slow — a Lane-B bench ran 25 minutes before being killed,
  looking merely "slow" rather than wrong. Fixed twice over: the parse failure now degrades to
  unsymmetric instead of returning 0 (which is what the warning always claimed), and the sweep script
  keeps the full log instead of grepping four lines out of it.
- **MKL 2026.1 renamed its runtime DLLs `.2.dll` → `.3.dll`**, so `build.bat`'s hardcoded staging list
  copied nothing and `import opensees` died with "DLL load failed". Written up as
  `BUILD_GOTCHAS.md §9` (fixed concurrently by #627 — base-name globbing + a stale-DLL purge; P1d only adds the new `mkl_avx10` kernel).
