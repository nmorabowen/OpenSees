---
title: "WP-194 — mixed-precision MUMPS (single-precision factor, double residual)"
project: Ladruno
type: performance work package (solve lane)
status: "prototype, PR open (draft) — not merged; default off; use gated on job 948"
priority: medium
owner: nmora
related:
  - "[[40_ladruno_performance_adr]] (lever 1: sparse direct solvers)"
  - "[[188_ladruno_factor_reuse_adr]] (factor reuse: the kept factor is the float one here)"
  - "WP-190 (PR #978, MUMPS -stats phase timings: same labels in the single-precision block)"
  - "paper repo: review/engstruct-2026/30_campaigns/timehistory/perf/GPU_MIXED_SCOPE.md (§3 SMUMPS bench, §5b the 948 gate)"
tags: [performance, solver, mumps, mixed-precision, mpi, wp-194]
updated: 2026-10-14
---

# WP-194 — mixed-precision MUMPS

> [!summary] The short version
> ```tcl
> system Mumps ... -precision single ?-refine k?      ;# or the alias -mixed
> ```
> factors the distributed tangent with **SMUMPS** (single precision) instead of DMUMPS.
> Assembly, the residual, the solve's right-hand side and solution vectors, and the Newton
> update stay double; only the LU is float. Off by default (byte-identical when off).
> Prototype built so it is ready if the real-4D gate (job 948) says go; **do not use it on a
> production run before that gate passes** (§6).

## 1. Why

On the 4D-sized proxy (n = 139 968, 3.04e8 factor entries) SMUMPS factored **2.1× faster**
than DMUMPS on gpu01 (2×8: 1.15 vs 2.42 s) with **about 0.6× the memory**. Its first solve
has a relative residual of 2e-2; each double-precision refinement step contracts it by about
1e-3 (2e-2 → 5e-5 → 7e-8 → 2e-11). Paper repo, `perf/GPU_MIXED_SCOPE.md` §3; driver
`perf/gpu_mixed/mumpsx.c`.

Newton already supplies that double-precision correction: every iteration evaluates the
exact residual in double and solves for a correction, so a float factor only turns Newton
into an inexact Newton (linear convergence at roughly cond·ε_single per iteration instead
of quadratic), exactly as BLR and WP-188 reuse do. With a 2.1× cheaper factor and at most
one extra iteration, the 4D projection is 56 h → 39–45 h per record on `computes`
(GPU_MIXED_SCOPE §4).

## 2. What

- `-precision single` (alias `-mixed`; `-precision double` = the default) on `system Mumps`,
  in both parsers: classic Tcl `SRC/tcl/commands.cpp` (OpenSeesSP.exe, OpenSeesMP.exe) and
  `OPS_MumpsSolver` in `SRC/interpreter/OpenSeesCommands.cpp` (openseespy MP).
- `-refine k` (default 0): k classical mixed-precision refinement steps inside each solve —
  r = b − A x in double against the assembled double matrix (distributed: x broadcast,
  local triplet matvec, `MPI_Reduce` to rank 0), correction through the float factor. Default
  0 because Newton refines; k ≥ 1 is for a linear (or modified-Newton) path that wants a
  double-accurate solve from a float factor.
- Coexists with `-ICNTL7`, `-ICNTL14`, `-BLR`/`-ICNTL35`/`-CNTL7` (all passed to the SMUMPS
  instance; BLR then compresses the float factor) and `-stats`.
- `-stats` in single mode prints the same block as double, with the header line unchanged (the
  existing log parsers match it with anchored regexes) and one extra line
  `  precision           = single (SMUMPS factor, double residual), refine k`, the WP-190
  timing lines (`t_analysis`, `t_factor`, `t_solve`) and, on the first solve after each
  factorization, the residual history `  mixed: ||b-Ax||inf/||b||inf  step i = …`. The double
  path's block is untouched.
- WP-188 `-reuseFactor`: works unchanged. Reuse keeps the SOE `factored` flag set, so the
  kept factor is the float one and a reused step runs SMUMPS JOB=3 only.
- `sendSelf`/`recvSelf` carry `precision` and `refine` in their own `ID(2)` after the existing
  controls, so OpenSeesSP subordinates instantiate SMUMPS too (a DMUMPS/SMUMPS pair across
  ranks would deadlock or mis-solve).
- Fail closed: a bad `-precision` value, a missing value, `-refine` without single,
  `-precision single` on the serial MUMPS path, or on a build without SMUMPS **refuse the
  command** (`Mumps ERROR …`, TCL_ERROR / null SOE); nothing silently factors in double.

## 3. How (code)

| File | Change |
|---|---|
| `SRC/system_of_eqn/linearSOE/mumps/LadrunoMumpsMixed.{h,cpp}` (new) | The SMUMPS instance (`SMUMPS_STRUC_C`), float copy of the local triplet values, analysis / factor / solve, RHS scaling by ‖b‖∞ (keeps the float RHS clear of overflow and underflow), `-refine`, the stats block. Whole file under `#ifdef _LADRUNO_SMUMPS`. |
| `MumpsParallelSolver.{h,cpp}` | `precision`, `refineSteps`, `LadrunoMumpsMixed *theMixed`; `setPrecision()`, static `singlePrecisionAvailable()`, `solveMixed()`. `solve()` diverts to `solveMixed()` before `initializeMumps()` when single, so the DMUMPS instance is never created (that is where the memory saving comes from). Teardown and `setCommunicator()` also terminate the SMUMPS instance. |
| `SRC/tcl/commands.cpp`, `SRC/interpreter/OpenSeesCommands.cpp` | The options and refusals. |
| `CMakeLists.txt` | `LADRUNO_MUMPS_SINGLE` = AUTO (default) / ON / OFF. AUTO compiles the option in when `${MUMPS_DIR}/smumps.lib` (`libsmumps.a`) **and** `smumps_c.h` exist; then the library is prepended to `MUMPS_LIBRARIES` and `_LADRUNO_SMUMPS` appended to `MUMPS_FLAG` (both already reach exactly the SP/MP/PyMP targets). ON makes their absence a configure error. `LadrunoMumpsMixed.cpp` is in the three targets' source lists. |

Design choices:
- **No DMUMPS factor at all in single mode.** Keeping a double factor "just in case" would
  cancel the memory saving. The double matrix the SOE already holds is the only double copy
  (the refinement matvec reads it); the float copy of the values is 0.5× nnz·8 bytes per rank.
- **Index shifting.** The double path shifts `rowA`/`colA` to 1-based around every MUMPS
  call and back; the mixed path does the same for JOB 1 and 2, and on every error path too.
- **Overflow guard.** An entry beyond `FLT_MAX` (or NaN) refuses the factorization
  collectively (`MPI_Allreduce`) with a message, rather than factor `inf`.

## 4. Tests

`tests/test_wp194_mumps_mixed.py` (zone_a; OpenSeesMP 2 ranks) on
`tests/wp194_mumps_mixed/mixed_frame.tcl`, which runs the WP-188 frame deck unchanged
(`frame_gm.tcl`: 3-story RC fiber frame, Concrete01/Steel02, KrylovNewton, NormDispIncr 1e-8)
and only wraps `system` (appends `$WP194_MUMPSARGS` to `system Mumps`) and `analyze` (sums
`testIter`).

| Test | Pins |
|---|---|
| default off | no option vs `-precision double`: same iterations, roof history to 1e-12 relative (exact equality is not testable on 2 ranks, see §7) |
| refusals | `-precision quad`, `-refine 2`, `-precision double -refine 1`, `-refine -1` refuse; a build without SMUMPS refuses `-precision single` |
| single vs double (`-precision single`, `-mixed`) | same steps, no failure, roof history within 1e-5·peak, iterations reported |
| `-refine 2 -stats` | header line unchanged, precision line, residual history decreasing |
| `-BLR 1e-5` | converges, BLR header |
| `-reuseFactor 10` | converges, fewer factorizations than steps |

Results (local, Windows, i7-13700KF P-cores, 2 ranks, 1 thread/rank, MUMPS 5.9.0 + METIS
with `BUILD_SINGLE=ON`): see §7.

## 5. Build

- **Windows (local):** the MUMPS 5.9.0 + METIS toolchain `C:\Users\nmb\toolchains\mumps-5.9.0-metis\`
  had only the `d` arithmetic. WP-194 built `s`+`d` into a separate prefix
  `C:\Users\nmb\toolchains\mumps-5.9.0-metis-sd\` with the paper repo's recipe
  (`pins/windows_mumps_metis/build_mumps.bat`) changed in three places only: `-DBUILD_SINGLE=ON`,
  separate build/stage/prefix directories, and `smumps.lib` added to the `mumps5_*` METIS rename
  and its verify (VERIFY_OK: `smumps.lib` references `mumps5_METIS_SETDEFAULTOPTIONS` only, like
  `dmumps.lib`). The shared d-only prefix is untouched.
- **Header trap:** the fork puts `OTHER/MUMPS/include` (MUMPS 5.5.1 headers) on every TU's
  include path ahead of `MUMPS_INCLUDE_DIR`. A 5.9.0 build therefore force-includes the
  5.9.0 `dmumps_c.h` (`/FI`, as the v3+METIS pin does), which also pins the 5.9.0
  `mumps_c_types.h`. Without it, a TU that includes both arithmetics would get the 5.5.1
  `mumps_c_types.h` first (through the in-tree `dmumps_c.h`) and its include guard would
  hide the 5.9.0 one from `smumps_c.h`. Today that is harmless (the two differ only by
  `MUMPS_OFF_T`, which neither `smumps_c.h` nor `dmumps_c.h` uses), but it is luck, so
  WP-194's local build force-includes both 5.9.0 headers.
- **esmeralda:** `perf/gpu_mixed/opt/lib` already holds an SMUMPS built for the bench, with the
  `mumps5_*` METIS rename applied to a private copy. A fork build there needs `libsmumps.a`
  in `MUMPS_DIR` next to `libdmumps.a`; CMake then picks it up (AUTO). Not done in this WP
  (no esmeralda builds, no pins).
- **CI (Zone-A, Linux):** no SMUMPS → AUTO compiles the option out; the tests then check the
  refusal and the default-off identity and skip the single-precision cases.

## 6. When to use it — the 948 gate

The proxy is a mild uniform solid (RINFOG(10) ≈ 3.5e3). The real 4D tangent couples shells
and soil under a cracking, non-associated tangent, with Penalty 1e14 constraints, and its
condition number is unknown. The gate is set in GPU_MIXED_SCOPE §5b, before the results:

- **GO** if SMUMPS (16×1) on the t17 tangent analyses and factors without error **and** its
  FP64 refinement contracts the relative residual by **≥ 1e-2 per step**, reaching ≤ 1e-8
  within 3 steps;
- **NO-GO** if SMUMPS fails, contracts by less than 1e-2 per step, or RINFOG(10) × 6e-8 ≳ 1e-1.

Even on GO, a production switch needs one more measurement this WP cannot make locally:
a short 4D window (a few hundred strong-phase steps) with `-precision single` against the
double run at the same NormDispIncr 1e-4 — iterations per step (budget: at most +1 on 3.5),
wall time per step, `INFOG(21)`, and the EDP drift against double. Until then it stays an
opt-in prototype.

## 7. Results

Local, Windows, i7-13700KF, 2 ranks pinned to P-cores 12/14, 1 thread per rank, MUMPS
5.9.0 + METIS (`s`+`d`), OpenSeesMP built from this branch with `-j 4` on P-cores.

**Tests.** `tests/test_wp194_mumps_mixed.py`: 10 passed, 1 skipped (the no-SMUMPS refusal)
on the SMUMPS build; on a `LADRUNO_MUMPS_SINGLE=OFF` build of the same tree: 6 passed (default
off, the four refusals, `-precision single` refused with `Mumps ERROR: -precision single needs
single-precision MUMPS (libsmumps)…`), 5 skipped.

**Default off.** The stock 2-rank MUMPS path is not byte-reproducible run to run (quirk
`WP-194-mumps-two-rank-runs-are-not-byte-reproducible`), so it is compared as distributions,
10 runs per arm, roof-history md5 (frame deck, 300 steps):

| Model | pre-WP-194 binary | WP-194, no option | WP-194, `-precision double` |
|---|---|---|---|
| elastic | modal md5 9/10, 820 iterations | same modal md5 9/10, 820 | same modal md5 10/10, 820 |
| fiber | modal md5 4/10, 1287 iterations | same modal md5 4/10, 1287 | same modal md5 3/10, 1287 |

The pre-WP-194 binary is the same tree with the four touched sources reverted to
`origin/ladruno` and `LADRUNO_MUMPS_SINGLE=OFF`.

**Single vs double, nonlinear.**

| Model | double: iterations | single: iterations | max \|Δroof\| / peak | notes |
|---|---|---|---|---|
| RC fiber frame, 300 steps, NormDispIncr 1e-8, KrylovNewton (n = 27) | 1287 | 1287 (+0 %) | 8.5e-13 | `-mixed` alias identical |
| same, `-reuseFactor 10` | 1371 | 1371 | 1.8e-12 | 160 float factorizations for 300 steps |
| 3D J2 block (`perf_block.tcl`, 20×20×10 H8, about 14 k equations), 40 steps, NormDispIncr 1e-6 | 155 | 155 (+0 %) | 4.8e-11 | `-refine 1`: 4.2e-16; factor MB/proc INFOG(21) 82 → 48 (0.59×) |
| same, `-reuseFactor 10` | 186 | 186 | 1.1e-6 against full-Newton double, the same as double with reuse | the reuse (WP-188) tolerance-scale shift, not the precision |

`-refine 2 -stats` on the elastic frame: ‖b−Ax‖∞/‖b‖∞ 3e-7 → 7e-14 → 6e-16 per solve.

Reading: on these well-conditioned small models a float factor costs **no** extra Newton
iteration (cond·ε_single ≈ 1e-7·cond ≪ the per-iteration contraction Newton already gets).
That says nothing about the 4D tangent, whose conditioning is exactly what job 948 measures.
Wall time on this shared desktop is noise-dominated (the same double arm ran in 18.6 s and
38.6 s); alternating pairs on the block gave single 31.9/34.0 s vs double 38.6/38.1 s, but
the per-factorization timing that matters is the gpu01 bench (2.1×, GPU_MIXED_SCOPE §3).

**Link note (Windows).** SMUMPS (and MKL ScaLAPACK's single-precision PBLAS) call `SAXPY`/
`SDOT`; MKL defines them in objects that also define the C names `saxpy`/`sdot`, which the
fork already defines (sparseSYM `utility.c`, in-tree METIS 4 `util.c`), so the first link
failed with LNK2005. `LadrunoMumpsMixed.cpp` defines reference `SAXPY`/`SDOT` (BLAS-1 only;
the factorization kernels stay MKL's) so those MKL objects are never pulled. A Linux build
with static MKL may need the same with `saxpy_`/`sdot_`; untested.
