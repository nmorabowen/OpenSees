---
title: WP-253 — MPM performance baseline and the explicit transfer-loop speedup (#1097)
project: Ladruno
tags:
  - mpm
  - performance
  - wp-253
---

# WP-253 — MPM performance baseline and the explicit transfer-loop speedup

Program #1075. Queued by the owner on 2026-11-23, ahead of 232m (phase 7, contact): contact and
the FE–MPM seam both add work to mode F's per-step transfer loop, which M0b measured at 91 % of an
elastic step at 10⁵ on one thread, with threads changing nothing (#1097). This WP measures, then
makes that loop fast. It is not phase 11: MPI is out of scope, and so are thread-safe soil models,
#1096 and #1098 (reported with numbers if they dominate; the owner decides).

**Status:** pre-registration (this section and the next four are committed before any number is
taken). The #1097 target is committed after the profile and before any code change, in
"The #1097 target". Results follow in "Results"; nothing above "Results" moves once a number
exists. A miss goes to the owner.

## 1. What is measured — the benchmark matrix

Two models, both after M0b (`Ladruno_scripts/mpm_m0b/models.py`), so the numbers continue M2's:

- **column (mode F)** — the granular column collapse on a frictional floor (`mpm wall zmin`): 3D
  is M2's quarter column (a × a × 2a, symmetry `fix xmin x`, `fix ymin y`), 2D plane strain is
  its half-column analogue (a × 2a, `fix xmin x`). ppc 2 per axis. Δt = 0.25 Δt_crit as M2.
- **block (mode A)** — M1's laterally confined block under self-weight, quasi-static, PARDISO,
  `renormalize` 5e-2; its 2D plane-strain analogue.

| Axis | Values |
|---|---|
> **Amendment (2026-11-23, before any timed number).** The two rules above conflict in 3D: M2's
> a × a × 2a quarter column at ppc 2 gives N = 8,192 or 11,664 near 10⁴, neither within 10 %.
> The harness keeps height = 2 × width (x) and lets the y extent differ, giving N = 9,800,
> 10,368, 98,496 and 998,400 (`Ladruno_scripts/mpm_bench/models.py`). The only numbers taken
> before this note are the harness's smoke runs (1 repeat, N ≤ 10⁴, a loaded machine), which are
> not results and set nothing.

| mode | F (explicit), A (implicit, quasi-static) |
| basis | `gimp`, `bspline2` |
| material | `ElasticIsotropic`; `DruckerPrager` (UW, as M2's column / M1's block); `LadrunoSANISAND` (IntScheme 2, the fork's test deck parameters, kPa deck) |
| dimension | 2D plane strain, 3D |
| N (particles) | 10⁴, 10⁵ (nmb); 10⁶ mode F only (esmeralda). The model's cell count is chosen so N is within 10 % of the nominal; the actual N is recorded |
| threads | 1, 4, 8 (nmb); 1, 4, 8, 24 (esmeralda) — `ladrunoThreads <n>`, with `OMP_NUM_THREADS` = `MKL_NUM_THREADS` = n |

Exclusions fixed now, with their reason:

- `LadrunoSANISAND` in mode F only on `bspline2` (232k's K0 gimp miss, owner's decision); its
  material update is never threaded (W14). At 1, 4, 8 threads its rows measure the transfer loops
  only (see §4).
- Mode A: N = 10⁴ for every cell of the matrix; N = 10⁵ only for 3D `gimp` elastic and DP (one
  mode-A run at 10⁵ is 30–50 min). Threads 1 and 8 only (mode A's threading is the solver's and
  WP-192's, not this WP's). SANISAND in mode A: 10⁴, 3D, `gimp`, 1 thread.
- Mode F at 10⁶: esmeralda only, 3D and 2D, every basis/material cell above, SANISAND excepted
  (its plastic call is ~0.6 ms; one step at 10⁶ would be minutes).

## 2. Metrics

- **Throughput T** [particle-steps/s] = N × n_timed / W, with W the sum over the timed steps of
  the profiler's `step` scope wall (coarse profile run), as M2. ns per particle-step = 10⁹ / T.
  Mode A also records the Newton iterations per step k and T per Newton iteration.
- **Phase split** (mode F): the share of W in each named child scope of the step. The scopes are
  added by this WP's instrumentation commit (profiler scopes only; they are no-ops with the
  profiler off and change no result — checked bit for bit before the baseline):
  `mpm.f.dtcrit` and `mpm.f.admissible` (host, before the advance), and under `mpm.advance`:
  `mpm.f.support` (support search and shape functions), `mpm.f.p2g` (mass, momentum and forces,
  image folds), `mpm.f.grid` (damping, grid momentum update, grid conditions), `mpm.f.g2p`
  (particle velocity and position, the trial admissibility check), `mpm.f.project2` (the second
  projection and its grid conditions), `mpm.f.update` (ΔF, the kernel, the constitutive call, the
  stress push), `mpm.f.commit` (material commit); then `commit` (Domain::commit, recorders). The
  remainder of `step` is reported as `other`. The material's own share inside `mpm.f.update` and
  `mpm.f.commit` comes from M2's deep run (`mpm.f.material` elem_by_type); deep and coarse runs
  are separate processes.
- **Speedup** S(n) = T(n threads) / T(1 thread), same build; **gain** G(n) = T_after(n) /
  T_baseline(1), the after build against the baseline build.
- **Determinism**: every particle array the engine keeps (x, v, V, F, R, eps, sigHat, sig) after
  K steps, compared at n ∈ {2, 4, 8} threads against 1 thread with exact float64 equality. Any
  tolerance instead would be pre-registered here, and none is.

## 3. Machines and builds

- **nmb**: Intel i7-13700KF (8 P-cores with HT + 8 E-cores), 64 GB, Windows 11 IoT Enterprise
  LTSC 2024; MSVC through `Ladruno_scripts\build.bat` (all five targets), Intel oneAPI MKL
  2026.1, Python 3.12. Runs pin the process to the P-cores (affinity `0xFFFF`, as M0b).
- **esmeralda**: a 24-core Slurm node (model recorded from `lscpu`), GCC and spack MKL as the
  fork's Linux guide; a fresh clone for this WP (never `~/ladruno_build_test/OpenSees`), every
  job `--exclude=node1,node2,node12,node13 --mem-per-cpu=…`, `--exclusive`, `OMP_PROC_BIND=close`
  `OMP_PLACES=cores`.
- **Builds**: the *baseline* is `origin/ladruno` `a3820d5ac` plus the instrumentation commit
  only; the *after* build is the branch head at the ready flip. Every result row records
  `ladrunoBuild` from both `OpenSees.exe` and `opensees.pyd` (they must agree with each other and
  with the commit measured), the machine, the MKL version and the thread count.

## 4. Noise protocol

- Every run is a fresh process; 3 repeats per row; the row's value is the **median**, its spread
  (max − min) / median. A row with spread > 10 % is re-run once with 5 repeats and the 5-repeat
  median replaces it (both kept in the data).
- Mode F: 5 warm-up steps (excluded), 20 timed steps in one `mpm analyze` under one profiler run.
  Mode A: 2 warm-up + 6 timed, as M1.
- Two medians differ only if their ratio differs from 1 by more than the sum of their spreads;
  otherwise the difference is reported as "within noise".
- No build or other heavy job runs during a timed run on nmb (the orchestrator sequences them).
- A timeout or crash is data (`ok: false`, reason, wall), as M0b. A crash at 8 threads gets a
  minidump and the build SHA posted on #1095.

## 5. Method for the speedup (fixed before the profile; the target is not)

- **The 1-thread path is the existing serial code, unchanged** (the WP-107 policy: at 1 thread
  the byte-identical serial path, not a one-thread parallel region). Every existing MPM record
  therefore stays bit for bit at 1 thread.
- **Threads come from the existing knob** (`ladrunoThreads`, `ladruno_loopThreads()`), run
  through `ladruno_parallelFor` (OPS_Utilities); the core stays OpenSees-free behind a hook.
- **Bit-identity at n threads, by construction:** particle-local loops (supports, G2P, ΔF, the
  kernel, the stress push) are independent per particle. Grid accumulations (P2G, the second
  projection) are computed as a node-centric *gather* over each node's contributors in ascending
  particle order — the serial scatter's exact summation order — so no reduction order changes.
  Grid-node loops are independent per node. The serial-order bookkeeping (the touched list, the
  image folds) is kept in serial order.
- **Materials:** the material calls (`setTrialStrain`, `commitState`) run threaded only when
  every particle's material is on the WP-107 allowlist (`ladrunoThreadSafeUpdate()`), with
  `ops_TheActiveElement = 0` on every worker thread (W9); otherwise they run serially in particle
  order, so a non-allowlisted material sees exactly the serial call sequence. The transfer loops
  thread in both cases — they touch no material. This is the reading of "soil models stay serial"
  used here; it is flagged in the handoff.
- **Failures keep the serial semantics:** the lowest failing particle is named and the materials
  touched are reverted, as the serial loop does.
- Route S (`advanceUP`, two-phase) and mode A are not threaded by this WP; mode A's re-pointing
  is threaded only if the mode-A profile shows it above 10 % of a step.

## The #1097 target

*Committed after the baseline profile, before any code change.*

## Results

*Empty until the baseline runs.*
