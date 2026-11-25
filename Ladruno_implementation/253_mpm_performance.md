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
  > **Amendment (2026-11-24, the critic of record, finding 5; wording only, no target moves).**
  > The threaded touched list is in *discovery* order, which depends on thread timing; only the
  > image folds are kept in the serial first-touch order (re-sorted by each node's lowest entry).
  > Bit-identity holds because every walker of the touched list is per node; a guard comment at
  > the list's construction says so.
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

*Committed after the baseline profile, before any code change (2026-11-23).* Baseline: build
`0efa8ca17` (origin `a3820d5ac` + scopes, bit-identical to origin with the profiler off, coarse
and deep), nmb, 1 thread, 3 repeats, `Ladruno_scripts/mpm_bench/results/baseline_F_t1/`.

**What the profile says (3D `gimp` elastic, N = 98,496, 725 ns per particle-step, spread 6 %;
M2 measured 722):** `mpm.f.update` 41.1 % (of which the material ≈ 9.4 % of the step — the rest is
ΔF, the polar decomposition and the stress push), `mpm.f.p2g` 21.1 %, `mpm.f.support` 15.8 %,
`mpm.f.g2p` 10.5 %, `mpm.f.project2` 5.7 %, `mpm.f.commit` 2.2 %, `mpm.f.dtcrit` 2.0 %,
`mpm.f.grid` 0.6 %, `mpm.f.admissible` 0.5 %, other 0.4 %. The dominant loop is the particle
update's kinematics (≈ 32 % of the step), then P2G and the support search. **No sub-loop is
inherently serial**: all are per particle or per node. What stays serial under §5 is the touched
list in first-touch order, the image folds and the gather's transpose, estimated at ≤ 5 % of
the step. DP: the material is 49.6 % of the step (3D gimp 10⁵); 2D elastic: 14.1 %
(`ElasticIsotropicPlaneStrain2D` is not `ladrunoThreadSafeForm`, so its calls stay serial).

**Amdahl bounds at 8 threads** (serial fraction s): 3D elastic s ≈ 0.05 → 5.9×; 2D elastic
s ≈ 0.19 (material + residue) → 3.4×; DP s ≈ 0.55 → 1.6×; SANISAND s ≈ 0.9 → 1.1×. The
transfer loops are memory-bound and the i7's 8 P-cores share one memory controller, so the
targets sit below the bounds:

| # | Row (mode F, nmb unless stated) | Target |
|---|---|---|
| T1 (#1097 headline) | 3D `gimp` elastic, 10⁵, 8 threads | **S(8) ≥ 4.0** (T ≥ 5.5 M particle-steps/s; non-material ns per particle-step from ~657 to ≤ 165) |
| T2 | 3D `gimp` and `bspline2` elastic, 10⁵, 4 threads | S(4) ≥ 2.5 |
| T3 | 3D `bspline2` elastic, 10⁵, 8 threads | S(8) ≥ 4.0 |
| T4 | 2D `gimp` and `bspline2` elastic, 10⁵, 8 threads | S(8) ≥ 2.5 |
| T5 | DP, 2D and 3D, both bases, 10⁵, 8 threads | S(8) ≥ 1.4 |
| T6 | SANISAND rows | reported only; S(n) ≥ 0.95 (threads never slow a run down) |
| T7 | every row at 1 thread, after vs baseline | G(1) within noise (§4) of 1.0, and bit-identical results |
| T8 | every row, every n | bit-identical to 1 thread (§2) |
| T9 | esmeralda, 3D `gimp` elastic, 10⁶, 24 threads | S(24) ≥ 8.0 |

T1–T5 and T9 are read on the medians, with §4's noise rule: a median under its target by less
than its spread is reported as "within noise of the target", not as a pass. A miss goes to the
owner.

## Results

**Binary of record:** `3d2cb6330` (exe and pyd `ladrunoBuild` checked before every run). Two
later commits change code. `f9964c50e` changes only GCC/Clang compile flags (`-ffp-contract=off`
appended; MSVC binaries are unchanged). `dc29414ae` is the review fix to the mode-F transpose.
Its per-block histogram was a dense nb x nTouched array, so memory and the per-thread clear/sum
work did not shrink with threads. Each block's row now spans only the slots of the blocks that
own the nodes it touches. Slots, `tPtr` and `tEntry` are unchanged, and incoherent particle order
degrades to the old size, never beyond it. The node claim became a compare-exchange: a losing
`exchange` overwrote the stored owner id, a heap corruption at 8 threads that the shuffled-order
case and a temporary self-check caught. See the quirk ledger. Confirmation on `dc29414ae`: the
full MPM suite passes (536, now including the shuffled-order bit-identity case), and the 3D
elastic 1e5 rows at 1 and 8 threads are within noise of `3d2cb6330`:
gimp 1.41/6.53 M vs 1.42/6.24 M; bspline2 1.04/4.44 M vs 1.06/4.80 M, with ratio 0.925 inside a
noise band of 0.101 (`results/after_F_dc29414ae_transpose/`). The tables below stay on
`3d2cb6330`. T9 (esmeralda job 10272) was queued on that build. Machine nmb: i7-13700KF (8 P + 8 E cores, 24 logical), Windows 11,
MKL 2026.1.0.0, affinity 0xFFFF. Data:
`Ladruno_scripts/mpm_bench/results/after_F_3d2cb6330/` (+ `_rerun`, the 5-repeat reruns of §4),
`after_A_3d2cb6330/` (+ `_rerun`), baselines `baseline_F_t1/` and `baseline_A_1e4/` (build
`0efa8ca17`). Each directory has a `summary.md`.

### Target verdicts

| # | Measured (median, spread) | Verdict |
|---|---|---|
| T1 | 3D `gimp` elastic 10⁵: S(8) = **5.00** (T = 7.18 M p-steps/s, 3 %); non-material ns per particle-step 657 → ≈ 127 (1-thread share applied to the 8-thread time) | **pass** |
| T2 | S(4) = 2.97 (`gimp`), 2.83 (`bspline2`) | **pass** |
| T3 | 3D `bspline2` elastic 10⁵: S(8) = 4.54 (4 %) | **pass** |
| T4 | 2D elastic 10⁵: S(8) = 3.50 (`gimp`), 3.55 (`bspline2`) | **pass** (10⁴, not a target row: 2.53 and 2.78, spreads 30 % and 23 % after the 5-repeat rerun) |
| T5 | 3D DP 10⁵: 1.54 (`gimp`), 1.71 (`bspline2`); **2D DP 10⁵: 1.27 (`gimp`), 1.31 (`bspline2`)** | 3D pass; **2D MISS, owner-accepted 2026-11-25** |
| T6 | SANISAND S(8) 1.02–1.11 (`bspline2`, 2D and 3D, 10⁴ and 10⁵) | **pass** (2D rows noisy, see below) |
| T7 | Mode F, every row at 1 thread: G(1) 0.99–1.06, each "faster" or "within noise"; bit-identical (bitcheck against `a3820d5ac` at 1/2/4/8 threads, and every digest). Mode A, 3D 10⁴: G(1) 1.02, bit-identical to `baseline_A_1e4` | **pass** |
| T8 | Mode F: every 4- and 8-thread digest equals 1 thread (240/240 runs) | **pass** |
| T9 | esmeralda, 3D `gimp` elastic 10⁶, 24 threads | *pending job 10272* |

### Misses for the owner

- **T5, 2D DruckerPrager (S(8) 1.27–1.32 against 1.4).** The target was set from the 3D
  material share (≈ 50 %). In 2D the DP material is 61–65 % of the step at 1 thread
  (`phase split`, `material share`), and those calls stay serial under §5. Amdahl at s ≈ 0.63
  gives 1.45 at 8 threads, so the target sat within 4 % of the bound; the transfer loops
  themselves scale as in 3D. Lifting it needs a thread-safe DP (`ladrunoThreadSafeForm`), which
  is out of scope (thread-safe soil models).
  **Owner-accepted on 2026-11-25** (relayed by the MPM orchestrator, #1075): merge on green with
  the target unchanged. WP-258, the thread-safe soil materials, opts DP, SANISAND and HS into
  the threaded material path.
- **T9** *(filled from esmeralda)*.

### What "soil models stay serial" means in the numbers

Only 3D `ElasticIsotropic` passes `ladrunoThreadSafeUpdate() && ladrunoThreadSafeForm()`; every
other material runs the split path (threaded kinematics, serial material calls in particle
order, threaded push). So the speedup of a row is capped by its material share: elastic 3D
(6–9 %) → S(8) 4.5–5.1; elastic 2D (11–14 %, `ElasticIsotropicPlaneStrain2D` not Form-safe) →
2.5–3.6; DP (42–65 %) → 1.3–1.7; SANISAND (85–95 %) → 1.0–1.1, never slower (T6).

### #1096, #1098 and #1095

- `mpm.f.admissible` (#1096) is 0.0–0.8 % of a step and `mpm.f.dtcrit` (#1098) 0.2–5.7 %
  (1 thread, every row). Neither dominates; no action under this WP.
- #1095: `test_cut_retry_across_connectivity_changes_3d` is green on `3d2cb6330`; no crash at 8
  threads in the 240 mode-F and 63 mode-A runs, so no minidump to post.

### Mode A (reported, no target)

Threads 1 and 8 only. S(8) 1.05–1.18 at 10⁴ (3D `bspline2` DP excepted, below) and 1.57 at 3D `gimp` elastic 10⁵ (PARDISO and
WP-192's element loops; not this WP's threading, and not bit-identical across thread counts,
which T8 does not ask of mode A). **DruckerPrager mode A fails** in 2D (10⁴, both bases) and 3D
`gimp` 10⁵: Newton does not converge after 4 step cuts near load factor 0.6. Build `c45af1407`
(before this WP) fails at the same step, t = 0.597656, so it is the block model's limit, not
this WP's. Kept as data (§4).

The one noisy mode-A row (3D `bspline2` DP 10⁴, 8 threads, 41 %) stays at 43 % with 5 repeats:
it is bimodal (234–240 against 333–338 p-steps/s), so it is reported within noise.

### Noisy rows

2D SANISAND at 10⁴ stays at 29–45 % spread even with 5 repeats, and its 10⁵ 1-thread row was
68 %. 2D elastic 10⁴ at 8 threads is 23–30 %. All are reported "within noise" where §4 says so;
none is a target row except through T6, which passes on every median.

### External reference: Kratos

`Ladruno_scripts/mpm_bench/external/kratos_column.py`, Kratos MPMApplication 10.4.4 (BSD, pip
wheels, isolated venv), the same quarter column, cells, ppc, E, ν, ρ, g and Δt; 5 warm-up + 20
timed steps, 3 repeats, nmb (`results/kratos_10.4.4_nmb/kratos_results.jsonl`). **Not matched:**
Kratos has linear hex shape functions (8 nodes per particle, against the wider `gimp` support; README item 1), a
NeoHookean law, a fixed floor instead of a Coulomb wall, and a Python-driven step loop (see
`external/README.md`). A reference point, not a like-for-like ratio.

| N (Kratos / fork) | Kratos T, t = 1 / 4 / 8 [M p-steps/s] | Kratos S(8) | fork 3D `gimp` elastic, t = 1 / 8 | fork / Kratos at t = 1 / 8 |
|---|---|---|---|---|
| 11,664 / 10,368 | 0.101 / 0.270 / 0.369 (≤ 4 %) | 3.65 | 1.40 / 6.85 | 14 / 19 |
| 93,312 / 98,496 | 0.097 / 0.299 / 0.385 (≤ 7 %) | 3.95 | 1.44 / 7.18 | 15 / 19 |

### Regression check

`Ladruno_scripts/mpm_bench/regress.py` against `baseline_nmb.json` (recorded on `3d2cb6330`,
quiet nmb); the fork checklist entry tells later MPM work packages to run it.
