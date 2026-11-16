# M0b measurements (WP-232g)

The harness behind ADR 232 §9 phase 2: every number of the pre-registered table in
`.claude/skills/ladruno-mpm/references/wp-briefs.md` ("232g"), reproducibly, from one build.
The definitions there bind this code; this file says how each one is measured and why each
model and setting was chosen. The worker never applies the table: the orchestrator does.

```bat
python3.12 Ladruno_scripts\mpm_m0b\run_all.py --smoke                                 :: N = 10^3
python3.12 Ladruno_scripts\mpm_m0b\run_all.py --plan full --dry-run                   :: list + estimate
python3.12 Ladruno_scripts\mpm_m0b\run_all.py --plan full                             :: the orchestrator's plan
python3.12 Ladruno_scripts\mpm_m0b\m1_rebuild.py --N 100000 --solvers Pardiso --treatments renormalize
```

`--pyd` selects the binary (default: this worktree's `dist\bin\opensees.pyd`); a relative path is
resolved in the parent, from the current directory, before any child starts (the children run
in their results directories). The child refuses to run if `import opensees` resolves anywhere
else. Subset flags (comma lists; a flag a script's configurations do not carry is ignored):
`--N`, `--modes` (M1 `qs,transient`; M3 `A,F`), `--mats`, `--solvers`, `--treatments`, M4's
`--cases r1,r2,r3` and `--cutoffs`; `--repeats`, `--repeats-by-N 100000=1`, `--warmup`,
`--timed`; `--timeout-s` per run; `--only key=value,...`; `--list` prints the selected runs.
**A timeout or a crash is data**: the run is kept with `ok: false`, its `reason` (`timeout` or
`error`) and its wall time, and every summary lists its `failures`.

**Affinity.** Each child sets its process affinity, before `opensees` (and MKL's threads) start,
through `SetProcessAffinityMask`: default `0xFFFF`, the eight P-cores with their hyperthreads
(logical 16–23 are the E-cores of the i7-13700KF); `--affinity none` leaves it. The requested
mask, the mask the OS reports back and the system mask are in each run's `conditions.affinity`.
The thread counts are unchanged.

**`--plan full`** is the orchestrator's run decision (2026-11-15; it moves no threshold), one
script invocation per item, each in `<out>/<item>/`:

| Item | What |
|---|---|
| `m1_pardiso_renormalize` | the decision row: N 10³, 10⁴, 10⁵ × {qs, transient} × {elastic, dp}, 3 repeats |
| `m1_umfpack_renormalize` | the same, 3 repeats at 10³ and 10⁴, 1 at 10⁵; 2700 s timeout per run |
| `m1_pardiso_ghost` | the same pairs, every N; 3 repeats at 10³ and 10⁴, 1 at 10⁵ |
| `m1_pardiso_exclude` | 10³ and 10⁴, 1 repeat (cost only) |
| `m2` | the full design, 3 repeats |
| `m3` | N 10³, 10⁴, 10⁵; modes A and F; elastic and dp; 2 repeats; slope on 10⁴ → 10⁵ |
| `m4` | r1–r3 × cutoffs {5e-2, 1e-6, 1e-10}, 1 repeat |

UmfPack × ghost / exclude is not run (it decides nothing). At N = 10⁵ M1 keeps at least 4 timed
steps. `--dry-run` prints every run and a wall-time estimate scaled from a smoke
(`--smoke-dir`, default the latest `results/run_all_smoke_*`): the stepping part of each smoke
run (its Python wall around `mpm analyze`) is scaled by (N / 10³)^p and by the ratio of step
counts, the rest of the child's wall is fixed; p = 1.5 for PARDISO and mode A (the pre-232g
plumbing build gave ×20–27 in step time for ×10.6 in N on the DP block), 2.0 for UmfPack, 1.0 for
explicit; the cantilever scales by its own particle count with p = 1.5; a run is capped at its
timeout. It does not know Newton's iteration count at 10⁵, so read it as ±2×.

Results: `results/<script>_<stamp>/<script>_results.json` (the aggregated file; the per-run JSON,
the profiler HDF5 and each child's log sit beside it and are git-ignored).

## Pinned conditions (every result JSON carries them)

- One build: each child prints `ladrunoBuild` from the binary it imported.
- Machine i7-13700KF (8 P + 8 E), 64 GB (recorded as pinned, and what the parent sees).
- Python 3.12. `OMP_NUM_THREADS` and `MKL_NUM_THREADS` set in the child's environment before
  `opensees` is imported: 8 for implicit (M1, M3 mode A, M4), 1 for explicit (M2, M3 mode F),
  plus one elastic M2 row with both unset ("default").
- **Every run is a fresh process** (`<script> --child <cfg.json>`).
- Median of 3 repeats with the spread (max − min) / median. Warm-up steps are excluded and
  stated: M1 2 warm-up + 6 timed (the smoke too: 10³ is cheap), M2 5 + 20 (smoke 2 + 5); the
  smoke runs one repeat.
- No build or other heavy job during a timed run: the orchestrator launches the full runs.

## Instrumentation (profiler scopes only)

`OPS_PROFILE_SCOPE` is a no-op when the profiler is off (`SRC/utility/profiler/ProfilerMacros.h:77`:
the name is passed to `ScopedTimerOpt` only when `theProfiler().enabled()`), so the scopes
change no result.

| Where | Scopes |
|---|---|
| `StaticAnalysis.cpp` (vanilla) | `domainChanged` around the call in `analyze()`; inside `domainChanged()`: `dc.clearAll`, `dc.handle`, `dc.numberDOF`, `dc.setSize` ⊃ `dc.s.graph`, `dc.integrator`, `dc.algorithm` |
| `DirectIntegrationAnalysis.cpp` (vanilla) | adds `dc.clearAll`, `dc.integrator`, `dc.algorithm` to its existing `dc.*`, so the rebuild splits with no remainder (Newmark re-reads every node in its `domainChanged`, W1) |
| `LadrunoMPMImplicit.cpp` | `mpm.a.step` (one mode-A step: prepare, attempts, cuts); `mpm.a.prepare` ⊃ `mpm.a.support`, `mpm.a.sp` (spec, SP removal, SP rebuild), `mpm.a.churn` (grid-node retire / create), `mpm.a.p2g`, `mpm.a.repoint` (element re-pointing and first add), `mpm.a.domainChange` (the raise); `mpm.a.g2p` (the batched G2P) |
| `LadrunoMPMPattern.cpp` | `mpm.f.material`: a **deep** scope around `eng.advance`; each particle material call (`setTrialStrain`, `commitState`) runs under an `ElemScope` that folds its wall into that scope's `elem_by_type` row for the material's class tag |

**Why mode F's material share is measured per call.** The core (`LadrunoMPMCore.h`) is
OpenSees-free by design (it is compiled standalone by `tests/test_mpm_engine_core_cpp.py`) and
calls the material inside the particle loop, interleaved with the kinematics, so no scope can
bracket "the material update" as one region. The adapter (`LadrunoMPMPattern::Mats`) times each
call instead, through the profiler's existing per-element mechanism (`ElemScope`, wall clock
only). It is inert unless the profiler runs with `-deep`, which is why M2 runs every
configuration twice: `coarse` for throughput (no per-call clock reads) and `deep` for the share.
The `mpm.f.material` node's own wall is the whole advance; the material time is the sum of its
`elem_by_type` rows.

**Threads.** The profiler keeps one tree per OS thread, each rooted at `root`, and the report
merges the trees by name: a scope opened on a worker thread would appear at the top level (not
under its parent), with its wall summed over threads. None of the measured scopes opens on a
worker: mode F's particle loop is serial in v1; `LadrunoMPMPoint` opts into none of the
`ladrunoThreadSafe*` paths, so assembly, `Domain::update` and the commit stay serial; PARDISO's
MKL threads open no scope. Every decomposed step lists any top-level scope outside
`{mpm.a.step, step}` as `foreign_top_level`, so a threaded scope would show, not hide.

## M1 — rebuild against one Newton iteration (`m1_rebuild.py`)

**Model: the laterally confined soft block** (`models.block_commands`). A 1 m cube of particles
(`ppc 2`, uGIMP) in a grid of `n × n × (n+2)` cells, `h = 1/n`, rollers on the four sides, bonded
base (`zmin xyz`), gravity 10 m/s², ρ = 1000. N = 8n³: n = 5, 11, 23 give N = 1000, 10648,
97336. Why this model:

- *it converges at every N in both modes*: confined compression has no limit load; the only
  free surface is the top, so no particle hangs on one support node;
- *it crosses cells during the timed steps*: E = 5×10⁴ Pa makes ρgB/E = 0.2 (top settlement
  about 0.074 m = 0.37 h at n = 5, 1.7 h at n = 23); particles start on a uGIMP support boundary
  (x₀ + h/4 with lp = h/4), so their support changes as soon as they move. Each step records
  whether `mpm info connectivityChanges` moved (a 0 / 1 per step);
- *one geometry, every N*: N scales by h alone, so a slope against N is a slope against resolution.

Loading: quasi-static, `LoadControl` ramp reaching 0.75 g after 8 steps (Δλ = 0.09375, whatever the
run length: the smoke takes the same increments and stops early). Transient, gravity applied
suddenly, Newmark 0.5 / 0.25, Δt = T/16 with T = 4B/c_p (the confined column's quarter-wave
period), so the 8 steps cover the first loading half-cycle.

Materials: ElasticIsotropic (E, ν = 0.3) and **DruckerPrager** (UW, K and G of the same E, ν;
σY = 400 Pa, ρ = 0.21 — the 232f footing's cone — and ρ̄ = 0, non-associated as the footing, so the
tangent is unsymmetric). Under K0 loading this cone yields in the lower part of the block (recorded as
`yield_fraction_end`: 0.4 quasi-static at 0.75 g and 0.8 transient, at 10³ on the plumbing
build); the top stays elastic, so no particle sits at the apex.

Why 0.75 g and σY = 400 Pa: on build `84ff95616` (a pre-232g binary, plumbing only), the
quasi-static DP block at σY = 200 Pa and a full-gravity ramp needed 6 → 47 Newton iterations
over the last steps and cut its last step three times (10³), and 28 iterations at 10⁴; at
0.75 g and σY = 400 Pa it takes 4 → 20 (10⁴) with no cut. The elastic block takes 4 → 7.

Analysis: `constraints Plain`, `numberer RCM`, `system Pardiso | UmfPack`, `test NormDispIncr
1e-6 50`, `algorithm Newton`, `mpm massCutoff 5e-2` (the mode-A default, set explicitly for all
treatments). The tolerance is 10⁻⁶ of the block size, 10⁻⁵ of the top settlement: the hypo kernel's
tangent is not the exact linearisation (232f), Newton converges linearly, and a tighter
tolerance only adds iterations — which lowers f. The looser choice is the one that does not
flatter mode A.

Empty-node treatments: `renormalize` (default), `exclude`, `ghost 1e-2` (the 232f study's γ), all
at cutoff 5e-2. M1-t's observable is `settlement_top`: the mean vertical displacement of the top
particle layer (initial z ≥ B − h/2), recorded after every step.

**Decomposition** (`common.m1_step`, per step, from the profiler tree; every `mpm analyze 1 dt`
is its own profiler run):

- `t_rb = mpm.a.prepare + domainChanged + Σ soe.symbolic outside domainChanged` — UmfPack's
  symbolic runs in `setSize` (inside `domainChanged`), PARDISO's in the first solve after a new
  pattern (inside `linearSolve`); each is counted once;
- `t_it = (formTangent + (linearSolve − soe.symbolic) + formUnbalance + update + convTest) / k`,
  `k` = the step's `linearSolve` calls (= Newton iterations; `testIter` is recorded beside it).
  Newton's first `formUnbalance` (before the loop) is spread over the k iterations;
- `R = t_rb / t_it`, `f = t_rb / wall(mpm.a.step)`;
- the rebuild's sub-phases: the six `mpm.a.*` children (+ `mpm.a.prepare.other`), `dc.clearAll`,
  `dc.handle`, `dc.numberDOF`, `dc.s.graph`, `dc.setSize.soe` (setSize minus graph and minus any
  symbolic), `dc.integrator`, `dc.algorithm`, `dc.other`, and `soe.symbolic`; plus
  `symbolic_share_of_t_rb` and `formTangent_share_of_t_it` (M4's item-15 rule reads it at 10⁵).

Per run, over the timed steps: t_rb and k are means, t_it = Σ iteration time / Σ k, R = t_rb /
t_it, f = Σ t_rb / Σ step wall. Across repeats: the median and the spread. `fits` gives the log-log
slope against N of R, t_rb, t_it and f per (mode, material, solver, treatment). A step that was cut
(`attempts` > 1, `substeps` > 1, `retries`) stays in the record and is flagged.

**Changed against unchanged connectivity** (for the C-cheaper-first row): every step carries
whether `mpm info connectivityChanges` moved during it; each run also reports `t_rb_changed_ns`,
`t_rb_unchanged_ns`, `symbolic_changed_ns`, `symbolic_unchanged_ns`, the counts `n_changed` /
`n_unchanged` and the rebuild sub-phases over the unchanged steps alone. On a step whose
connectivity did not change, every nanosecond of t_rb is avoidable without mode C.

Smoke (`--smoke`): N = 10³, one repeat, a covering design of 10 children — every solver × mode ×
material on `renormalize`, plus `exclude` and `ghost` on (qs, elastic, Pardiso).

The script's own default is the whole cross product (3 N × 2 modes × 2 materials × 2 solvers × 3
treatments); `run_all --plan full` runs the orchestrator's subset of it.

## M2 — explicit throughput (`m2_explicit.py`)

**Model: a 3D granular column collapse**, a quarter column a × a × 2a (a = 1 m) on symmetry
rollers (`xmin x`, `ymin y`), a Coulomb floor (`mpm wall zmin 0.4`), in a 4a × 4a × (2a + 2h)
grid, ρ = 2000, g = 9.81, ppc 2, uGIMP. N = 16 m³ with m cells across a: m = 4, 9, 18 give
1024, 11664, 93312. Δt = 0.25 Δt_c (`mpm criticalTimeStep`; the host refuses above 0.5). Three
twins: ElasticIsotropic (E = 1 MPa, ν = 0.3), **DruckerPrager** (same K, G; φ = 30° triaxial match
ρ = 2√2 sin φ / (√3 (3 − sin φ)), ρ̄ = 0, σY = 10 Pa — practically cohesionless), LadrunoJ2
(same K, G, σy = 2 kPa, perfectly plastic). The fraction of particles on their yield surface is
recorded at the start and end of the timed window.

Each configuration runs as `coarse` and `deep` (above). The warm-up steps are one `mpm analyze`
call; the timed steps are ONE `mpm analyze` call under one profiler run with `-perStep` (the
per-step wall series), so no report interrupts the stepping. Throughput = N × steps / Σ `step`
wall; also from the per-step median (robust to one descheduled step) and from the Python wall.
Material share (deep) = Σ `mpm.f.material` rows / Σ `step` wall; non-material share = 1 − it.

**The SANISAND projection.** SANISAND is not admitted at particles in v1, so its row is composed:

- *the single-point driver* (`models.point_build`): one `stdBrick` with every DOF prescribed
  (zero free equations), on the fork's confine-first SANISAND deck (`tests/test_ladruno_sanisand.py`,
  `_build_confined`), in three legs: `confine` (isotropic, to the column's p ≈ 10 kPa, 50 steps),
  `pre` (isochoric, past first yield, 50 steps), `measure` (200 isochoric steps of a fixed strain
  increment). The same deck runs ElasticIsotropic, DruckerPrager, LadrunoJ2 and LadrunoSANISAND;
  only the confinement and pre-yield amplitudes differ per material. Per-call cost =
  `elem.update` wall (deep, per class tag) / (updates × 8 points) + `commit.elements` /
  (commits × 8). A first pass classifies every step elastic or plastic (DP, J2: the point on its
  yield surface; SANISAND: stage 0 / 1); a second pass on a rebuilt point profiles each run of
  same-class steps of one leg as one profiler run. The **plastic** per-call cost is taken from the
  `measure` leg only; the elastic one from every elastic step;
- *the strain increment per call*: the `measure` leg's increment norm is matched to the DP
  column's **median per-step particle strain-increment norm** (`strain_increment_norm.median`:
  per particle, the norm of epsHat(n+1) − epsHat(n) with engineering shear,
  √(Σ ε_ii² + ½ Σ γ²), over one unprofiled step right after the timed window; the median over
  the yielded particles is recorded beside it), from the largest N whose DP coarse column ran.
  An isochoric step of axial strain a has norm a √1.5. Every driver material also runs at ×0.5
  and ×2 of that increment (`driver_plastic_per_call_band_ns`, `sanisand_over_dp_ratio_band`):
  a per-call cost that moves with the increment is an increment artefact, one that does not is
  the material's. Each leg's increment norm per step is in the result (`path.increment_norm_per_step`);
- *the composition*: c(X, φ) = (1 − φ) c_el(X) + φ c_pl(X), with φ(DP) the DP column's measured
  yield fraction and φ = 1 for SANISAND (a bounding-surface model loads plastically; φ(DP) is
  given beside it). **literal**: t(X) = non-material(elastic column, deep) + c(X, φ).
  **difference**: t(X) = t(elastic column, coarse) + c(X, φ) − c_el(elastic) — the stdBrick's
  own per-point strain work, inside every driver figure, cancels. Both forms are validated on DP
  (`dp_miss`, the rule's 25 %); both are reported, the orchestrator picks.

Why the two-pass, block-profiled driver: on the plumbing build, profiling every driver step as its
own run (a report and a reset between steps) doubled the elastic per-call cost (85 → 162 ns): the
report evicts what the next call needs. Blocks of steps measure the call, not the cold start.

## M3 — memory per particle-element (`m3_memory.py`)

Mode A: the M1 block (quasi-static, PARDISO, renormalize). Mode F: the M2 column. Elastic and DP,
N = 10³, 10⁴ and 10⁵, slope on the two largest (smoke: 10³ and about 2 × 10³). Private bytes at s0 (model), s1 (every command:
particles, material copies, analysis objects), s2 (after step 1), s3 (after step 2); three
profiler `-memory` windows (the commands, step 1, step 2) give the Matrix / Vector / ID bytes each
window still holds and the live-object census per class tag. The result is the slope of private
bytes against N after the first step (the brief's definition), with the same slope for every stage
and every counted item: that is the attribution (element and FE_Element matrices in the step-1
window's Matrix bytes; the material copies in s1 − s0; the grid share through `gridNodes`;
the solver's factor is the private bytes the profiler does not count).

**The analytical count** (`analytic`, per run; `analytic_per_particle` and
`measured_minus_analytic` beside each slope) adds, per particle: the material copy
(ElasticIsotropicThreeDimensional: 2 Vector(6); DruckerPrager3D: 8 Vector(6), 4 Matrix(6,6), 26
scalars), the particle state held twice (state and `mpm reset` snapshot, 57 doubles + ids), and
- mode A: the `LadrunoMPMPoint` (its ID, node pointers, S, ∇S, v̂, â, ghost springs, g, f_int,
  u_T, two Vectors of 3·nen; the tangent is a class-wide pool), the `FE_Element` (two IDs, and
  its OWN numDOF² tangent plus residual only when numDOF = 3·nen > 64 — `FE_Element.cpp:188-219`,
  `MAX_NUM_DOF` 64; else the class pool), at the run's measured support-size histogram
  (`support_histogram`, from `eleNodes`), and the grid share (`gridNodes` / N × Node + DOF_Group +
  host map entries);
- mode F: the support arrays (node id, S, ∇S per support node) at the counted support sizes
  (from each particle's position and l_p), the trial arrays (43 doubles), and the dense grid's
  share (19 doubles + 2 flags per grid node, ghost layer included).
Object and heap-block sizes are estimates for MSVC x64 (`ANALYTIC_NOTE` in the result). What the
count leaves out is the disagreement: in mode A the SOE matrix and PARDISO's factor (MKL-internal),
in both modes allocator slack. On the smoke (10³ → 1.7 × 10³) mode F agrees within 3–6 %; mode A's
measured slope is about twice the count, the rest being the solver.

## M4 — the cantilever data (`m4_cantilever.py`)

`cantilever_qs_r1..r3` exactly as the 232f harness builds them (`tests/ladruno_mpm_implicit_harness.py`
`model_commands` / `analysis_commands`: Plain, RCM, UmfPack, NormDispIncr 1e-10 30, Newton,
LoadControl), under `mpm emptyNodes renormalize` with cutoffs 5e-2, 1e-6 and 1e-10. Per step:
Newton iterations, retries (a divergence and its step are data), the tip (`tip_disp_y`: mass-weighted
mean of y − Y over X ≥ x_from), formTangent's share of t_it. `tip_comparison` gives each
cutoff's tip against 5e-2's relative to the 232f gate tolerance (5 %). On the plumbing build r1
reproduces the 232f tip at 5e-2 (−2.37476) and stops at step 7 at both 1e-6 and 1e-10.

## Files

`common.py` (binary, children, memory, profiler HDF5, the M1 decomposition), `models.py` (block,
column, single-point driver), `m1_rebuild.py`, `m2_explicit.py`, `m3_memory.py`,
`m4_cantilever.py`, `run_all.py`. The parent reads the profiler HDF5 with h5py; the child never
imports h5py (opensees.pyd and h5py each bring an HDF5 library).
