# M0b measurements (WP-232g)

The harness behind ADR 232 §9 phase 2: every number of the pre-registered table in
`.claude/skills/ladruno-mpm/references/wp-briefs.md` ("232g"), reproducibly, from one build.
The definitions there bind this code; this file says how each one is measured and why each
model and setting was chosen. The worker never applies the table: the orchestrator does.

```bat
python3.12 Ladruno_scripts\mpm_m0b\run_all.py --pyd dist\bin\opensees.pyd            :: full
python3.12 Ladruno_scripts\mpm_m0b\run_all.py --smoke                                 :: N = 10^3
python3.12 Ladruno_scripts\mpm_m0b\m1_rebuild.py --only N=100000,solver=Pardiso       :: a subset
```

`--pyd` selects the binary (default: this worktree's `dist\bin\opensees.pyd`); the child refuses
to run if `import opensees` resolves anywhere else. `--only key=value,...` filters the
configurations; `--repeats`, `--warmup`, `--timed` override the pinned counts; `--out` the
directory. Results: `results/<script>_<stamp>/<script>_results.json` (the aggregated file; the
per-run JSON, the profiler HDF5 and each child's log sit beside it and are git-ignored).

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

Smoke (`--smoke`): N = 10³, one repeat, a covering design of 10 children — every solver × mode ×
material on `renormalize`, plus `exclude` and `ghost` on (qs, elastic, Pardiso).

Full matrix: 3 N × 2 modes × 2 materials × 2 solvers × 3 treatments = 72 configurations × 3
repeats × 8 steps. At 10⁴ a quasi-static PARDISO step takes 1–12 s on the plumbing build; 10⁵ is
an order of magnitude more and UmfPack slower still, so the full matrix is many hours: trim with
`--only` (e.g. UmfPack and the non-default treatments at 10⁵) if the orchestrator wishes.

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
  `_build_confined`): an isotropic leg, then an isochoric deviatoric ramp. The same deck runs
  ElasticIsotropic, DruckerPrager, LadrunoJ2 and LadrunoSANISAND; only the two strain amplitudes
  differ, so each material sits at the column's stress level (p ≈ 10 kPa) and yields early in the
  deviatoric leg. Per-call cost = `elem.update` wall (deep, per class tag) / (updates × 8 points)
  + `commit.elements` / (commits × 8). A first pass classifies every step elastic or plastic
  (DP, J2: the point on its yield surface; SANISAND: stage 0 / 1); a second pass on a rebuilt
  point profiles each run of same-class steps as one profiler run;
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
N = 10⁴ and 10⁵ (smoke: 10³ and about 2 × 10³). Private bytes at s0 (model), s1 (every command:
particles, material copies, analysis objects), s2 (after step 1), s3 (after step 2); three
profiler `-memory` windows (the commands, step 1, step 2) give the Matrix / Vector / ID bytes each
window still holds and the live-object census per class tag. The result is the slope of private
bytes against N after the first step (the brief's definition), with the same slope for every stage
and every counted item: that is the attribution (element and FE_Element matrices in the step-1
window's Matrix bytes; the material copies in s1 − s0; the grid share through `gridNodes`;
the solver's factor is the private bytes the profiler does not count).

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
