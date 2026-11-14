---
title: "WP-232 — the explicit MPM engine (mode F): user guide"
project: Ladruno
tags: [mpm, wp-232, guide]
---

# The explicit MPM engine (`mpm`) — user guide

Slice 1 of WP-232 ([[232_ladruno_mpm_general_driver_adr]]): a native explicit Material Point
Method engine, hosted as a `LoadPattern` (`LadrunoMPMPattern`, class tag 33027) and driven by
the `mpm` command in classic Tcl (`OpenSees.exe`) and Python (`opensees.pyd`). The numerics are
frozen by `tests/mpm_interface.md` (WP-232b); this guide says how to drive them.

## A minimal run

```tcl
model basic -ndm 2 -ndf 2                      ;# ndm 2 = plane strain (thickness 1), or 3
nDMaterial ElasticIsotropic 1 1.0e7 0.3
mpm grid 1 -origin 0 0 -h 0.1 -cells 50 10     ;# the host is load pattern 1
mpm body 1 box 0.5 0.0 1.5 0.4 -ppc 2 -mat 1 -rho 1000
mpm wall ymin 0.3                              ;# stationary Coulomb wall on the grid face
mpm gravity 5.0 -8.66
timeSeries Path 2 -time {0 0.2 1e9} -values {0 1 1}
mpm timeSeries 2                               ;# body force = factor(t) * gravity
set dtc [mpm criticalTimeStep]
mpm analyze 5000 2.0e-4
puts [mpm particleData x]
```

Python is the same with `ops.mpm('grid', 1, '-origin', 0.0, 0.0, ...)`.

## Commands

| Command | What it does |
|---|---|
| `mpm grid <tag> -origin <x0..> -h <h> -cells <n..> [-basis gimp\|linear]` | creates the host (load pattern `<tag>`); one per Domain in v1. `gimp` (uGIMP) is the production basis; `linear` is the cell-crossing contrast |
| `mpm body <id> box <lo..> <hi..> -ppc <n> -mat <tag> -rho <ρ>` | lattice seeding (interface §5): `n` particles per cell per axis at `x0 + (a+½)h/n`, `V0 = (h/n)^ndm`, `lp = h/(2n)`; also `disc cx cy R` (2D) and `ball cx cy cz R` (3D). Ids run 0, 1, 2, … in seeding order. Returns the count |
| `mpm shift <body> <dx..>` | translate a seeded body, before the first step |
| `mpm particle <id> <body> <mat> <ρ> <V0> -x .. -lp .. [-v ..] [-s ..]` | one explicit particle (`-s`: initial stress, Voigt, tensor shear) |
| `mpm velocity <body> uniform <v..>` / `sine <comp> <axis> <v0> <β> <x_ref>` | initial velocity of a body |
| `mpm stress <body> uniform <s..>` / `geostatic <y_top> <ρ> <g> <K0>` | initial stress of a body (2D geostatic: σyy = −ρg(y_top − y), σxx = K0 σyy) |
| `mpm fix <face> <comps>` | zero the listed grid-velocity components on a face (`xmin … zmax`; `comps` ⊂ `xyz`). A roller is the face's normal component |
| `mpm wall <face> <μ>` | stationary rigid wall with Coulomb friction μ (no penetration, friction once per step) |
| `mpm gravity <g..>` | body force per unit mass, signed |
| `mpm timeSeries <tsTag>` | the factor applied to gravity (default 1) |
| `mpm damping local <α>` | local non-viscous (FLAC) damping, α ∈ [0, 1) — for quasi-static runs |
| `mpm flip <α_F>` | FLIP fraction: 1 pure FLIP (default), 0 pure PIC |
| `mpm massCutoff <c_m>` / `mpm safety <s>` | node active iff `m_i > c_m · max m_p` (default 1e-10); Δt refused iff `Δt > s · Δt_c` (default 0.5) |
| `mpm criticalTimeStep` | `Δt_c = h / max_p (c_p + |v_p|)`, `c_p = sqrt(M_p/ρ_p)` from the material's INITIAL (elastic) tangent (interface A6): a softened current tangent would overestimate it |
| `mpm analyze <n> <Δt>` | n explicit steps; recorders run at each `Domain::commit`, then flush |
| `mpm reset` | particles and materials back to the state before the first ACCEPTED step; counts an epoch (`mpm info epoch`); clears an abort |
| `mpm particles -file <path>` | load a Contract-1 particle set (232a), all or nothing; the path is opened as given (apeGmsh's emitters anchor it to the deck's folder) |
| `mpm units <L> <F> <T> [-system <id>]` | the units the `.ladruno` group carries when no particle set declared them (default `unspecified`) |
| `mpm output <path> [-every N \| -dt Δt] [-fields F..] [-compress k] [-flush s]` / `mpm output -close` | the `.ladruno` particle group (Contract 2); fields `POSITION` (always) `VELOCITY STRESS STRAIN VOLUME GIMP_HALF_LENGTH DEFORMATION_GRADIENT`, default `POSITION VELOCITY STRESS VOLUME`; no `STATE_*` token is admitted in v1 |
| `mpm info numParticles\|epoch\|time\|summary` | `summary` = [mass, momentum.., max \|v\|, kinetic energy] |
| `mpm particleData <field>` | flat list in storage order: `id body mat m V0 V x v lp F R sig sigHat epsHat sig0` |
| `mpm debug shape <x..> -lp <lp..>` / `shape1d <h> <lp> <r>` / `activeNodes` | unit-level hooks (interface Q4) |

## One step

1. refuse if a particle centre is outside the real grid (interface A2; the error names the lowest such id) or `Δt > s·Δt_c`;
2. `Domain::applyLoad(t+Δt)`: the domain time advances and the host reads its series factor;
3. P2G of mass, momentum and forces over the particles' support nodes; nodes below the cutoff drop;
4. damping, grid momentum update, grid conditions (walls, then fixes);
5. MUSL: particle velocity (PIC/FLIP blend) and position from the grid — a step whose TRIAL
   positions would carry a centre out of the real grid is refused with nothing changed (interface
   A7), so an accepted step never leaves a particle outside — then a second projection,
   ΔF and the incremental hypo-kernel entry (Hughes–Winget, unrotated frame), the material
   update — with no active element (`ops_TheActiveElement = 0`);
6. commit the particle materials (a declared refusal aborts the run);
7. particles take their trial state; 8. `Domain::commit()` (recorders); 9. profiler rows.

The engine grid carries one ghost node layer beyond every face (interface A5), so a particle's
GIMP support is never cut: the partition of unity holds right up to the face, and a body slides
along a face without drag. A face's `fix` or `wall` acts on its face nodes and its ghost layer
alike; a face with no condition just has free ghost nodes. Node ids count the ghost layer
(interface §2): `id = (i+1) + (nx+3)·((j+1) + (ny+3)·(k+1))`, with the ghost layer at index −1.

## Materials

| Material | v1 |
|---|---|
| `ElasticIsotropic` | admitted; reproduces the interface's hypoelastic law, with initial stress as an offset |
| `LadrunoJ2` | admitted; `-autoRegularization` refused (it reads an element length; there is none at particles) |
| `DruckerPrager` (UW) | admitted; stage switch `materialState` through `parameter` (below) |
| anything else | refused at `mpm body` / `mpm particle`, naming the class tag |

**Initial stress** (`mpm stress`, `-s`) is accepted on `ElasticIsotropic` only. On LadrunoJ2 or
DruckerPrager it is refused: an additive offset would bypass the yield check, and v1 has no
yield-consistent seeding (interface Q2, 232c's design).

**Stages.** `updateMaterialStage` walks Domain elements only and never reaches particles
(ADR W12). Use the host as the parameter target — the spelling is `loadPattern` in both engines:

```tcl
parameter 1 loadPattern 1 materialState
updateParameter 1 0           ;# DruckerPrager elastic
```

## What `reset` does

`reset` (the Domain's `revertToStart`) reaches a load pattern only as a re-`applyLoad`, so the
particles **stay where they are**. `mpm reset` puts them back. For a clean restart use both.
`wipe` and `remove loadPattern <tag>` delete the host and its particles.

- **The snapshot** `mpm reset` restores is the state before the first ACCEPTED step: every
  `mpm analyze` retakes it until a step is accepted, so bodies added after a refused or zero-step
  analyze are part of it. Before any accepted step, `mpm reset` keeps the particles as they are
  (they ARE the initial state, bodies and initial conditions included) and only reverts the
  materials and clears an abort.
- **What it restores:** kinematics (x, v, F, R, V, strain feed, stress) and each material's
  committed state (`revertToStart`). **Not** stage flags or moduli set through
  `parameter … loadPattern` — DruckerPrager's `materialState` and its moduli survive a reset
  (`DruckerPrager.cpp`, `revertToStart` leaves them). Re-issue `updateParameter` after a reset when
  the run must start from the original stage.
- **Parameters** that reach particle materials are detached when the host is deleted (`wipe`,
  `remove loadPattern`) or when a reset drops particles: a later `updateParameter` does nothing
  instead of touching freed memory.
- **Abort (ADR W17):** a material refusal at commit aborts the run with the particles unmoved;
  every later `mpm analyze` is refused until `mpm reset`.
- **Initial conditions** (`mpm body`, `mpm particle`, `mpm particles`, `mpm stress`,
  `mpm velocity`, `mpm shift`) are refused after the first step, until `mpm reset`.

## Refused in v1

- a Domain with a free FE DOF (no FE–MPM seam yet): fix every node or remove the FE model;
- database `save` / `restore` while the host exists in mode F (mode A has a restart path);
- `mpm` in OpenSeesSP / OpenSeesMP;
- adding particles after the first step (until `mpm reset`).

## Output

Stock recorders cannot see particles. `mpm output` writes a `.ladruno` file of its own (232a
Contract 2): `FORMAT_VERSION` 1, `FORMAT_MINOR_VERSION` 2 (only a file that carries particles is
2; `recorder ladruno`'s FE-only files stay at 1). Each stage is particle-only — the recorder's empty
stage (`EMPTY_PARTITION = 1`, empty nodes, elements and result groups) plus `MODEL/PARTICLES` and
`RESULTS/ON_PARTICLES` — with `KIND = "transient"` and `MODE = "explicit"`.

- The file is created when the first stage opens: `mpm output` followed by `-close` with no step
  leaves no file.
- A stage opens when the first step of an epoch has passed its refusals, and its output 0 is the
  state before that step; so bodies and initial conditions set after a refused or zero-step analyze
  reach the file. `mpm reset` starts a new epoch, hence a new stage (`PARTICLE_EPOCH` + 1). The
  particle set is fixed once a stage is open; a mismatch is refused, never written short.
- **Flush:** by default the file is flushed after EVERY output (`-flush 0`): particle outputs are
  already throttled by `-every` / `-dt`, and an HDF5 file is crash-consistent only at a flush —
  between flushes a kill can corrupt it (unreadable headers, or deflated chunks that fail their
  filter), not merely truncate it. `-flush s` relaxes it to a wall-clock cadence, at that risk. The
  data are laid out so a flushed output is never rewritten (one chunk per output for `ID` and every
  `DATA`, uncompressed `STEP` / `TIME` / `OFFSETS`), and every output is written behind a barrier: a
  flush after its rows and before its `OFFSETS` entry, so an entry the crash rule counts always has
  its rows on disk. The file is also written in **SWMR mode** (latest HDF5 format) from the end of
  each stage's layout to close: without it, a kill inside the flush where a dataset's chunk index
  (a v1 B-tree) splits leaves a parent node pointing past the stored end of the file, because
  HDF5 writes the parent before the superblock, and the committed outputs the split moved become
  unreadable ("addr overflow"; 1 kill in 137 behind the barrier alone, the same with
  metadata-cache eviction disabled). A new stage closes and reopens the file to lay itself out.
  A killed run leaves the SWMR flags set: read it by the contract's reader rule (a plain
  read-only open, else `h5py.File(path, "r", swmr=True)`; never `h5clear`), then apply the crash
  rule (232a §2.5). Readers need HDF5 ≥ 1.10. The stage layout itself (at the start, and at each
  `mpm reset`) is written without SWMR: a kill inside it leaves a file that neither open accepts,
  and the reader rule's step 3 recovers a copy (`recover_particle_file`); every earlier stage
  survives. The `-dt` cadence re-arms when the Domain time goes back
  (`reset`, `setTime`).
- Outputs are written after `Domain::commit` (positions and stresses of the same step), every `N`
  steps or every `Δt` of time, and the last step of every `mpm analyze` is always written.
- Ids are int64 (≤ 2⁵³ − 1), sorted within an output, stable within an epoch. A second particle
  set whose ids meet the epoch's is refused.
- When a particle set was loaded, its `SET_ID`, units, thickness and body names are copied; lattice
  bodies are named `body<id>`.
- v1 writes FE output and particles to separate files: a mixed FE + particle stage is later work.
- Moving 232a's Contract 2 (§2) into `ladruno_schema_v1.md` (the contract's §2.9) is deferred to the
  232a follow-up.
- Classic Tcl returns `mpm criticalTimeStep`, `mpm particleData` and `mpm info time|summary` as
  full-precision doubles (the elementAPI's `%35.20f` does not round-trip).

## Mode A — implicit dynamic and quasi-static (WP-232f)

Phase 1 of the ADR (§5.3; frozen interface `tests/mpm_implicit_interface.md`). `mpm mode
implicit` or `mpm mode quasiStatic` makes every particle a `LadrunoMPMPoint` element (class tag
33026) and its active grid support Domain nodes, and **your own analysis** solves each step.
Explicit stays the default (`mpm mode explicit`); the mode is fixed at the first step of an epoch
(`mpm reset` to change it).

```tcl
model basic -ndm 2 -ndf 2
nDMaterial ElasticIsotropic 1 6.0e7 0.0
mpm grid 232 -origin 0.0 -4.0 -h 0.25 -cells 44 22
mpm mode quasiStatic
mpm body 1 box 0.0 0.0 10.0 1.0 -ppc 2 -mat 1 -rho 1000.0
mpm fix xmin xy
mpm gravity 0.0 -10.0
timeSeries Path 9 -time {0 1 1e9} -values {0 1 1}
mpm timeSeries 9
constraints Plain; numberer RCM; system UmfPack
test NormDispIncr 1e-10 30; algorithm Newton
integrator LoadControl 0.05; analysis Static
mpm analyze 20 0.05
```

For `implicit`, build `integrator Newmark 0.5 0.25` (or `HHT`, `GeneralizedAlpha`) and
`analysis Transient` instead. Python is the same through `ops.mpm(...)`.

**One step.** Supports at the committed positions; the active nodes (`m_a > c_m max m_p`) become
Domain nodes — a node that stays active keeps its tag, one that leaves is removed, one that enters
gets a NEW tag (a retired tag is never reused); P2G writes displacement 0 and the projected
velocity and acceleration and commits them; every particle element is re-pointed in place; the
grid conditions become this pattern's `SP_Constraint`s; the C++ `Domain::domainChange()` is raised
every step (the integrator re-reads the nodes only then); then `analyzeStep(dt)` (transient) or
`StaticAnalysis::analyze(1, false)` (quasi-static). G2P runs once, in the last particle element's
`commitState()`: positions move by Σ S u, F / stress take the converged trial state, velocity and
acceleration take FLIP increments (`mpm flip`, default 1); quasi-static leaves v and a alone.

- **Quasi-static dt is a pseudo time.** `mpm analyze <n> <dt>` steps your `LoadControl` with
  Δλ = dt (and its clamp to [dt, dt]) and puts your Δλ, min and max back afterwards; it says so once
  per call when that changes your value. The host's time series is evaluated at that pseudo time.
- **Retry.** A step that fails before any commit (Newton divergence, an inverted particle, a trial
  centre leaving the grid) is cut by `mpm cut <factor> <maxCuts>` (default 2, 4) and retried — no
  regrowth: the cut dt holds to the end of that `mpm analyze`, which still ends at t₀ + n·dt.
  `mpm info dtHistory` lists the accepted steps. A failure after a commit began aborts the run
  (`mpm reset`); a material refusal at commit leaves the particles unmoved (ADR W17).
- **Empty support.** `mpm emptyNodes renormalize` (the default, ADR §5.3 amendment 232f): a
  support node with m_a <= c_m max m_p is dropped and the kept shape functions are renormalised
  to a partition of unity, so the dropped node is absent rather than held. Mode A's default
  cutoff is **5e-2** (mode F's stays 1e-10); `mpm massCutoff` sets both. With a small cutoff,
  nodes held by one particle's GIMP tail make quasi-static Newton diverge. For M0b:
  `mpm emptyNodes exclude` (dropped nodes held at zero increment — it anchors free surfaces) and
  `mpm emptyNodes ghost <γ>` (kept, tied to their step-start position by a spring γ M₀ V / h²).
- **Tags.** Grid nodes and particle elements take tags above every tag in use at the first
  mode-A step (`mpm tagBase <node> <ele>` to choose). Create your FE components first.

| Command | What it does |
|---|---|
| `mpm mode explicit\|implicit\|quasiStatic` | the engine mode |
| `mpm gridvel <name> box <lo..> <hi..> v <value\|free ..> [start <t0>] [track x\|y\|z]` | prescribed grid velocity on every grid node in the closed box: Δu = v·dt per step; absent before `start`; `track` moves the box's lower bound with v. A non-zero value needs `constraints Transformation`, `Penalty` or `Lagrange` (Plain assumes homogeneous: refused) |
| `mpm condForce <name>` | the grid internal force Σ V σ·∇S over the gridvel's nodes at the last accepted step (resisting-force sign) — e.g. a footing load |
| `mpm pforce <body> box <lo..> <hi..> <f..>` | a constant force per particle (initial centre in the box), scaled by the time series |
| `mpm rayleigh <alphaM> <betaK>` | Rayleigh damping on the particle elements |
| `mpm emptyNodes`, `mpm cut`, `mpm tagBase` | above |
| `mpm info g2p\|acceptedSteps\|retries\|domainChanges\|connectivityChanges\|gridNodes\|mode` | counters (`g2p` counts G2P passes, one per accepted step) |
| `mpm info dtHistory\|gridNodeTags` | lists |
| `mpm debug failAttempt <k> [<m>]` / `mpm debug reseed default\|always\|skip` | test hooks: fail the first m attempts of step k; recreate every node each step, or skip `domainChange()` on a step with no support change (wrong on purpose: W1) |

**Refused in mode A (loudly):** any transient integrator but Newmark, HHT, GeneralizedAlpha
(`LadrunoHHT`, `LadrunoGeneralizedAlpha` included, until a test admits them);
`VariableTransient`; any static integrator but `LoadControl` (every `DisplacementControl` — never
on a grid node, and phase 1 has no FE node to drive — and every arc-length); a prescribed velocity
under the Plain handler; `betaK0` / `betaKc` Rayleigh on particles (they name the stiffness of an
earlier support); ground-motion patterns; `mpm wall` and `mpm damping local` (mode F only); free FE
DOFs (no FE–MPM seam yet); element loads on particles.

**Recorders.** No stock recorder and no `recorder ladruno` may watch grid nodes or particle
elements: their tags start above yours and churn every step (ADR W4, amendment 232f), and
`recorder ladruno` would open a model stage per step. Particle output is `mpm output` (MODE
`implicit` / `quasi-static`, KIND `transient` / `static`).

**Restart.** `database File <name>`, `save <n>` between `mpm analyze` calls (after the first
mode-A step), then `restore <n>` — in the same process after `wipe` and `model basic …`, or in a
new one — rebuild your analysis and continue. The particle state, the elements and their
materials, the grid nodes and the counters come back; re-issue `mpm output`. Mode F still refuses
a save.

**Verification (mode A).** `tests/test_ladruno_mpm_implicit_engine.py` (Zone-A, Python and classic
Tcl): the implicit self-weight column, the Newmark bar (the re-seed leg), the failed-step retry,
the no-cell-crossing re-seed (normal = forced re-creation, skipped differs), recorder synchrony,
restart both ways in both engines, every refusal, engine parity, and the data rows (cantilever QS
tip at two cutoffs, exclude against ghost, Newton iterations per step).
`tests/test_ladruno_mpm_implicit_gates.py` holds the fine-resolution gates (cantilever and footing
r2 / r3, the plate with a hole) — `LADRUNO_MPM_GATES=1`, too long for Zone-A. Both read
`tests/mpm_implicit_cases.txt` (232f-V); a `TBD` expected value is skipped loudly. The suite
joins the MPM mutation family.

## Verification

- `tests/test_mpm_engine_core_cpp.py` — the pure core (`LadrunoMPMCore.h`) and the kernel entry,
  standalone with g++, against every case of `tests/mpm_cases.txt` and
  `tests/hypo_incremental_cases.txt` (232b's independent numpy oracle).
- `tests/test_ladruno_mpm_contracts_io.py` — Contract 1 in (232a's fixture; every malformed
  variant refused with its rule code), Contract 2 out (every file passes `ladruno_format.validate`
  and `tests/mpm_contracts.py`), the round trip, output synchrony, id stability, epochs, Tcl.
- `tests/test_ladruno_mpm_engine.py` — every case through `OpenSees.exe` **and** `opensees.pyd`,
  and the two engines against each other. The MPM mutation family's suite.
- `tests/test_ladruno_mpm_lifecycle.py` — refusals, lifecycle (the reset snapshot, parameters
  after `remove loadPattern`, initial conditions after the first step, a kernel refusal changing
  nothing), the stage path, the critical Δt from the initial tangent, and a two-leg check: every
  plastic particle's stress (LadrunoJ2, DruckerPrager) against perfect J2 plasticity driven through
  that particle's own strain path. The v1 particle loop is serial, so thread-count determinism holds
  by construction; its test arrives with the threaded loop.
- The ADR W17 paths (a kernel refusal fails the step; a commit refusal aborts with particles unmoved;
  a stale refusal count does not abort) are tested on the core with a mock material
  (`mpm_engine_check.cpp`): no admitted v1 material declares a refusal.
- **Mutation gate (ADR-87 D2), family `MPM`:** `ZERO` drops the P2G internal force, `SCALE`
  multiplies the G2P position increment by 1.5. The suite is `tests/test_ladruno_mpm_engine.py`.
  Run conditions: the mutant `opensees.pyd` on `PYTHONPATH`, the default tests directory (not
  `--module-dir`, which runs pytest where the test files are not), and `LADRUNO_TCL_EXE` pointing
  at a nonexistent path in BOTH the baseline and the mutant run — a mutant build has no
  `OpenSees.exe`, so the Tcl half would otherwise run the normal binary and score as survivors.
  Measured on the Python half: ZERO 18/38 = 0.474 (floor 0.45), SCALE 22/38 = 0.579 (floor 0.55).
  The survivors are the unit hooks, the case orderings, the stress-free closed forms under ZERO,
  and the cases that never test a position increment under SCALE.
- Known: the 20%-compression self-weight column misses its frozen 0.03 (0.072, as in the oracle;
  uGIMP cell crossing) — a strict xfail per interface A4.
