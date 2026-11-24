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
| `mpm grid <tag> -origin <x0..> -h <h> -cells <n..> [-basis gimp\|linear\|bspline2\|bspline3]` | creates the host (load pattern `<tag>`); one per Domain in v1. `gimp` (uGIMP) is the production basis; `linear` is the cell-crossing contrast; `bspline2` / `bspline3` are the quadratic / cubic B-splines (WP-232i, below; ≥ 2 cells per axis) |
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
| `LadrunoSANISAND`, `LadrunoHS` | admitted (WP-232k): see "Soil particles" below |
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
| `mpm gridvel <name> box <lo..> <hi..> v <value\|free ..> [start <t0>] [track x\|y\|z]` | prescribed grid velocity on every grid node in the closed box: Δu = v·dt per step; absent before `start`; `track` moves the box's lower bound with v. For a rigid footing pushed into soil, start the lower bound half a cell below the footing base (`lo = −h/2`): a bound at the base loses the soil before each node row joins (WP-244, #1112). A non-zero value needs `constraints Transformation`, `Penalty` or `Lagrange` (Plain assumes homogeneous: refused) |
| `mpm condForce <name>` | the grid internal force Σ V σ·∇S over the gridvel's nodes at the last accepted step (resisting-force sign) — e.g. a footing load |
| `mpm pforce <body> box <lo..> <hi..> <f..>` | a constant force per particle (initial centre in the box), scaled by the time series |
| `mpm rayleigh <alphaM> <betaK>` | Rayleigh damping on the particle elements (applied at once to live ones; `0 0` clears it) |
| `mpm emptyNodes`, `mpm cut`, `mpm tagBase` | above |
| `mpm info g2p\|acceptedSteps\|retries\|domainChanges\|connectivityChanges\|gridNodes\|mode` | counters (`g2p` counts G2P passes, one per accepted step) |
| `mpm info dtHistory\|gridNodeTags` | lists |
| `mpm debug failAttempt <k> [<m>] [afterTrial]` / `mpm debug reseed default\|always\|skip` / `mpm debug tangentCheck <scale> <eps>` | test hooks: fail the first m attempts of step k (`afterTrial`: after the material update, so the retry starts from a reverted trial state); the element tangent against central differences; recreate every node each step, or skip `domainChange()` on a step with no support change (wrong on purpose: W1) |

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

**`reset` in mode A.** A plain `reset` does NOT reach the particles (it warns once): their state is
the host's, and reverting only their materials would drop plastic history under deformed
positions. Use `mpm reset` (particles, materials and the grid together). A plain `analyze` cannot
run the particle elements either: they move only inside `mpm analyze`.

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

## B-spline bases, images, the particle body force (WP-232i)

Frozen interface: `tests/mpm_bspline_interface.md` (B1–B12); ADR 232 §8 P3.

**The bases.** `-basis bspline2` (quadratic, support `|r| < 3h/2`, 3 nodes per axis) and
`-basis bspline3` (cubic, `|r| < 2h`, 4 nodes per axis): node-centred uniform B-splines on the same
vertex lattice as GIMP, with the same single ghost layer and node ids. `S(0) = 3/4`,
`S(±h/2) = 1/2`, `S(±h) = 1/8` (quadratic); `S(0) = 2/3`, `S(±h) = 1/6` (cubic); C¹ (cubic C²).
Admissibility is the centre rule (`lp` is ignored, as for `linear`). A B-spline grid needs at least
two cells on every axis (`mpm grid` refuses fewer). In mode A every particle element has the full
3ᵈ (quadratic) or 4ᵈ (cubic) support, so a 3D quadratic run sits at the top of mode A's memory per
particle-element (27 nodes).

**The boundary treatment (K4): images on `fix` faces.** A B-spline's support reaches the ghost
layer from inside the grid, so zeroing the ghosts (what `fix` does on GIMP and linear grids) would
leave a first-order boundary layer. On a B-spline grid every ghost behind a face with `fix` is the
**image** of its mirror node: for a component the face fixes, its velocity is the mirror's negated
(odd); for a free component, the mirror's (even). This is the method of images — the body plus its
reflection across a symmetry plane — so a roller face (fix the normal component) is reproduced
exactly, and a clamped face (fix all) to second order. A face that fixes only tangential
components is an **antisymmetry** plane: the image is odd in the tangential components and even in
the normal one. The face-layer node lies on that plane: its fixed components carry no momentum or
force, and the active-node cutoff tests its mass times 2^k (k = the number of `fix` faces through
it), the mass the body plus its reflections would give. With both, a half model with a roller on
its symmetry plane equals the full symmetric model to round-off (G8 in mode F, G8-A in mode A at
the default cutoff, G8-F-cut at `massCutoff 0.05`). In mode A the ghosts behind a `fix` face are never
Domain nodes: their shape functions are folded into their mirrors per component (no MP constraint,
so any constraint handler works). GIMP and linear grids are unchanged. **`wall`** faces keep the
slice-1 treatment (face and ghost layers alike) and are first order on a B-spline grid; a face with
both `fix` and `wall` is refused there, and so is a `gridvel` box reaching the ghost layer behind a
`fix` face. A prescribed B-spline node velocity (`gridvel`) is a control-point value, not the
material velocity at that node.

| Command | What it does |
|---|---|
| `mpm particleBodyForce <ndm·N values> \| -clear` | a body force per unit **mass** for every particle, one flat list in `mpm particleData` order, **not** scaled by the load factor; set it between any two `mpm analyze` calls (each step, for a manufactured solution); held until replaced; cleared by `mpm reset`; carried by a mode-A database restart; a stale length (particles added since) is refused at `mpm analyze`. Mode F adds it to gravity in P2G; mode A to the particle element's external force. A step that mode A cuts and retries keeps it constant over its sub-steps |
| `mpm particleVelocity <ndm·N values>` | the initial velocity of every particle, before the first step (kept by the `mpm reset` snapshot) |
| `mpm info basis` | 0 gimp, 1 linear, 2 bspline2, 3 bspline3 |
| `mpm debug shape <x..> [-lp ..] -folded` | per own node: id, then per component `S′`, `∇S′` (the images' fold) |
| `mpm debug shape1d <h> <lp> <r> -basis bspline2\|bspline3` | the 1D factor (default `gimp`) |

**Restart format.** Since WP-232i the mode-A host's restart pack is version 2 (it carries the
particle body force). A database saved by an earlier build (version 1) is refused loudly at
`restore`; re-run the model from its deck instead.

**Tcl.** Both list commands read through the shared argument cursor, so classic Tcl accepts them
(`mpm particleBodyForce {*}$b`), but computing a body force per particle per step is a Python idiom;
the manufactured-solution harness is Python only.

**Contracts.** The `.ladruno` particle group writes `GRID/BASIS = "bspline2" | "bspline3"` (232a §2
as amended; no minor bump — a reader older than the amendment refuses the value by name). Contract 1
is unchanged: `mpm particles -file` into a B-spline grid is refused (the set's basis differs from
the deck's).

**Verification (WP-232i).** `tests/mpm_bspline_cases.txt` holds every row with its pre-registered
gate. `tests/test_mpm_engine_core_cpp.py::test_cpp_mpm_bspline_core` (g++, `tests/mpm_bspline_core_check.cpp`):
G1, G1-K4, G8 and the body force on the core. `tests/test_mpm_mms_reference.py`: the manufactured
solution's derivation (`tests/mpm_mms.py`) by centred differences, never against MPM.
`tests/test_ladruno_mpm_bspline_basis.py`: G1, G1-K4, G7, G8 through the binaries, classic Tcl too.
`tests/test_ladruno_mpm_bspline_mms.py`: the MMS rows on 1/h = 8, 16, 32 (Zone-A; in the MPM
mutation suite). The full ladders are the warrant run, `tests/mpm_bspline_warrant.py`, results in
`tests/data/mpm/bspline/` (checked by `tests/test_mpm_bspline_warrant_results.py`); method and
tables in [[232i_mpm_bspline_mms]].

## Anti-locking: the cell-centre F-bar (WP-232j)

Frozen interface: `tests/mpm_antilocking_interface.md` (J1–J13); ADR 232 §9 phase 4.

**What it does.** With `mpm antiLocking fbar`, every particle's deformation-gradient increment
of the step is rescaled so its volume change is the one of its background cell:
`ΔF̄ = (det ΔF_c / det ΔF)^{1/d} ΔF`, where `ΔF_c` is the increment of the bilinear (2D) /
trilinear (3D) field of the cell's corner nodes at the cell centre and `d = ndm` (plane strain:
the in-plane 2×2 block, `F_zz = 1`). That is de Souza Neto et al. (1996)'s F-bar, on the step
increment as in implicit MPM (Coombs et al. 2018), with the background cell as the averaging
domain: every particle of a cell gets the same volumetric increment, one constraint per cell, and
nearly incompressible flow (J2 plasticity, non-dilatant Drucker–Prager, ν → 0.5) stops locking at
2×2 particles per cell. Only the increment is modified; the carried `F`, its volume `V = det F V₀`
and the stress follow from it. The internal force is the one work-conjugate to that kinematics,
`f_a = V [σ g_a + p̄ (ĝ_a − g_a)]` with `p̄` the in-plane (2D) or full (3D) mean stress and `ĝ`
the cell centre's corner gradients — a B-bar force (Hughes 1980) beside the F-bar kinematics — so
the internal power is `V σ : ℓ̄` and the scheme keeps the energy balance (gate row L9). Without it
(the 232j round-2 build) a soft free body under FLIP gained energy without bound. One rule for `gimp`, `linear`, `bspline2` and `bspline3`, in mode F
(after the grid update, before the material call) and in mode A (inside the particle element).

**Where it does nothing.** A homogeneous (affine) field is unchanged (to 1e-12). A particle keeps its
own increment (the fallback) only when a corner of its cell is not an active node: a corner below
the mass cutoff (mode F's 1e-10, mode A's 5e-2 of the largest particle mass) or an inactive `ghost`
node. The corners of a particle's cell always carry that particle's own mass, so this almost never
happens: it did not happen once in the 232j warrant records (the free surface's cells keep the
averaging). `mpm particleData antiLocking` (0 off / 1 averaged / 2 fallback, last accepted step) and
`mpm info antiLockingCounts` (averaged, fallback particle-steps this epoch) report it.

**Mode A.** The corners of a particle's cell always lie in its own support, so the particle
elements keep their node sets, sizes and memory (about 0.4 kB more per element, on or off). The
tangent gains the consistent anti-locking terms; it is symmetric for small-strain elasticity but
**unsymmetric in general** (finite strain, plasticity): use an unsymmetric system (UmfPack,
SparseGEN, FullGeneral, `system Pardiso` with matType 0, `system Mumps -matrixType 0`).
`mpm analyze` refuses every symmetric-storage system with `fbar`: `BandSPD`, `ProfileSPD`,
`SProfileSPD`, `SparseSPD`, the distributed SPD systems, and PARDISO or MUMPS with a symmetric
matrix type. The restart pack is version 3 (it carries the setting); a version-2
database (232i) restores with anti-locking off.

| Command | What it does |
|---|---|
| `mpm antiLocking off \| fbar` | default `off` (every answer bit for bit as before); fixed for an epoch like `mpm mode` (change it after `mpm reset`) |
| `mpm info antiLocking` | 0 off, 1 fbar |
| `mpm info antiLockingCounts` | [averaged, fallback] particle-steps over the accepted steps of this epoch |
| `mpm particleData antiLocking` | per particle, the last accepted step's path: 0, 1 averaged, 2 fallback |
| `mpm debug affineCheck <A: ndm²> <b: ndm>` | L1 hook: the affine field A x + b on every support node; returns [max rel change, max \|Δε̂\|, averaged, fallback, refused] |
| `mpm debug tangentCheck <scale> <eps> -antiLocking off` | the 232f tangent check with the F-bar suppressed on the same elements (L2) |
| `mpm debug packVersion 2\|3` | the restart layout of the next database save only (test hook) |
| `mpm debug failDump <path> [<k>] \| off` | on every mode-A attempt that fails before commit, append the k particles with the largest residual at their support nodes: time, attempt, id, x, residual, trial stress, tangent norm / initial, apex flag (< 1e-2), path, det ΔF_c / det ΔF |

**When to use it.** Plasticity with isochoric or nearly isochoric flow (J2, Tresca, Drucker–Prager
with small dilatancy), and elasticity near ν = 0.5, at 2×2 (3×3×3) particles per cell. It is off by
default; flipping the default is the owner's call. It changes nothing for a homogeneous
deformation, so a compressible elastic body converges at the same order (gate row L5).

**Footing runs (K6).** The footing cases run the same solution ladder as the FE references, per
accepted step: KrylovNewton (NormUnbalance 1e-6, 40), NewtonLineSearch (1e-6, 60), KrylovNewton
(1e-5, 80, "relaxed"), then halve the step (floor 5e-4 of pseudo-time). The harness drives it
(`tests/ladruno_mpm_antilocking_harness.py::run_ladder`): `mpm cut 2 0`, so a failed try reverts to
the last accepted step and the next try (`algorithm`, `test`, `mpm analyze 1 h`) starts clean.

**Verification (WP-232j).** `tests/mpm_antilocking_cases.txt` holds every row with its
pre-registered gate. `tests/test_ladruno_mpm_antilocking_kinematics.py`: L1 on the header (g++,
`tests/mpm_antilocking_core_check.cpp`) and through the binaries, both modes, every basis, 2D and 3D;
L8 (the setting's round trip, classic Tcl = Python, refusals, the restart pack).
`tests/test_ladruno_mpm_antilocking_tangent.py`: L2. `tests/test_ladruno_mpm_antilocking_switch.py`:
L7's switch oracle (in the MPM mutation suite). The long rows (L3, L4, L5) are the warrant run,
`tests/mpm_antilocking_warrant.py`, results in `tests/data/mpm/antilocking/` (checked by
`tests/test_ladruno_mpm_antilocking_results.py`); L4 also replaces 232f's footing xfail in
`tests/test_ladruno_mpm_implicit_gates.py`. Method and tables: [[232j_mpm_antilocking]].

## Soil particles — LadrunoSANISAND and LadrunoHS (WP-232k)

Interface `tests/mpm_soil_interface.md` (D1–D23); results [[232k_mpm_soil_models]].

| Material | Admitted | Options |
|---|---|---|
| `LadrunoSANISAND` | 3D and PlaneStrain. Mode A on every basis (quasi-static runs only measured); **mode F on `bspline2` / `bspline3` only** (K0, below) | IntScheme 1, 2, 129 (others refused); `-implex`; floors `-Pmin`, `-Presidual`, `-pRe` at their defaults unless you choose (ADR-93: the owner's) |
| `LadrunoHS` / `LadrunoHSSmall` | 3D and PlaneStrain, both modes, every basis | `-smallStrain`, `-sigmat`, `-implex`; **`-sig0` / `-sig03` refused** |

- **Serial.** Neither model is on the thread allowlist; with `ladrunoThreads N` mode A prints the
  usual "not on the thread-safe allowlist … SERIAL" warning and the answer is unchanged.
- **No initial stress.** Soil particles take no `mpm stress` offset (interface Q2). The geostatic
  state comes from staging: stage 0 (elastic) under the gravity ramp, then the flip.
- **The stage switch is `mpm stage <matTag> <s>`**, both modes. It walks every particle material
  of the tag through the material's own `materialState`, so each particle gets the full switch
  (SANISAND's `Elastic2Plastic` reset included). `s` passes through unchanged: SANISAND 0/1, HS
  0/1, DruckerPrager 0/1/2 (UW DruckerPrager's 1 is still elastic; 2 is elastoplastic). A tag with
  no stage is refused. `parameter <p> loadPattern <mpm> materialState|updateMaterialStage <tag>` +
  `updateParameter <p> <s>` and, in mode A, `updateMaterialStage -material <tag> -stage <s>` route
  into the same walk.
- **SANISAND's stage is class-wide** (one static for every SANISAND tag in the process, FE and MPM).
  If it moves by a route that missed the particles (an FE element answering `updateMaterialStage`
  first, another tag's flip), `mpm analyze` is refused until `mpm stage` walks them. HS keeps its
  stage per tag and flips every instance lazily with the full transition.
- **Mode F needs a B-spline basis for SANISAND.** On `gimp` / `linear` the explicit cell-crossing
  noise excites the two-rows-per-cell stress pattern of ppc 2, and SANISAND locks it in (alpha
  moved in 70 of 160 particles; plasticity is not separated from a p-dependent elastic instability):
  the K0 meter gave A = 7.95 against DruckerPrager (bspline2: 0.66). `mpm analyze` refuses it,
  naming K0. Mode A is admitted on every basis (A 0.99 gimp, 1.06 bspline2), on quasi-static runs
  only. The admitted K0 passes engaged no SANISAND plasticity (alpha never moved); the plastic
  meter is K0p.
- **Add all bodies before `mpm stage`.** A SANISAND or HS particle added after its tag was staged to
  a non-zero stage would miss the walk, and is refused. `mpm stage` takes 0 or 1 for SANISAND and HS
  (DruckerPrager 0/1/2). Restart packs older than 4 are refused with soil particles.
- **Symmetric systems** (ProfileSPD, BandSPD, SparseSYM, PARDISO/MUMPS symmetric) are refused in
  mode A with a SANISAND or HS particle unless `-implex` (its tangent is the frozen isotropic
  elastic one).
- **Refusals.** A trial refusal fails the step (mode F) or cuts it (mode A), naming particle,
  material and tag. A commit refusal (SANISAND `-implex` companion at its `-maxSubsteps` cap)
  aborts the run with the particles unmoved (W17) until `mpm reset`. HS has no commit refusal at
  particles.
- **Queries:** `mpm particleData sigZZ`, `mpm particleResponse [-size] <args>` (e.g. `state`:
  SANISAND alpha at 6–11, e at 24), `mpm info materials`.
- **Particle stress under `-implex` is the trial stress.** `particleData sig` / `sigHat` keep the
  extrapolated stress the step was equilibrated with, not the material's committed (corrected)
  stress. In the S1 one-point legs the gap is 0.035–0.233 relative from the first plastic step.
  Read the committed stress with `mpm particleResponse stress` (quirk
  `WP-232k-implex-particle-stress-is-trial`).
- **Restart** (mode A): the restored host re-asserts SANISAND's stage flag (pack v4).
- **Cost** (M0b upper bounds): a SANISAND plastic call 0.6–1.2 ms; serial.
- **Grid motion for drivers:** `mpm gridvel <name> box <lo> <hi> affine <L> [ref <x>]`, changed
  per step with `mpm gridvel <name> setAffine <L>`; `mpm gridvel … profile <t0> <Tr> <t1>
  translate <axis>` (a trapezoidal velocity, both bounds moving). Mode F now applies `gridvel`
  and refuses `mpm pforce` (#1109: both were silently ignored before); `mpm condForce` is mode
  A's and is refused in mode F. On a B-spline grid a moving box is checked each step: the step at
  which it reaches the ghost layer behind a `fix` face is refused (B10).
- **Local damping is not Galilean invariant** (`mpm damping local` acts on the absolute grid
  velocity): turn it off for a rigid translation and back on to settle.

## Two-phase particles — saturated Biot u-p (WP-232l)

Interface `tests/mpm_twophase_interface.md` (P1–P13, §7a, §7b, §7c); results [[232l_mpm_two_phase]].

```tcl
mpm twoPhase <body> -Kf <Kf> -poro <n> -rhoF <rhoF> -perm <k1> <k2> [<k3>] \
    [-alpha <a>] [-Ks <Ks>] [-stab off | auto <a0> | <alpha>] [-fluidBody <bx> <by> [<bz>]]
```

- **The flags are LadrunoUP's.** `-perm` is k̄ = k_hydraulic/γ_w per axis. The body's `-rho` is the
  saturated mixture density. The default stabilization is `auto 0.25`, with α_st = α₀h²/(K_d + 4G_d/3)
  and h the grid spacing. `-dynSeepage on` is refused.
- **Signs.** p is compression-positive excess pore pressure. A particle's `sig` is σ′. The fluid drive
  is λ·(`-fluidBody`, by default `mpm gravity`).
- **Route by mode.**
  - **Route E, mode A (`mpm mode implicit`):** a monolithic u-p step. Every grid node carries
    ndf = ndm + 1 in a two-phase epoch. **Newmark only**: HHT and GeneralizedAlpha are refused, as
    are `quasiStatic` (no S or Qᵀ terms, so it cannot consolidate) and symmetric-storage systems.
  - **Route S, mode F (`mpm mode explicit`):** an explicit skeleton with an implicit pressure solve
    (Jacobi-PCG) and a velocity correction. Two-phase particles take the skeleton strain from the
    corrected grid velocity v₁ (USL, §7c). Particles carry p only.
- **Conditions.**
  - `mpm pfix <face> <p0>`: a drained face. On B-spline faces only p0 = 0 is accepted.
  - `mpm plate <face> <F>`: a rigid, frictionless plate, one normal unknown. F is the total force,
    positive into the body, × λ.
  - `mpm traction <name> <face> <lo..> <hi..> <t..>`: consistent nodal forces, GIMP / linear only.
  - A face without `pfix` is impervious.
- **Cutoffs.** `mpm pressureCutoff <c>` (default 5e-2) sets p-activity: a node's p is active only if
  its two-phase nodal mass is above c·m_max and some two-phase particle has S^p > 0 there.
  `mpm massCutoff` is the u side's cutoff.
  - Mode A refuses `massCutoff` > `pressureCutoff` with a two-phase body.
  - In a route-S step of a two-phase model, a u node below max(massCutoff, pressureCutoff)·m_max is
    dropped, and each particle's u weights are renormalised over its remaining nodes (mode A's
    `renormalize`), so a free block translates exactly (T13).
- **Route S's time step.** The step must satisfy dt ≤ s·Δt_d, where Δt_d is the drained pencil
  (`mpm criticalTimeStep drainedPencil`) and s is `mpm safety`.
  - s defaults to 0.5. Values in (0.5, 0.7] run with a warning once per epoch, and s > 0.7 is
    refused (T7's measured USL limit: gimp 0.7×, bspline2 0.8×).
  - Δt_d is computed once per epoch and refreshed when a particle changes cell, on a stage walk
    or a parameter update, and at least every 100 steps.
  - `mpm criticalTimeStep undrainedPencil` reports the undrained limit. Route S is not bound by it,
    but at the default cap it runs below it anyway: on the gated soil (T2's), Δt_d / Δt_u = 1.414,
    so s 0.5 admits 0.707 × the undrained pencil.
  - `mpm debug routeSLimit off` lifts the limit and its cap (debug only; T7's scan uses it).
- **Queries and state.**
  - `mpm particleData p | pd | pdd`.
  - `mpm particlePressure <p..> [<pd..> <pdd..>]` and `mpm porePressure <body> uniform <p>` set the
    initial field; setting p by command re-arms A′.
  - `mpm info twoPhase <body>` returns [route, constants, α_st, pressureCutoff, max |π_⊥|].
  - `mpm info gridPressure` returns [id, p, p̂] for route S's last step.
  - `mpm kinematics small | hypo` (default `hypo`); `mpm pressureBlend <α_p>` (default 1).
  - Output: the `PORE_PRESSURE` particle field (additive). Restart: mode A pack 5. Mode F keeps
    refusing a database save with particles.
- **Serial.** No new threaded code. Mode A prints the existing serial warning, and the route-S
  pressure assembly and PCG are serial.
- **Refused.**
  - `mpm antiLocking fbar` with a two-phase body.
  - A non-zero `mpm stress` offset on a two-phase particle (stage under gravity instead).
  - Porous materials at particles (PDMY, PIMY, FluidSolidPorous).
  - Walls in route S; restart packs < 5 with a two-phase body.
  - A singular pressure projection, for example one particle per cell on B-spline.
- **Known limits.**
  - **A7:** route S refuses fix 2's term. An impervious B-spline face under a normal fluid drive,
    such as gravity over a fixed base, needs mode A, gimp / linear, or a drained face.
  - **T7 is a recorded miss:** with the USL skeleton, route S is stable only to 0.7× (gimp) / 0.8×
    (bspline2) the drained pencil, hence the cap above.
  - **No drainage switch:** a face is drained or impervious for the whole epoch.
  - Porosity is constant (no J-update).
- **The particle-pressure motion residual.** Particles carry p with FLIP increments on a consistent
  L2 projection (P5a). When supports move, the part of the particle field outside the new supports'
  range stays on the particles and is never damped. The grid field and the field at the particles
  are accurate; the particle's own p can differ by a small residual (6e-6 γ_wH on the 0.044 m
  settlement patch column). Read the field at the particles when that matters.

## Threads in mode F (WP-253)

`ladrunoThreads <n>` (Python `ops.ladrunoThreads(n)`, or `LADRUNO_THREADS`) now applies to mode F's
explicit step. It is the same knob as the element loops (WP-107 / WP-192); the default is 1, and at 1
thread the step runs the serial code unchanged.

- **What threads.** The transfer loops always thread: the support search, P2G, the grid update, G2P,
  the second projection and the particle kinematics (dF, the polar kernel, the stress push), plus the
  admissibility check. The material calls (`setTrialStrain`, `commitState`) and the critical-Δt
  scan (`getInitialTangent`) thread only when every particle material is on the allowlist for
  both `ladrunoThreadSafeUpdate()` and `ladrunoThreadSafeForm()`. Today, among the materials MPM
  admits, that is the 3D `ElasticIsotropic` alone. The 2D plane-strain elastic material is not
  Form-safe (class-wide `sigma` / `D` statics), and DruckerPrager, LadrunoJ2, SANISAND and HS are not
  on the list. For those, the material calls run serially in particle order (the exact serial
  call sequence) between a parallel kinematics pass and a parallel stress push.
  `mpm debug threads` returns the last step's [transfer-loop, material-loop] thread counts.
- **Same answer at every thread count, bit for bit.** P2G and the second projection are computed
  per node as a gather over the node's contributions in ascending (particle, entry) order, which
  is the serial loop's summation order. The slaved ghosts fold in the serial order. Reductions are
  max or lowest-id. Every particle array after any number of steps equals the 1-thread run exactly
  (`tests/test_ladruno_mpm_wp253_threads.py`).
- **Failures.** They are the serial ones: the lowest failing particle is named with the serial
  message (A7, kinematic or material refusal), the materials touched are reverted, and a commit
  refusal still aborts the run (W17).
- **Not threaded.** Route S (two-phase `advanceUP`) and mode A are unchanged by WP-253; mode A's
  own threading is the solver's and WP-192's.
- **Profiler.** The coarse profile splits `mpm.advance` into `mpm.f.support`, `mpm.f.transpose`
  (threaded runs only: the node-to-entry transpose the gather needs, part of P2G's cost),
  `mpm.f.p2g`, `mpm.f.grid`, `mpm.f.g2p`, `mpm.f.project2`, `mpm.f.update` and `mpm.f.commit`.
  The step itself has `mpm.f.admissible` and `mpm.f.dtcrit`. The phase walls are wall time on the
  calling thread. With `-deep` and threaded material calls, the `mpm.f.material` elem_by_type rows
  sum each call's time on its own thread, so they show summed thread time, not wall time.

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
  that particle's own strain path. Thread-count determinism of mode F is
  `tests/test_ladruno_mpm_wp253_threads.py` (WP-253).
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
