# The MPM driver — the frozen design, for implementers

Condensed from the ADR (§4–§7). The ADR is the authority; anchors are on `origin/ladruno` @
`aa6229bf9` — line numbers drift, so re-open them before relying on one.

## Settled (do not re-litigate)

Particles own all history and the grid owns none; the grid is reset every step. Solve over active
nodes only. Implicit MPM is Newton on grid increments with particles frozen at step-start
positions and G2P at convergence; quasi-static is the implicit loop minus inertia. Linear hats are
not production grade (cell crossing): GIMP at least. Contact between bodies is a multi-velocity
field with Coulomb friction.

## Command surface (names provisional until the 232c PR and its guide fix them)

One `LADRUNO_COMMAND` row for `mpm` in `SRC/interpreter/LadrunoCommandTable.h`; it feeds classic
Tcl, DL Tcl and Python (`ci/check_ladruno_commands.py` fails anything registered elsewhere).
Sub-commands: `mpm grid`, `mpm body`, `mpm particles -file {path}`, grid conditions (`mpm fix`,
planes, the frictional wall), `mpm damping local <α>`, `mpm criticalTimeStep`,
`mpm analyze <nSteps> <dt>`, `mpm reset`, `mpm output` (cadence, fields). Classic Tcl keeps its
analysis in different state from the `cmds` singleton, so `mpm analyze` needs two bridges, as
`LadrunoStaggeredAnalyze` has.

## Mode F — the explicit engine (slice 1)

**Host:** a `LoadPattern` subclass (provisional `LadrunoMPMPattern`; a class tag from the PATTERN
registry's fork band — `PATTERN_TAG_LadrunoPorousOverlay` is 33022), after ADR-73's
`LadrunoPorousOverlay`. Not an element: an element with no DOFs makes `FE_Element` call `exit(-1)`
at the first analysis set-up (W15). The host gives, with no Domain edit: deletion by `wipe`
(`Domain::clearAll`, `Domain.cpp:1114`) and `remove loadPattern`; the `parameter … pattern <tag>`
target (`OpenSeesParameterCommands.cpp:156-176`); its own `TimeSeries`. It does not get a
`revertToStart` (only a re-`applyLoad`, `Domain.cpp:2503`), hence `mpm reset`. Outside
`mpm analyze` its `applyLoad` does nothing.

**Each step, in order** — everything a particle carries is trial until step 7:

1. Refuse Δt above the critical-Δt estimate (cell size over the largest particle wave speed plus
   particle speed, from each admitted material's current tangent) times a safety factor.
2. `Domain::applyLoad(t + Δt)`: sets the domain time and `ops_Dt` (`Domain.cpp:2170`) and calls the
   host's `applyLoad`, which reads its time-series factor. Domain loads, `fix` and `sp` never
   reach the engine grid; grid conditions are the engine's own.
3. Grid reset; P2G of mass, momentum, external and internal force over active nodes; drop nodes
   below the mass cutoff.
4. Grid momentum update; grid conditions.
5. MUSL: particle velocity from the grid (PIC / FLIP blend); a second projection of momentum; grid
   velocities; ΔF and rotation through the new kernel entry, with shape functions at the
   step-start positions; `ops_TheActiveElement = 0` on every worker thread; the constitutive
   update. Threads only for allowlisted materials; identical results at any thread count.
6. Commit particle materials. If `ladrunoPendingCommitRefusals() > 0`
   (`LadrunoMaterialStatus.h:162-165`): revert the trial states, clear the counter, fail the step.
7. G2P kinematics: position, velocity, F and volume take their trial values.
8. `Domain::commit()`: FE nodes and elements (none free in v1), committed time, commit tag, then
   recorders — which therefore see positions and stresses from the same step.
9. Replicate the profiler rows: `OPS_PROFILE_SCOPE("step")` and `OPS_PROFILE_STEP`
   (`DirectIntegrationAnalysis.cpp:194`, `:273`).

At the end of `mpm analyze`, call `Domain::flushRecorders()`. **v1 refuses:** free FE DOFs, any
non-admitted material, and database save or restore while particles exist.

## Mode A — particle-as-element (phase 1, not slice 1)

Particle element `LadrunoMPMPoint` (provisional): one integration point, one `NDMaterial` copy;
its nodes are the support of the cell it occupies. Every step: re-point connectivity in place
(the element refetches its `Node*`), raise `Domain::domainChange()` **unconditionally**
(`analyzeStep` re-seeds only when the stamp moved, `DirectIntegrationAnalysis.cpp:207-217`), and
let the analysis rebuild. Transient: `analyzeStep(dt)` (never `analyze`, which sub-steps).
Quasi-static: `StaticAnalysis::analyze(1, false)`, then flush recorders yourself.

- Integrators, checked by `getClassTag()`: Newmark, HHT and Generalized-α only. `LadrunoHHT`
  (33013) and `LadrunoGeneralizedAlpha` (33014) stay refused until a test admits them;
  `VariableTransient` is refused (detect it with `dynamic_cast`).
- G2P in `commitState()`, which must not fail (`Domain::commit` drops its return code). A material
  refusal fails the step through the WP-99 counter; the driver then restores every particle's
  previous kinematics before retrying. Keep xₙ and Fₙ until the step is accepted.
- Stay off the WP-192 threaded commit path (do not opt into `ladrunoThreadSafeForm`).
- Output only through the particle group: `LadrunoRecorder` does not follow in-place re-pointing
  (its topology fingerprint is unchanged).

## Particle state and the kernel entry

State: x, V, m, v, a (implicit), F, R, the accumulated strain feed, σ, GIMP half-lengths, one owned
`NDMaterial` copy, a stable particle id, body id, material tag, flags.

The new hypo-kernel entry goes in `SRC/element/solidTransformation/LadrunoHypoKernel.h` beside
`hypoIncrement` (which needs nodal displacements a reset grid does not have). Its expected shape
is (Fₙ, F₁) → (strain increment in the unrotated frame, R₁, J₁), reusing `polar3`, the
Hughes–Winget midpoint and `pushStress6`. 232b freezes the exact signature; the plane-strain
embedding is new.

## Materials (v1)

| Material | Status | Notes |
|---|---|---|
| `ElasticIsotropic` | admitted | on the WP-107 allowlist (3D, PlaneStrain2D) |
| `LadrunoJ2` | admitted | serial; refuse `-autoRegularization` at particles — its only length read (`LadrunoJ2.cpp:352-356`) |
| `DruckerPrager` (UW, built after `31322a47a`) | admitted | serial; stage switch `materialState` (`DruckerPrager.cpp:1045`); `updateMaterialStage` returns −1 (`:1063-1064`) |
| `J2Plasticity` | refused | process statics (`J2Plasticity.h:137-139`) |
| concrete and steel length consumers | refused | they read the element-length global |
| SANISAND, HS, NorSand, PDMY, PM4Sand | refused in v1 | later phase, starting with the SANISAND cell-crossing meter |

Stage path in mode F: `parameter <t> pattern <mpmTag> …` then `updateParameter`; the host forwards
to its particles' materials. `updateMaterialStage` never reaches mode F — say so in the guide.

## The two contracts

**Contract 1 — the particle set (apeGmsh → fork), HDF5.** Attributes: format version, provenance,
units, `ndm`, plane-strain thickness, grid (origin, cell size, extents). Datasets of length N:
particle id (int64, never reused), position (float64), volume (positive), mass or density (one,
declared), GIMP half-lengths, initial velocity, initial stress (Voigt order and sign stated), body
id, material tag (resolved by apeGmsh's emitter from a label — charter P4), flags. Chunked and
compressed. Read by `mpm particles -file {…}`, with the path brace-quoted in Tcl.

**Contract 2 — particle results (fork → apeGmsh): a new group inside `.ladruno`.** Written by
`LadrunoRecorder` or a writer sharing its file; `FORMAT_VERSION` goes from 1
(`LadrunoRecorder.cpp:989`) to 2, with a format-document PR. Per output: time, step, particle ids,
positions (float64), velocities, stress, named state variables; grid metadata once per stage.
Never every step — output every N steps or every Δt_out, with field selection. apeGmsh reads it
in `src/apeGmsh/results/readers/_ladruno.py` and `apeGmshViewer/src/reader/`.

## Walls — one line each

| # | Wall | Anchor |
|---|---|---|
| W1 | transient integrators read the Nodes only in `domainChanged()` | `Newmark.cpp:206-208`, `:398-414` |
| W2 | one ConstraintHandler per model; the contact handler is Plain-only | `LadrunoContactHandler.h:22-29` |
| W3 | a DOF_Group for every Domain node | `PlainHandler.cpp:123-132` |
| W4 | stock recorders cache pointers; `LadrunoRecorder` re-maps only on a topology-fingerprint change | `LadrunoRecorder.cpp:664-694` |
| W5 | removal cleans nothing | `Domain.cpp:1198-1241` |
| W6 | `PartitionedDomain`: repartition stub, dead subdomain loop | `PartitionedDomain.cpp:1147-1156` |
| W8 | the hypo kernel entry is reference-configuration and 3D-only | `LadrunoHypoKernel.h:166-176`; `OPS_LadrunoUP.cpp:511-514` |
| W9 | `ops_TheActiveElement` is written by the Element constructor, `addElement` and the update loops | `Element.cpp:131`; `Domain.cpp:493` |
| W10 | `analyze()` sub-steps when `numSubLevels != 0` | `DirectIntegrationAnalysis.cpp:176-177` |
| W11 | recorders run last in `Domain::commit`; element `commitState` return codes are dropped | `Domain.cpp:2294-2423` |
| W12 | stage and parameter routing walk Domain elements only | `MaterialStageParameter.cpp:60-81` |
| W13 | the GP layout of `.ladruno` is mesh-bound | `LadrunoRecorder.cpp` |
| W14 | no soil model on the WP-107 thread allowlist | `ManzariDafalias.cpp:1235` |
| W15 | an element with no DOFs → `exit(-1)` at analysis set-up | `FE_Element.cpp:126-129`; `PlainHandler.cpp:291-327` |
| W16 | `StaticAnalysis` has only `analyze(numSteps, flush)` | `StaticAnalysis.h:66` |
