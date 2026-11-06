---
wp: LEGACY
title: "Multi-element Concrete3D softening under implicit Newton cut-crawls even where single elements pass — use -implex (uniform LoadControl) for band/localization r…"
legacy_seq: 154
---
### Multi-element Concrete3D softening under implicit Newton cut-crawls even where single elements pass — use `-implex` (uniform LoadControl) for band/localization runs
- **Bites:** implicit displacement-ramped runs of a MULTI-element `LadrunoConcrete3D` specimen through localization: plain Newton + step-cutting (the recipe that walks a SINGLE element through its limit point) converges only at micro-steps once several elements carry the indefinite softening tangent simultaneously — the G5 coarse band (4 elements!) took ~2400 micro-steps / 13 min wall; refinement makes it worse.
- **Why:** the indefinite/non-symmetric tangent cluster around the band throws the global Newton into cut/recover oscillation; no single step fails permanently, so nothing surfaces except wall time.
- **Workaround/status (2026-07-06, ADR 66 P5.2 G5):** put the MATERIALS in `-implex` and drive with CONSTANT-dlam LoadControl (the uniform-pseudo-time regime IMPL-EX wants; kinematic sp ramp): the SPD-ish secant lets Newton track the full localization at the planned step size (13 min -> 16 s on the G5 coarse mesh; committed states stay implicit-exact). The ADR 66 risk register lists exactly this toolbox row; the `-implex`+DisplacementControl limit-point trap (Concrete3D ledger) does NOT apply because LoadControl dlam is uniform.
**`LysmerTriangle` stage-3 `getResistingForce()` MUTATES state on every call — any recorder that
reads element forces perturbs it (2026-07-05, ADR-69).** At stage 3 ("preserve elastic spring
forces after gravity") `getResistingForce()` executes `internalForces -= springForces` on EVERY
invocation — it is not idempotent. The EnergyBalanceRecorder (v1 AND v2) calls it once per record
per element, so each record subtracts `springForces` again from the member the residual path also
serves (rebuilt only at the next `getResistingForceIncInertia`). Consequences: (a) stage-3 Lysmer
energy readings are untrustworthy; (b) anything else querying element forces between residual
formations (nodal reactions, other recorders) compounds it. The ADR-69 leak publisher deliberately
RECOMPUTES `R_inj = getDamp()*v_gnd` in `commitState` instead of reading the member, so E_inject is
immune. Upstream-origin behavior — left unfixed (vanilla change budget); avoid stage 3 + per-step
force recorders in the same model, or accept the drift.
**Tcl `eleLoad` SILENTLY ACCEPTS unknown `-type` flags (returns TCL_OK, no warning) — a
no-op that looks like success (2026-07-06, ADR-69 P0.5).** The eleLoad handler's tail falls
through to `return 0` when no `-type` branch matches, so `eleLoad -type -fooBarBazNotALoad`
"succeeds". Any deck relying on a loader that is not actually wired (LysmerVelocityLoader was
exactly this for 15 years) runs unloaded with zero diagnostics. Left unfixed (changing the
return could break decks); when a load seems dead, FIRST verify the `-type` string exists in
TclModelBuilder.cpp before debugging the physics.

**Stage-0 `LysmerTriangle` under implicit Newmark realizes only ~HALF the dashpot energy
(DW = 0.50*ULW measured; 2026-07-06, ADR-69 P0.5 F2).** The element's
`getResistingForceIncInertia` uses `0*v_node + gnd_velocity` — the node-velocity damping
force C*v NEVER enters the element residual; under Newmark damping then acts only through the
a1*C term in the effective tangent, giving an energy-inconsistent solve (the recorder's
DW = int v'Cv dt books the full ideal dashpot power and RES exposes the ~0.5*W gap). EXPLICIT
integrators assemble the damping force from getDamp() directly and are consistent. For
implicit absorbing runs prefer ASDAbsorbingBoundary; for Lysmer prefer explicit. NOT a
recorder bug — the recorder is the instrument that surfaced it.

**`UniformExcitation` input work HIDES INSIDE the EnergyBalance recorder's IE column (IE = -DW
exactly, RES accidentally closed; 2026-07-06, ADR-69 P1.6).** Elements that implement
`addInertiaLoadToUnbalance` (FourNodeQuad etc.) store `-M*ug''` in their element load vector
`Q`, and `getResistingForce` returns `K*u - Q` — so the recorder's IE integral
(`int F.v dt`) silently accumulates MINUS the seismic input work. The balance then "closes"
with IE the negative mirror of the genuine absorbed/damped energy and ULW = 0. Consequence
for closure gates: NEVER drive an energy-balance validation model with UniformExcitation —
the input-work pollution drowns whatever leak the gate is trying to isolate (use initial
velocities: no patterns, Q = 0, IE = pure strain energy). Not a recorder bug per se — a
consequence of OpenSees folding element loads into the resisting force.

**`ASDAbsorbingBoundary2D/3D::addInertiaLoadToUnbalance` is a deliberate NO-OP ("we don't
need this!") — free-field columns are NEVER driven by `UniformExcitation` (2026-07-06,
ADR-69 P1.6).** The FF masses (`addMff`) receive no `-M*ug''` effective load, so under
uniform excitation the FF column rides rigidly in relative coordinates: zero strain, zero
`addRffToSoil` transfer, dead lateral boundary. The intended input path for ASD boundaries
is the BOTTOM compliant base (`"B" -fx/-fy` time series); lateral elements take no
time-series args at all (parser rejects them for non-bottom). Also note the lateral
`addClk` dashpot writes only SOIL rows (one-way coupling): the FF column is UNDAMPED unless
element Rayleigh (`addCff`, alphaM) is set — an undamped FF column rings forever.

**openseesmp flat-per-rank nodal-term summation is convention-dependent — the split-mass idiom
sums CORRECTLY, only full-mirror emits double-count (2026-07-06, ADR-69 P2, measured).** The
upstream MPI example declares `mass(4, m, m)` for a shared node on BOTH ranks; the parallel
diagonal assembly SUMS duplicate contributions, so the assembled system has `2m` and each
rank's EnergyBalance recorder books only its own share — the naive cross-rank sum of the nodal
columns equals the serial (2m) reference EXACTLY (gate `energy_v2/p2_mpi_owned_nodes.py`, G2).
The "nodal terms multiply-counted on shared boundary nodes" hazard (ADR-69) applies only to
FULL-MIRROR conventions: PartitionedDomain-style external-node mirrors, or an emitter writing
the full nodal mass/load on every touching rank (which also changes the assembled physics
unless the solver dedups). Consequence: don't "fix" per-rank energy sums blindly — first
determine which convention the model uses; `-ownedNodes <regionTag>` is the dedup tool for
mirror conventions and a no-op burden otherwise. Also note per-rank output files
(`stem.part-<rank>.ext`) are auto-suffixed under a detected MPI launcher since ADR-69 P2 —
ranks racing a single recorder file was the previous (corrupting) behavior.

**Modal-damping energy is published only by integrators using the BASE
`IncrementalIntegrator::commit()` (2026-07-06, ADR-69 P2).** 35 integrators override
`commit()` without chaining (HHT family, `*_TP` explicit): they still APPLY modal forces (via
`TransientIntegrator::formUnbalance` / their own `addModalDampingForce` calls) but never reach
the publish site, so their modal dissipation stays in RES and no `E_modal` column appears
(declare-on-first-publish prevents a silent zero column). Newmark does not override commit and
is fully covered. If you need E_modal under HHT: either chain the override to the base commit
(vanilla edit, ledger it) or accept the documented RES drift.

**Recorders NEVER receive domainChanged() — any recorder caching pointers is a
use-after-free waiting for `remove element`/`remove node` (2026-07-06, ADR-69 P2.1).**
`Domain::removeElement` calls `domainChange()` (sets a flag) but `Domain::record()` invokes
recorders directly with no invalidation hook, and the analysis propagates domain changes only
to its own components (handler/numberer/integrator/algorithm). Worse, a recorder CANNOT poll
`Domain::hasDomainChanged()` — it is STATEFUL (consumes the flag, increments currentGeoTag,
resets graph-built flags) and belongs to the Analysis; calling it from a recorder would eat
the analysis's own change detection. `getDomainChangeFlag()` is a pure read but is usually
already consumed by the time record() runs. The working patterns (both in-tree): (1)
re-resolve entities BY TAG on every emit (LadrunoMonitorRecorder, #489); (2) structural
re-validation per record — compare `getNumElements()/getNumNodes()` (O(1)) against cached
sentinels + verify cached tag→pointer identity before any virtual call through a cached
object, rebuild on mismatch (EnergyBalanceRecorder P2.1). Key ALL membership maps by TAG,
never by pointer — a freed pointer can be REUSED by a new allocation and silently inherit the
old binning. Note the count sentinel alone is blind to remove-then-readd-same-tag (counts
restore); the tag→pointer identity check is what catches it.
