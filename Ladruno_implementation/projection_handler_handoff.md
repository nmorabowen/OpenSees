---
title: "LadrunoProjectionHandler (ADR-30) — handoff / v1 + prescribed motion + Noh–Bathe adoption complete"
project: Ladruno
type: handoff
status: v1 (P0–P3) + P4a recorder + P4b/P4c prescribed-motion + P5 ExplicitBathe/LNVD adoption + P6 condition/staleness guards (closes O1) — P0–P5 merged to ladruno (P5 #334), P6 pending, 2026-06-21
related:
  - "[[30_ladruno_explicit_constraint_projection_adr]]"   # the spec
  - "[[_adr30_p1_design]]"                                 # P1 pseudocode + Gate-A resolutions
  - "[[_adr30_p5_design]]"                                 # P5 ExplicitBathe/LNVD adoption design
  - "[[_adr30_loop_state]]"                                # the full phase-by-phase build log
  - "[[05_robust_central_difference]]"                     # CentralDifferenceLadruno (host integrator)
tags: [adr, constraints, explicit-dynamics, projection, handoff, ls-dyna, abaqus]
---

# LadrunoProjectionHandler — handoff

**What it is.** A constraint handler (`HANDLER_TAG 33001`, the fork's first) that enforces
`MP_Constraint`s and `EQ_Constraint`s in **explicit dynamics** by **M-orthogonal projection
of the solved accelerations**, instead of penalty (eats `dt_cr`) or `Transformation`+`Diagonal`
(silently drops the condensed off-diagonal mass — falsified in P0). It is the explicit
constraint enforcement OpenSees lacked: **exact, momentum-conserving, Δt-neutral, and it keeps
the O(n) `system Diagonal` recipe**. = LS-DYNA `*CONSTRAINED_LINEAR` / Theory §28 nodal-set
`a = ΣMᵢaᵢ/ΣMᵢ`; the kinematic-contact half of Abaqus/Explicit.

**Status: v1 + PRESCRIBED MOTION COMPLETE.** P0 #300 · P1 #305 · P2 #307 · P3 #312 · P4a
tie-force recorder #317 · **P4b prescribed motion (non-homog SP / imposedMotion on free DOFs)
#323** · **P4c prescribed MASTER (drives constrained slaves) #329** — all merged to `ladruno`.

---

## 1. How to use it (user-facing)

```python
ops.constraints("LadrunoProjection")          # <-verbose> <-projectICs> <-icTol $tol>
ops.system("Diagonal")                          # REQUIRED (the explicit lumped recipe)
ops.integrator("CentralDifferenceLadruno")      # or ExplicitBathe / ExplicitBatheLNVD and their
                                                #   SMS/Consistent variants (P5); the whole explicit family
ops.algorithm("Linear")                         # one solve/step
# ... equalDOF / rigidLink / rigidDiaphragm / equationConstraint as usual ...
ops.analyze(nsteps, dt)

f = ops.ladrunoProjectionTieForce(nodeTag, dof) # the constraint tie force M(a_raw - a_proj)

# native HDF5 recording (P4a): writes the CONSTRAINT_TIE_FORCE field per node per step
ops.recorder("ladruno", "out.ladruno", "-N", "constraintTieForce")

# prescribed motion (P4b/P4c): base excitation + supports that drive ties, all under projection
ops.sp(supportNode, dof, value, "-const")                 # non-homogeneous SP (prescribed disp)
ops.pattern("MultipleSupport", p); ops.groundMotion(g, "Plain", "-accel", ts)
ops.imposedMotion(supportNode, dof, g)                     # ground-motion support (disp+vel+accel)
```

**Supported constraints:** `equalDOF`, `equalDOF_Mixed`, `rigidLink -bar`/`-beam`,
`rigidDiaphragm`, `equationConstraint` (EQ). Homogeneous `fix` is handled (equation exclusion).
**Prescribed motion (P4b/P4c):** non-homogeneous `SP` / `imposedMotion` works on a free support
DOF (Tier 1, P4b — the structure feels it through `F_int`/damping) AND as a constraint **master**
that drives its slaves (Tier 2, P4c — the slave is KINEMATICALLY imposed `u_c=ΣC_k u_{m_k}+delta`,
exact). A plain `SP` supplies disp only (vel/accel = 0, same as `Transformation`); `imposedMotion`
supplies all three. A DOF that is BOTH prescribed and a constraint slave, or a slave tied to BOTH
a free and a prescribed master (mixed), or a group with zero free equations, is refused (named).

**HARD REQUIREMENTS / gotchas (all enforced with named errors — see LEDGER_quirks):**
- **`system Diagonal`** only (the projector reads the assembled lumped-mass diagonal). A
  non-Diagonal SOE is refused.
- **Every TIED DOF needs nonzero lumped mass** — the projection keeps slave DOFs in the
  equation set (it does NOT eliminate them like `Transformation`). `rigidLink -beam` and 3D
  `rigidDiaphragm` tie the slave's **perpendicular rotation**, so those rotational DOFs need a
  small rotational mass (~0.01–0.1 % of the node's translational mass). A massless tied DOF is
  refused with a recipe-bearing message.
- **Do not SP-fix a DOF a constraint controls** (e.g. a diaphragm slave's `rz`) → SP-on-slave
  refusal.
- **Initial state must be on the constraint manifold** (committed `u`,`v` satisfy the MP
  increment relation `u_c−Uc0 = C(u_r−Ur0)`). Off-manifold ICs ABORT the analysis; pass
  `-projectICs` to snap them compliant instead.
- **Refused (named errors):** MP chains (a slave retained elsewhere), a DOF constrained twice,
  partition-crossing / missing-node constraints (OpenSeesSP/MP — v1 is partition-interior only).
- **Frozen small-rotation `Ccr`** (same limitation `Transformation` carries): the transport
  lever arms are captured once; not valid for accumulated rotation ≳0.1 rad. Use the RBE2 element
  (`LadrunoKinematicCoupling`) or `Transformation` for finite-rotation rigid offsets.

---

## 2. Architecture (where the pieces live)

Fork-owned, `SRC/analysis/handler/`:
- **`LadrunoConstraintProjector.{h,cpp}`** — the math core. Per connected-component group it
  stores `L=[I;C]` (rows = all group DOFs retained-first; constrained rows = `C`) and computes
  `a_proj = L(LᵀML)⁻¹LᵀM a_raw`. `buildMass()` reads `diag(M)` from the Diagonal SOE; `project()`
  is the in-place per-group solve + tie-force cache `f = M(a_raw−a_proj)`; `checkIC`/`snapICs`/
  `enforceIC` handle the IC manifold check.
- **`LadrunoProjectionHandler.{h,cpp}`** — Plain-style assembly that does NOT eliminate slave
  DOFs; union-find grouping over `(node,dof)` vertices from `getMPs()` + `getEQs()`; the
  diagnostics (chain / double-constraint / SP-on-slave / zero-free / massless / consistent-mass /
  partition-crossing); `doneNumberingDOF()` freezes equation indices and pushes the projector;
  `getTieForce(node,dof)`.
- **`LadrunoProjectionConsumer.h`** — a 1-method interface (`setConstraintProjector`). The handler
  `dynamic_cast`s the Integrator to this, so it never names a concrete integrator → ExplicitBathe/
  LNVD adopt projection in P4 by implementing this one method.

CDL hook (fork-owned, `SRC/analysis/integrator/CentralDifferenceLadruno.{h,cpp}`): implements the
consumer; **buildMass between formTangent and solve** (the DiagonalDirectSolver overwrites the
diagonal with 1/m on the factor pass); **project a0** in the starter before seeding `v_{−1/2}`;
**project a_{n+1}** in `update()` — the `*Aprev = Aproj` write is what carries manifold invariance.
P4a adds a tie-force **scatter in `commit()`**: `M(a_raw−a_proj)` → `Node::setProjectionTieForce`
(before `commitDomain()`), guarded by `theProjector && isMassReady()` so a non-projection run pays
nothing.

Recorder (fork-owned): `Node::projTieForce` slot (vanilla Node, P4a); the `CONSTRAINT_TIE_FORCE`
node-result = `Ladruno_Types.h` enum + `Ladruno_NodeResults` `ConstraintTieForceSource` +
`LadrunoRecorder.cpp` keyword/switch.

Vanilla seams (all `// Ladruno`, in LEDGER_vanilla): `classTags.h` (33001); `DiagonalSOE::getDiagonalA`;
the `constraints LadrunoProjection` branch in `OpenSeesCommands.cpp` + `tcl/commands.cpp`; the broker
case; the `ladrunoProjectionTieForce` command (OpenSeesCommands.h + OpenSeesOutputCommands.cpp +
PythonWrapper.cpp + TclWrapper.cpp); `DirectIntegrationAnalysis` + `TransientDomainDecomposition
Analysis` error-contract checks; **`Domain::clearAll()` EQ-leak fix** (upstreamable).

Command: `SRC/interpreter/OpenSeesOutputCommands.cpp::OPS_LadrunoProjectionTieForce`.

---

## 3. Why it is correct (the validation trail)

- **P0 (#300)** falsified the premise: `Transformation`+`Diagonal` silently drops the condensed
  off-diagonal mass (rz-coupled `uy` → 0), matched to <1% by an OpenSees-free analytic modal
  solution. Also pinned §2.5: a massless DOF makes the SOE untrustworthy (Diagonal aborts opaquely;
  Full/Band swallow the singular factor and report success) → the handler MUST scan up front.
- **Gate-A** (pre-code math panel): the projection algebra (M-orthogonality, momentum, equalDOF→
  mass-weighted average, Δt-neutrality via the Galerkin-reduced pencil) is SOUND and **general-C**
  — which is why P2 needed no new projector code.
- **Gate-B** (P1 code) + **Gate-C** (P2 transport) + **Gate-D** (P3) + a **general pre-merge review**:
  implementation verified exact (tie to machine zero, momentum to 4e-16, transport gap ~1e-12).
- **Tests:** `tests/test_adr30_projection_p{0,1,2,3,4}.py` (~30 cases): T1 momentum, T2 vs
  `Transformation` ref, T3 diaphragm, T4 Δt-neutrality, T5 conflict battery, T6fix (fixes the P0
  mass-drop), T6b analytic anchor, T7 rigidLink free-vib, T8 tie-force, T9 energy closure, EQ,
  IC-violation, UAF, missing-node, P4a recorder HDF5 readback. Full Zone-A green (942).

**Bugs this work caught and fixed** (the high-value byproducts): a latent upstream silent-swallow
of handler/integrator error returns in `DirectIntegrationAnalysis` (now honored across the transient
drivers), and a real upstream `Domain::clearAll()` bug that **leaked `EQ_Constraint`s across `wipe`**
(invisible until a handler read `getEQs()`).

---

## 4. P4 — status & remaining backlog

**P4a — native tie-force recorder: DONE (#317, merge 61341837).** Chosen route = a **dedicated
nodal buffer** (not the reaction slot, to keep it separate from reactions): `Node` gets a lazily-
allocated `projTieForce` slot (`get/setProjectionTieForce`, mirrors `reaction`, transient/not-
serialized); `CentralDifferenceLadruno::commit()` scatters `M(a_raw−a_proj)` onto the nodes before
`commitDomain()`; a `ConstraintTieForceSource` (`NodalResultType::ConstraintTieForce`, keyword
`constraintTieForce`/`tieForce`) writes the `CONSTRAINT_TIE_FORCE` field to the `.ladruno` HDF5.
Test `tests/test_adr30_projection_p4.py` (h5py readback `DATA == ±F·m₂/M`). Two adversarial reviews
(focused + general) both SAFE-TO-MERGE — Node core has zero blast radius, the scatter is free when
projection is inactive. Optional test nits left: a plain-CDL-no-projection zeros-recording test;
`h5py` is `importorskip` (silently skips on a box without it).

**P4b — prescribed motion on free DOFs (Tier 1): DONE (#323, merge 1d445ddc0).** A non-homogeneous
`SP`/`imposedMotion` DOF stays SP-excluded (`eqn=-1`); a new `LadrunoProjectionHandler::applyLoad()`
imposes its prescribed DISPLACEMENT each step (`setTrialDisp(getValue()+getInitialValue(),dof)`,
mirrors `TransformationDOF_Group::enforceSPs`); `ImposedMotionSP` supplies vel/accel. The integrator
is UNCHANGED — free DOFs feel the support via `F_int`/damping (`setNode*` skips `eqn<0`). Handler-only,
ZERO vanilla. A prescribed DOF used as a constraint master/slave was refused here (master refusal
lifted in P4c).

**P4c — prescribed MASTER (Tier 2): DONE (#329, merge d9108c84).** A slave driven PURELY by prescribed
master(s) is KINEMATICALLY imposed: excluded from the equation set and `u/v/a` set from the masters
each step in `applyLoad()` (`classifyDerivedSlaves()` detects+excludes+refuses-mixed; `doneNumberingDOF()`
drops prescribed masters + skips derived slaves). EXACT (drift 0, matches `Transformation`).
**Key lesson (LEDGER_quirks):** ADR §2.4's known-RHS *acceleration* projection was implemented first
then ABANDONED — a prescribed master's disp is externally imposed while the slave leap-frogs, so
accel-only projection leaves a NON-converging O(1e-3) displacement-tie drift. A kinematic tie to an
externally-imposed DOF must be imposed kinematically. Mixed (free+prescribed master on one slave) and
zero-free-DOF groups are refused (named). Handler-only, projector UNTOUCHED, ZERO vanilla.

**P5 — ExplicitBathe / LNVD adoption: DONE.** The `LadrunoProjectionConsumer` is now implemented by the
**Noh–Bathe** family, so `constraints LadrunoProjection` works under it too (previously CDL-only). Only TWO
base classes were edited — `ExplicitBathe` and `ExplicitBatheLNVD` — and the 4 SMS/Consistent subclasses
(33009–33012) inherit it (all chain their `domainChanged()` to the base). The entire CentralDifference
family (incl. SMS 33007/33008) already had projection via inheritance. Each base reads `diag(M)` once per
stage (a one-shot `formTangent`+`buildMass` in `newStep`, before the algorithm's solve factors the SOE),
projects the committed a0 at the starter, and projects **both** Noh–Bathe sub-step accelerations in
`update()` (after any consistent-mass `refineAccel`, mirroring the CDL contract); the load-bearing carry is
`A_t = A_tdt` (projected) at commit. The `commit()` tie-force scatter is replicated so the P4a
`constraintTieForce` recorder works under these integrators. Handler error message broadened to name the
family. Design: `_adr30_p5_design.md`. Tests `tests/test_adr30_projection_p5.py` (TP5-1..8: tie-exactness,
projection==Transformation under EB and LNVD, FLAC-damped tie, the two inherited SMS subclasses, a
no-projection regression, recorder readback). ZERO vanilla (all six are fork-authored integrators).

**P6 — condition gate + frozen-Ccr staleness guard: DONE (closes ADR O1).** Two warn/refuse-only
robustness guards (NO change to the projection math). **(3a)** `LadrunoConstraintProjector::buildMass()`
replaced its single exact-pivot test-solve with a real **condition-number gate**: a self-contained
symmetric **Jacobi** eigensolve (`ladrunoSymEigJacobi`, no LAPACK dependency — sidesteps the `dsygv_`
gap) estimates `cond = λmax/λmin` of each SPD `LtML`; **refuse** above `1e12` (also catches the exact
singular `λmin≤0`), **warn** above `1e8`. Catches the near-singular set (barely-dependent retained
direction, disparate tied masses, near-redundant constraint) the old check missed. **(3b)** Handler
**frozen-Ccr staleness warn**: at transport, `flagRotMonitor()` records masters whose Ccr couples a
**rotational** master DOF into a **translational** slave DOF (the rigidLink-beam / rigidDiaphragm lever
arm — equalDOF rotation-ties and bar links are NOT flagged), capturing `θ0`; `applyLoad()` warns ONCE
per master when the rotation drift crosses `0.1 rad`, pointing to RBE2 / `Transformation` for finite
rotation. Design: `_adr30_p6_design.md`. Tests `tests/test_adr30_projection_p6.py` (TP6-1..5: gate
passes a healthy group / refuses cond~1e13 / warns-but-runs cond~1e10; staleness NOTE fires on a
>0.1 rad rigidLink-beam master and not on a small rotation). ZERO vanilla (projector + handler are fork).

**Remaining backlog (priority order, all independent, none required for the core feature):**
1. **MP-chain composition** (currently refused): substitute `C` matrices (Abaqus-style).
2. **RBE2/RBE3 eliminable-block routing** through the projector (retire bipenalty where eliminable).
3. **SOE-cooperative massless-slave elimination** — would relax the "every tied DOF needs mass"
   requirement (lets diaphragm slaves keep zero rotational mass like `Transformation`).
4. **Prescribed-motion deferred sub-cases:** a *plain* `SP` master supplies disp-only (vel/accel=0,
   documented); a slave tied to BOTH free and prescribed masters (mixed) is refused rather than
   solved — a full mixed treatment would need a genuine displacement-level projection.

## 5. Pointers
- Spec + decisions: `30_ladruno_explicit_constraint_projection_adr.md` (P0–P3 marked DONE inline;
  §2.4 prescribed motion). Prescribed-motion design: `_adr30_p4b_design.md` (Tier 1) and
  `_adr30_p4c_design.md` (Tier 2 — includes the abandoned known-RHS approach + the kinematic pivot).
- The full phase-by-phase build log (every gate, finding, fix): `_adr30_loop_state.md`.
- Gotchas: `LEDGER_quirks.md` (Diagonal off-diagonal drop · tied-DOF-needs-mass · frozen-Ccr ·
  near-singular LᵀML · EQ-clearAll leak). Provenance: `LEDGER_vanilla_files.md`,
  `LEDGER_implementations.md` (LadrunoProjectionHandler row).
