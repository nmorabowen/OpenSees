---
title: "LadrunoProjectionHandler (ADR-30) — handoff / v1 complete"
project: Ladruno
type: handoff
status: v1 SHIPPED (P0–P3 merged to ladruno, 2026-06-20)
related:
  - "[[30_ladruno_explicit_constraint_projection_adr]]"   # the spec
  - "[[_adr30_p1_design]]"                                 # P1 pseudocode + Gate-A resolutions
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

**Status: v1 COMPLETE.** P0 #300 · P1 #305 · P2 #307 · P3 #312, all merged to `ladruno`.

---

## 1. How to use it (user-facing)

```python
ops.constraints("LadrunoProjection")          # <-verbose> <-projectICs> <-icTol $tol>
ops.system("Diagonal")                          # REQUIRED (the explicit lumped recipe)
ops.integrator("CentralDifferenceLadruno")      # the only projection-aware integrator in v1
ops.algorithm("Linear")                         # one solve/step
# ... equalDOF / rigidLink / rigidDiaphragm / equationConstraint as usual ...
ops.analyze(nsteps, dt)

f = ops.ladrunoProjectionTieForce(nodeTag, dof) # the constraint tie force M(a_raw - a_proj)
```

**Supported constraints (v1):** `equalDOF`, `equalDOF_Mixed`, `rigidLink -bar`/`-beam`,
`rigidDiaphragm`, `equationConstraint` (EQ). Homogeneous `fix` is handled (equation exclusion).

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
- **Tests:** `tests/test_adr30_projection_p{0,1,2,3}.py` (~26 cases): T1 momentum, T2 vs
  `Transformation` ref, T3 diaphragm, T4 Δt-neutrality, T5 conflict battery, T6fix (fixes the P0
  mass-drop), T6b analytic anchor, T7 rigidLink free-vib, T8 tie-force, T9 energy closure, EQ,
  IC-violation, UAF, missing-node. Full Zone-A green (938).

**Bugs this work caught and fixed** (the high-value byproducts): a latent upstream silent-swallow
of handler/integrator error returns in `DirectIntegrationAnalysis` (now honored across the transient
drivers), and a real upstream `Domain::clearAll()` bug that **leaked `EQ_Constraint`s across `wipe`**
(invisible until a handler read `getEQs()`).

---

## 4. Deferred — P4 backlog (priority order)

1. **Native `LadrunoRecorder` tie-force source** — *user asked for this.* OPEN DESIGN DECISION:
   the `LadrunoRecorder` is node-based (iterates Domain nodes) and **cannot reach the handler-owned
   projector**, so recording needs either (a) write the tie force into the node **reaction** slot
   (free query via `nodeReaction` + free recording via the existing reaction source, but it lives in
   the reaction slot), or (b) a **dedicated nodal tie-force buffer** (integrator writes it each step)
   + a new `Ladruno_NodeResults` source. v1 ships the lean query only; pick (a) or (b) for P4.
2. **Prescribed-motion overwrite** (non-homogeneous SP / `imposedMotion`): overwrite `a` on those
   DOFs before projecting the rest.
3. **MP-chain composition** (currently refused): substitute `C` matrices (Abaqus-style).
4. **ExplicitBathe / LNVD adoption** — implement `LadrunoProjectionConsumer` + the per-sub-step hooks.
5. **RBE2/RBE3 eliminable-block routing** through the projector (retire bipenalty where eliminable).
6. **Near-singular `LᵀML`** condition-number gate (the rank check catches only an exact zero pivot;
   ADR O1) and a **frozen-`Ccr` runtime staleness guard** (warn at large accumulated rotation).
7. **SOE-cooperative massless-slave elimination** — would relax the "every tied DOF needs mass"
   requirement (lets diaphragm slaves keep zero rotational mass like `Transformation`).

## 5. Pointers
- Spec + decisions: `30_ladruno_explicit_constraint_projection_adr.md` (P0–P3 marked DONE inline).
- The full phase-by-phase build log (every gate, finding, fix): `_adr30_loop_state.md`.
- Gotchas: `LEDGER_quirks.md` (Diagonal off-diagonal drop · tied-DOF-needs-mass · frozen-Ccr ·
  near-singular LᵀML · EQ-clearAll leak). Provenance: `LEDGER_vanilla_files.md`,
  `LEDGER_implementations.md` (LadrunoProjectionHandler row).
