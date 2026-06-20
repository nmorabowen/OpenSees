---
title: "ADR-30 P4b — prescribed-motion overwrite (design / pre-code)"
project: Ladruno
type: design
status: PRE-CODE — design gate pending
related:
  - "[[30_ladruno_explicit_constraint_projection_adr]]"   # §2.4 + §7 P4 are the spec
  - "[[projection_handler_handoff]]"                        # P4 backlog item 1
  - "[[_adr30_loop_state]]"                                 # phase ledger
tags: [adr, constraints, explicit-dynamics, projection, prescribed-motion, imposedMotion, ground-motion]
---

# P4b — prescribed-motion overwrite

**Goal.** Make `imposedMotion` (ground-motion multi-support / uniform excitation) and
non-homogeneous `SP_Constraint` (prescribed displacement) WORK under
`constraints LadrunoProjection`. Today the handler **warns and silently drops** every
non-homogeneous SP (`LadrunoProjectionHandler.cpp:127-131`) and marks all SP DOFs `-1`
(out of the equation set) — so a base-excitation model run with the projection handler
gets the *wrong* answer (the prescribed motion is ignored, the supports act fixed).

ADR §2.4 / §7-P4 locked the mechanism: *"overwrite `a` on those DOFs with the prescribed
acceleration before projecting the rest."*

## 0. Code facts established (verified 2026-06-20)

- `ImposedMotionSP::applyConstraint(time)` sets the node's **trial vel & accel** directly
  from `GroundMotion::getDispVelAccel(time)` (a `Vector(3)` = `[disp,vel,accel]`);
  the **displacement is left to the constraint handler** (`ImposedMotionSP.cpp:127-141`).
- Plain `SP_Constraint` exposes only the prescribed **displacement** via `getValue()`
  (`= valueC`, the reference value × loadFactor) — **no vel/accel** (`SP_Constraint.cpp`).
- `TransformationConstraintHandler` imposes the prescribed disp by
  `myNode->setTrialDisp(getValue()+getInitialValue(), dof)` in `enforceSPs`
  (`TransformationDOF_Group.cpp:1048-1059`), called from its `applyLoad()`.
- `DOF_Group::setNodeDisp/Vel/Accel` **skip entries with `myID(i) < 0`** — a value set
  directly on the node for an `eqn<0` DOF survives the integrator's `setResponse`
  (`DOF_Group.cpp:565-642`).
- `AnalysisModel::updateDomain(time,dt)` order: `myDomain->applyLoad(time)` (→ each
  `SP::applyConstraint`, incl. `ImposedMotionSP` setting node vel/accel) **then**
  `myHandler->applyLoad()` (`AnalysisModel.cpp:604-606`). So a handler `applyLoad()` that
  sets the prescribed *disp* runs AFTER `ImposedMotionSP` set vel/accel — no clobber.
- CDL applies loads via `theModel->updateDomain(time,dt)` inside `newStep()` (its
  `setResponse(Ut,Vhalf,Azero)` precedes `updateDomain`, so free DOFs are scattered first,
  then prescribed DOFs imposed on their nodes).

## 1. Two tiers

The behavior splits on whether a prescribed DOF participates in an MP/EQ group.

### Tier 1 — prescribed DOF NOT in any constraint group (the common case)
Multi-support / uniform seismic excitation: support nodes carry imposed ground motion and
connect to the structure through **elements** (not MP ties). The free DOFs feel the support
through `F_int(u)` and the damping `C·v` — exactly how `Transformation`+`CentralDifference`
already works. **No projection interaction.**

Mechanism (handler-only; integrator unchanged):
1. Keep prescribed SP DOFs OUT of the equation set (`eqn = -1`, as today) — no mass needed.
2. Stop warning/dropping non-homogeneous SP. Record each prescribed `(node,dof)` + its
   `SP_Constraint*`.
3. Implement `LadrunoProjectionHandler::applyLoad()` to set the prescribed **displacement**
   on the node each step: `node->setTrialDisp(sp->getValue()+sp->getInitialValue(), dof)`
   (mirror `TransformationDOF_Group::enforceSPs`). `ImposedMotionSP` supplies vel/accel
   itself; plain SP leaves vel/accel = 0 (same limitation `Transformation` carries —
   documented, not a regression).
4. The free DOFs solve and integrate normally; their residual sees the imposed support
   state because `setResponse` skips the `eqn<0` support DOFs.

This alone makes base excitation correct under the projection handler — the dominant use.

### Tier 2 — prescribed DOF IS a master in an MP/EQ group (the literal "overwrite a")
e.g. a ground motion imposed on a `rigidDiaphragm` master, or a non-homogeneous SP on a
`rigidLink` master. The slaves must follow `a_c = C·a_master` with the master's PRESCRIBED
acceleration — not the current behavior, which drops a SP-fixed master and forces the
slave to 0 (`doneNumberingDOF` → `fixedEqns`). This is the case that needs the projector to
carry a known, non-zero master acceleration.

## 2. Tier-2 math (the known-RHS projection)

Partition a group's DOFs into: free-retained `f` (adjustable), **prescribed** `p` (accel
fixed each step to `a_p`), and constrained slaves `c` with `a_c = C_f a_f + C_p a_p`
(constant `C`, so the accel relation has no offset). The M-orthogonal projection minimizes
`Σ_i m_i (a_i − a_raw,i)²` over admissible `a`. The prescribed rows' term is constant
(`a_p` fixed) → drops out of the minimization over `a_f`. Moving the known `C_p a_p` to the
RHS of the slave rows:

```
ã_f = a_raw,f                       (free-retained rows)
ã_c = a_raw,c − C_p a_p             (slave rows; subtract the prescribed contribution)
a_f = (L_f^T M̃ L_f)^{-1} L_f^T M̃ ã     L_f = [I_f ; C_f],  M̃ = diag(m over {f,c} rows only)
a_c = C_f a_f + C_p a_p
a_p = prescribed (left unchanged)
```

The existing `project()` is exactly the `C_p` = ∅ special case. Generalization is
mechanical: `buildMass()` factors `L_f^T M̃ L_f` over **free** retained columns only
(prescribed columns excluded), and `project()` subtracts `C_p a_p` from slave targets and
leaves prescribed rows untouched. **Momentum:** generalized momentum is conserved on the
free directions; the prescribed DOF carries the reaction `m_p(a_raw,p − a_p)` (its tie
force) — exactly the support reaction, which is physically correct (momentum is NOT
conserved across an external prescribed-motion boundary).

### How the projector learns `a_p` each step
The prescribed master stays a numbered DOF? No — to avoid the mass requirement and match
Tier 1, prescribed masters ALSO stay `eqn = -1` and their `a_p` is imposed on the node
(`ImposedMotionSP` accel, or 0 for a constant SP). The handler resolves and stores `Node*`
+ `dof` for each prescribed row in the projector `Group`; `project()` reads
`node->getTrialAccel()(dof)` at call time. At `update()`-time `project()` the node accel is
the prescribed value (set during `newStep`'s `updateDomain`→`applyLoad`; `update`'s
`setResponse` hasn't run yet and would skip `eqn<0` anyway). This keeps the integrator hook
unchanged (still just `project(a)`) and localizes all prescribed logic in handler+projector.

> **Alternative considered & rejected:** keep prescribed masters IN the equation set and
> have the integrator overwrite `a[eqn]=a_p` before `project()`. Rejected: forces every
> prescribed DOF to carry nonzero mass (the wasted M⁻¹ solve), needs an integrator-side
> list of prescribed eqns + a way to fetch `a_p`, and diverges from Tier 1's out-of-eqn
> treatment. The node-read keeps ONE treatment for all prescribed DOFs.

## 3. Scope decision for THIS PR

- **Tier 1: IN.** The high-value, lower-risk core. Makes `imposedMotion` + non-homog SP
  correct for the dominant multi-support/base-excitation use.
- **Tier 2: ATTEMPT, gate-gated.** Implement the known-RHS generalization if the gate
  confirms the math + the node-read timing. If risk/complexity balloons, ship Tier 1 and
  defer Tier 2 to P4c — but in a Tier-1-only PR, a prescribed DOF used as a constraint
  MASTER must be **refused with a named error** (today it silently forces the slave to 0).

## 4. Seams

- `LadrunoProjectionHandler.{h,cpp}`:
  - `handle()`: replace the non-homog warn/drop with: record prescribed `(node,dof,SP*)`;
    keep `eqn=-1`. Refuse a prescribed DOF that is also an MP/EQ **slave** (overconstraint,
    like SP-on-slave). In a Tier-1-only build, refuse a prescribed **master**.
  - new `applyLoad()` override: set prescribed disp on each prescribed node.
  - `buildGroups()`/`doneNumberingDOF()` (Tier 2): tag prescribed master rows, pass
    `Node*`+dof + `C_p` columns to the projector instead of routing them to `fixedEqns`.
- `LadrunoConstraintProjector.{h,cpp}` (Tier 2): per-group prescribed rows (Node*,dof,
  free-vs-prescribed column split); `buildMass` factors over free columns; `project`
  subtracts `C_p a_p` (read from nodes) and leaves prescribed rows untouched; tie force
  on prescribed rows = `m_p(a_raw,p − a_p)` (reaction) — but `m_p` is on an `eqn<0` node
  with no diagonal-SOE mass, so the reaction tie-force is reported via the node mass or
  left 0 with a note (decide at gate).
- IC: `checkIC`/`snapICs` must treat a prescribed slave-row offset; prescribed DOFs start
  at `getInitialValue()`.
- CDL integrator: **no change for Tier 1.** Tier 2 also no hook change (node-read design).

## 5. Tests (planned `tests/test_adr30_projection_p4b.py`)
- **TP1** uniform base excitation (1-DOF-per-node frame, ground accel via `imposedMotion`)
  under LadrunoProjection+Diagonal vs the SAME model under Transformation+BandGen → match.
- **TP2** constant non-homog SP (prescribed support settlement) → static-offset response
  matches Transformation reference.
- **TP3** prescribed DOF used as an MP slave → named refusal.
- **TP4 (Tier 2)** ground motion on a `rigidLink` master → slave tracks `C·u_master`;
  vs Transformation reference.
- **TP5** the recorded support accel == the imposed ground-motion accel (round-trip).
- Regression: full Zone-A green.

## 6. Open questions for the gate
- **Q1** Tier-2 node-read of `a_p` at `project()` time — is the node accel guaranteed to be
  the prescribed value at every `project()` call (starter `a0` projection in `newStep`'s
  first-step block runs BEFORE that step's `updateDomain`/`applyLoad` — so at the starter
  the node may not yet hold `a_p`. Resolution candidate: starter reads `a_p` from the SP/GM
  at `t0` directly, or skips prescribed contribution at the starter).
- **Q2** Does `Domain::applyLoad`/`ImposedMotionSP::applyConstraint` get called for the
  starter and at the right time relative to the starter `project(a0)`?
- **Q3** Plain non-homog SP gives 0 vel/accel — confirm the free-DOF response is correct
  (matches Transformation, which has the same limitation) and document it.
- **Q4** Tier-2 reaction/tie-force reporting on a massless (`eqn<0`) prescribed node.
- **Q5** Interaction with `-projectICs` and the IC manifold check when prescribed DOFs move.
