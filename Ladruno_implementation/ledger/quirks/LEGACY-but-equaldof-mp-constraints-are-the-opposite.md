---
wp: LEGACY
title: "…but equalDOF/MP constraints are the OPPOSITE — added mid-stage they PRESERVE the offset (no snap)"
legacy_seq: 35
---
### …but `equalDOF`/MP constraints are the OPPOSITE — added mid-stage they PRESERVE the offset (no snap)
- **The asymmetry (this is the surprising part):** unlike SP, an MP constraint
  (`equalDOF`, `rigidLink`, `rigidDiaphragm`) added after a node has displaced does
  **not** snap the constrained node onto the retained one. It ties their *future
  increments* together while **preserving the relative offset that existed at tie-time**,
  with **zero initial constraint force**. This is the "install at the current deformed
  state" behavior you'd *wish* `fix` had — and for MP it's the default, no flag needed.
- **Why:** `MP_Constraint::setDomain()` captures BOTH nodes' current disps at add-time
  — `Uc0` (constrained), `Ur0` (retained) (MP_Constraint.cpp:294-313) — and every
  MP-capable handler enforces the relation on the **offset-removed** displacements,
  *unconditionally* (no `retZeroInitValue` equivalent): Penalty/Lagrange build the
  residual from `Uc - Uc0` and `Ur - Ur0` (PenaltyMP_FE.cpp:230-238 → equilibrium is
  `(Uc - Uc0) = C·(Ur - Ur0)`, not `Uc = C·Ur`); the Transformation handler under
  `TRANSF_INCREMENTAL_MP` transforms only the **increment** (`modUnbalance -=
  modTrialDispOld`, TransformationDOF_Group.cpp:525) and applies it via `incrTrialDisp`
  (line 560), so the standing offset is never overwritten. At tie-time `Uc=Uc0`,
  `Ur=Ur0` ⇒ constraint satisfied with zero force and zero jump.
- **Net rule:** SP (`fix`/`sp`) defaults to enforcing the **absolute** value → snaps to
  reference; MP (`equalDOF`/rigid) defaults to enforcing the **increment** → preserves
  the current offset. Same "capture init disp at add-time" machinery underneath,
  **opposite defaults** (MP was hardened for staged construction; SP kept its legacy
  absolute-value default and never exposed the incremental toggle to a command flag).
  Caveat: holds for the MP-capable handlers (Transformation / Penalty / Lagrange); the
  Plain handler isn't the one to use for nontrivial MP staged ties. Learned 2026-05-31.
  Full write-up with source trail: [[constraints_reference_position]].
