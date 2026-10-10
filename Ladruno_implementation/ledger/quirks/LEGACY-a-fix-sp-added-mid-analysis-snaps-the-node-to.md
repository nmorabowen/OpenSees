---
wp: LEGACY
title: "A fix/sp added mid-analysis snaps the node to the REFERENCE frame, not the current deformed shape"
legacy_seq: 34
---
### A `fix`/`sp` added mid-analysis snaps the node to the REFERENCE frame, not the current deformed shape
- **Bites:** in a staged analysis (deform under stage 1, hold load, then constrain
  part of the model), calling `fix nodeTag dof` on a node that has already displaced
  by `d` drives that DOF back toward `u = 0` on the next `analyze`, dragging the node
  to its **original undeformed location** and dumping spurious strains/forces into the
  attached elements. People expect the new support to "catch" the structure at its
  current deformed shape; it does the opposite.
- **Why — the conceptual trap:** an SP constraint prescribes the **total value of the
  displacement DOF**, `u = value` (`fix` ⇒ `u = 0`), and `u` is *always* measured from
  the original mesh at t=0. There is **only one displacement frame and it never
  re-zeros** — not at a stage boundary, not ever. Since `position = X_ref + u`, the
  statements "fix the deformation to zero" and "send the node back to its original
  position" are *identical* (`u=0 ⟺ position = X_ref`). The constraint is an algebraic
  equation on absolute `u`, not an incremental/ratchet condition on the change-from-now,
  so adding it later does **not** rebase `u` to the current state. "Constraints fix
  deformation, not position" is true but misleading — it's deformation *measured from
  the undeformed configuration*.
- **Source mechanics (this build):** the current displacement at constraint-add time
  *is* captured — `SP_Constraint::setDomain()` records `initialValue = U(dofNumber)`
  (SP_Constraint.cpp:380) — and all three handlers are wired to subtract it (Penalty
  `resid = alpha*(constraint - (nodeDisp - initialValue))`, PenaltySP_FE.cpp:139;
  Lagrange LagrangeSP_FE.cpp:143; Transformation under `#define TRANSF_INCREMENTAL_SP`
  in TransformationDOF_Group.h:44 → TransformationDOF_Group.cpp:1055). BUT that
  "stay-in-place" path only fires when `retZeroInitValue == false`, and
  `getInitialValue()` returns `0` whenever it's `true` (SP_Constraint.cpp:317).
  **`fix` and `sp` both default `retZeroInitValue = true`**, and the `sp -subtractInit`
  flag *also* just sets it `true` (OpenSeesPatternCommands.cpp:1065) — so the
  incremental/stay-in-place branch is compiled in but **not cleanly reachable from
  script**. Default behavior across all handlers = snap to reference.
- **Workaround:** to install a support that holds the *current* deformed position with
  zero initial force, prescribe the current displacement explicitly rather than `fix`:
  `d = ops.nodeDisp(n, dof); ops.sp(n, dof, d, '-const')` (needs an active pattern or
  `-pattern N`). To genuinely return the DOF to its t=0 position, `fix` is correct and
  the forces are physical. In dynamics, any sudden BC change also injects an impulse;
  ramp the prescribed value via a timeSeries. Learned 2026-05-31.
