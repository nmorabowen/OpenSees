---
wp: ADR-58
title: "Transformation handler is INCREMENTAL (TRANSF_INCREMENTAL_MP); an element's update() setTrialDisp on a constrained slave SURVIVES (ADR-58 P2-S2)"
legacy_seq: 137
---
## Transformation handler is INCREMENTAL (`TRANSF_INCREMENTAL_MP`); an element's `update()` setTrialDisp on a constrained slave SURVIVES (ADR-58 P2-S2)

`TransformationDOF_Group::setNodeDisp` is compiled with `TRANSF_INCREMENTAL_MP` (defined in the
header), so it imposes a constrained (MP-eliminated) node as `slave += T·δu_retained` —
**incremental**, off the retained node's per-step increment, homogeneous (no constant term; `Uc0/Ur0`
are never read on this path). Two consequences exploited for the rigid body's finite-rotation slaving:
1. A linear/time-varying MP can carry only what `T·u_retained` spans — it **cannot** place a slave at
   a nonlinear `(R−I)d⁰` offset whose source (the body orientation `q`) is not a retained DOF. A
   `setTrialDisp` done *inside* a custom MP's `getConstraint()` is **overwritten** by the subsequent
   `incrTrialDisp(T·δu)` (the MP_Joint3D "length-correction" is likewise a slow cross-step nudge only).
2. BUT a `setTrialDisp` done in the **owning Element's `update()`** SURVIVES into residual formation:
   per step the order is `newStep` → `applyLoad`/enforceSPs → `Domain::update()` (element `update()`
   imposes the absolute slave position) → `formUnbalance` (reads it). `setResponse` in the post-solve
   `update(U)` does the handler's incremental write, then `Domain::update()` re-imposes — the element
   always wins the **last write before the next residual**. Verified: slaves track `(R−I)d⁰` to 5e-16.
   This is why the rigid body imposes slaves directly (mechanism "C3") instead of via a custom MP.
   Caveat: the committed slave then lags `q` by one step (the element's `update()` uses the pre-commit
   `qTrial`; `q` advances later in `commitState`) — a standard explicit half-step offset, and the slave
   **velocity/accel are not re-imposed**, so velocity-dependent elements on a slave see an inconsistent
   `v` (keep incident elements displacement-only, or also impose `v_i = v_R + ω×(R·d⁰)`).
