---
wp: LEGACY
title: "Mortar tangential SLIP must come from DISPLACEMENTS, not positions — the closest-point projection makes the weighted relative POSITION purely normal"
legacy_seq: 124
---
### Mortar tangential SLIP must come from DISPLACEMENTS, not positions — the closest-point projection makes the weighted relative POSITION purely normal
- **Bites:** the ADR-41 C3.1 mortar friction. The natural-looking weighted relative position
  `r_I = Σ_J D_IJ x_s,J − Σ_K M_IK x_m,K` (= `∫N_I(x_s − x_m(ξ̄)) dΓ`) is **purely NORMAL**: `n·r_I = g̃_I`
  (the weighted normal gap) and its TANGENTIAL part is ≈ 0, because the closest-point projection `ξ̄`
  places the master point directly "under" the slave point (`x_s − x_m(ξ̄) ∥ n` by construction). So a
  friction slip built from positions is ~0 even when the slave has slid a finite tangential distance —
  the return map sees `gTeff ≈ 0`, stays in STICK, and assembles ZERO friction force (symptom: a driven
  block accelerates at the frictionless `a = Q/m`, friction silently inert, to 1e-13). Verified by a
  stderr probe: `gTeff=(2.8e-17, …)` at a step where the slave had displaced `x=1e-3`.
- **Why:** mortar inherits the NTS lesson — the ADR-39 `SEGMENT` path's `segmentActive` ALREADY documents
  "the closest-point projection makes (x_s − x̄) ∥ n, so POSITIONS carry NO tangential information; the slip
  is the slave DISPLACEMENT minus the interpolated master DISPLACEMENT at the projection: `d = u_s − Σ N_i u_i`."
  The C3.1 first draft re-made the position mistake the NTS path had already solved.
- **Fix (shipped C3.1):** build the slip from DISPLACEMENTS — `r_I = Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K`
  (`u = getTrialDisp()`), tangential part `/a_I`, minus the engagement origin `gT0_I`. This is the `D/M`-
  weighted generalisation of the NTS `u_s − Σ N_i u_i`. The normal gap still uses positions (it IS the
  normal projection); only the tangential slip switches to displacements. General lesson: in any
  closest-point-projected contact, the normal gap is a POSITION quantity and the tangential slip is a
  DISPLACEMENT quantity — they are not interchangeable. Found while bringing up C3.1 (the driven-block
  gate caught it). See [[_adr41_c3_design]] §mechanics step 1, #377.
