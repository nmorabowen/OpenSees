---
wp: LEGACY
title: "Tension-stiffening pinning tangent must span ALL 6 columns of the floored rows, not just in-plane"
legacy_seq: 83
---
### Tension-stiffening pinning tangent must span ALL 6 columns of the floored rows, not just in-plane
- **Bites:** the tension-stiffening floor pins `n^Tσn` to `σ_ts(ε1)` independent of the bare stress. The bare in-plane normal stress depends on **out-of-plane** strain too (eps_zz via the elastic λ coupling), so pinning it removes that dependence. A consistent tangent that subtracts the bare sensitivity `d(n^Tσn)/dε` only over the in-plane columns `{0,1,3}` leaves the **eps_zz column** (and any out-of-plane shear columns) at the now-pinned-away bare value → a forward-difference tangent check shows ~0.98 rel-error at `D[*][2]`.
- **Fix (proven):** compute `row[c] = ts_meas·Dtan[in-plane rows][c]` for ALL c in 0..5 and subtract over all 6 columns; the `dσ_ts/dε1·de1/dε` add-back is membrane-only (nonzero only on `{0,1,3}`). Note the BASELINE W_B secant tangent is itself only approximate in the softening regime (it omits `d(dt_bar)/dε`), so a global FD-tangent gate "fails" for both on and off; the right TS tangent checks are (a) the PINNED direction `D[0][0]==dσ_ts/dε1` (exact, since `σ0=σ_ts(ε1)` is smooth) and (b) TS does not DEGRADE the worst FD rel-error vs the off baseline. Learned 2026-06-18, [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3a (`tests/_testbed/rc_tensstiff_gpp.cpp`).

#### Follow-up (WP-241)
The `-interlock` tangent had the same defect and kept it until WP-241: the 2a-capped and the `-cyclic` blocks removed the smeared crack-shear sensitivity `m_σ·Dtan` over columns `{0,1,3}` only, so the ε_zz column kept the pinned-away value (≈3.5e3 for E = 30 GPa) in every condensed view (PlaneStress, PlateFiber). Fixed the same way: `row[j]` over all 6 columns; the friction stiffness `G m_ε ⊗ m_ε` stays in-plane. Gate: `tests/_testbed/rc_spatial_gpp.cpp` G5 (with `betaSrMin = 0`, `m_σ·D = 0` in all columns and the FD of the pinned shear w.r.t. ε_zz ≈ 0). The general rule: **any pinning of a stress projection must remove its sensitivity to every strain component, not only to the components of its own plane.**
