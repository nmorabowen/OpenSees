---
wp: ADR-78
title: "u-p storage coupling under corot: velocity contraction is chord-poisoned; incremental-coupling-only is a pump (ADR 78)"
legacy_seq: 243
---
## u-p storage coupling under corot: velocity contraction is chord-poisoned; incremental-coupling-only is a pump (ADR 78)

- **Bites (two ways, both measured):** (1) contracting the damp p-row `(R̄Q)ᵀ`
  against integrator velocities of a ROTATING body picks up the chord's
  apparent volumetric rate `2(1−cosΔθ)/Δt` per step — one-signed (a systematic
  dilation in the current frame; never averages out),
  amplified by `Q̄ ≈ K_f/n` undrained (order-1 spurious p at bearing-mechanism
  rotations); no velocity-linear operator can remove it. (2) Fixing ONLY the
  coupling to the incremental `QᵀΔu_d/Δt` while leaving `S·ṗ` on the Newmark
  velocity breaks the skew-symmetry of the discrete coupling pair and PUMPS the
  structural ringing (consolidation column grew a ±100·q p-oscillation at
  Δt=0.02 where linear decays).
- **Workaround/status:** the WHOLE p-row rate block goes incremental under
  corot — `QᵀΔu_d/Δt + (S+αH̃)Δp/Δt` (the Book's GN22/GN11 pairing);
  first-order convergent (measured 2.3e-2→4.2e-3 over Δt 0.16→0.02) and
  rigid-motion exact. Full analysis: ADR 78 §3.3. *2026-07-28 (ADR 78).*
