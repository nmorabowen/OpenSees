---
wp: LEGACY
title: "RC fixed-crack cyclic PINCHING is a panel phenomenon, not a material-point one"
legacy_seq: 77
---
### RC fixed-crack cyclic PINCHING is a panel phenomenon, not a material-point one
- **Bites:** expecting a pinched (waisted) `τ_nt`–`γ_nt` hysteresis loop from a single material point / one element under homogeneous cyclic shear at constant normal strain. You get a FAT loop instead.
- **Why:** with the crack frozen and `en` (hence `v_ci,max`) constant, the only nonlinearity is the `±v_ci,max` clamp; the elastic interlock band has width `2·v_ci,max/G ≈ 5e-4` in slip — sub-step at normal resolution — so the loop is essentially a `±v_ci,max` rectangle (maximum dissipation, zero pinch). A pinched waist requires the crack to OPEN/CLOSE during the cycle (`v_ci,max` low near slip-reversal, high at slip-peaks), which needs the principal direction to ROTATE relative to the fixed crack — i.e. a real panel / non-homogeneous stress field. So Phase-2b.1 material-point tests assert the MECHANISM (reversal re-cap, unload=G, closure cap-recovery, energy>0), and the pinching-shape + hysteretic-energy acceptance is a panel/experiment (Tran–Wallace) gate deferred to 2b.2. Learned 2026-06-16.
