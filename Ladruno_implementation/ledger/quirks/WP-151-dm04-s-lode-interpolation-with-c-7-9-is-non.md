---
wp: WP-151
title: "DM04's Lode interpolation with c < 7/9 is NON-CONVEX on the extension meridian: the axisymmetric extension path is unstable, and a 1e-9 perturbation decides wh…"
legacy_seq: 532
---
### DM04's Lode interpolation with c < 7/9 is NON-CONVEX on the extension meridian: the axisymmetric extension path is unstable, and a 1e-9 perturbation decides whether a CTXu test liquefies (WP-151)
- **Bites:** g(θ) = 2c/((1+c) − (1−c)cos3θ) is a convex polar curve only for c ≥ 7/9. At θ = 60°, g = c, g′ = 0 and g″ = 4.5c(1−c), so r² − r·r″ ≥ 0 needs c ≥ 7/9 ≈ 0.78.
  - Both the TIMs campaign set (c = 0.71) and DM04's own Toyoura set (c = 0.712) are below it.
  - In undrained cyclic triaxial (campaign e0 0.6944, CSR 0.2), every model, DM04 included, breaks axisymmetry in the first extension half-cycle. |σ_yy − σ_zz| grows from round-off to 27–40 kPa.
  - A σ_zz perturbation of ±1e-9 decides between 5 % DA at N = 8 and no liquefaction by N = 20. Without a perturbation, round-off decides (rtol, build, any model change).
  - At c = 0.80 the path stays axisymmetric (|σ_yy − σ_zz| ≤ 1e-7 kPa) and the test is well conditioned.
  - **It is NOT the cause of the footing wall.** At c = 0.71 every `loadingNonPosH` refuser has n on the extension side, and those five states are regular at c = 0.80 (102/320 → 0/320, `sanisand_reseat_r1/fan_c080.py`).
  - But a c = 0.80 footing (`C080_EB_off`) walls anyway at s/B 0.048, on compression-side states (cos3θ ≈ +0.65). Exact DM04 runs the same Zeno re-seat sequence there: 119/576 trials fail, R1 0/576 (`fan_c080_bvp.py`).
  - The non-convexity only selects where the singular set is met first. Do not recalibrate c to fix the wall.
- **Rule:** With c < 7/9:
  - Do not read a single axisymmetric-extension element test (CTXu, TE) as the model's answer. Report both branches, or perturb explicitly.
  - A comparison of two model variants on such a test is decided by round-off unless the SAME perturbation is imposed on both.
  - A test harness that perturbs a state by wrapping a shared function must wrap the PRISTINE function. Re-wrapping per task compounds the perturbations inside a pool worker; WP-151's first gate run had this defect.
- **Workaround/status:** Not a code change: it is a calibration item for TIMs (keep c ≥ 0.78, or accept the ill-conditioning), raised in #887. [[151_sanisand_reseat_singularity]] §6.3.
