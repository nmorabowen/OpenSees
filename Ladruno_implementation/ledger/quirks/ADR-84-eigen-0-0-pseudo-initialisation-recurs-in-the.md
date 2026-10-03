---
wp: ADR-84
title: "Eigen *= 0.0 pseudo-initialisation recurs in the YF/PF headers — the ADR-84 constructor trap, four more times (ADR-94 B4)"
legacy_seq: 390
---
### Eigen `*= 0.0` pseudo-initialisation recurs in the YF/PF headers — the ADR-84 constructor trap, four more times (ADR-94 B4)
- **Bites:** `VoigtVector pressure_part; pressure_part *= 0.0;` in `DruckerPrager_YF.h:64-65` and `DruckerPrager_PF.h:68-69`; an uninitialised `VoigtVector zero;` returned by `NullHardeningTensorPolicy`; AF's saturation-branch `derivative`. `EIGEN_INITIALIZE_MATRICES_BY_ZERO` is defined nowhere in the tree, so stale non-finite heap bits survive and DP is the only YF that NaNs at its apex.
- **Rule:** `setZero()` or `VoigtVector::Zero()`, never `*= 0`. Grep the idiom before trusting any new Eigen-backed component.
- **Platform-dependent:** Ubuntu CI (fresh heap) commits a clean finite history on the same path; only the dirty Windows pytest heap reproduces the NaN — which is the signature of UB, not of a deterministic defect.
- **Status (wp/94a):** FIXED at all six sites found by the grep — `DruckerPrager_YF.h` / `DruckerPrager_PF.h` (`pressure_part` and the degenerate `dev_part` branch), `NullHardeningTensorPolicy`, `ArmstrongFrederickPolicy`'s saturation `derivative`, `DuncanChang_EL.h`'s `EE_MATRIX`, plus 27 `*= 0` statements on ASDPlasticMaterial3D's own class-static `depsilon`/`dsigma`/`intersection_*`. The **rule is the durable part**: grep `\*= 0` before trusting any new Eigen-backed component; the idiom has now been found three times in this one subsystem (ADR-84 P0, ADR-94 B4, ADR-94 wp/94a's own sweep).
