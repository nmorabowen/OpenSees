---
wp: LEGACY
title: "Eigenvalue-floored 2x2 condensation: scale the floor by the DIAGONAL only, not the off-diagonal"
legacy_seq: 90
---
### Eigenvalue-floored 2x2 condensation: scale the floor by the DIAGONAL only, not the off-diagonal
- **Bites:** the biaxial hinge condenses a coupled 2x2 `K_aa` via an eigenvalue-floored inverse (`ladrunoFlooredInv2x2`), with the floor magnitude `1e-8*(|bulk_zz|+|bulk_yy|+|Czz|+|Cyy|)`. Folding the off-diagonal coupling `|Kaa_zy|` into that sum (to "scale with the full matrix") **over-floors** the modes near the elliptical onset — the floored inverse is then too damped, the radial inner Newton oscillates and hits `maxIter`, and 4 coupled tests break. The diagonal-only scale is the validated choice; the off-diagonal already enters the eigenvalues themselves, so it must NOT also inflate the floor. Learned 2026-06-18 (same review).
