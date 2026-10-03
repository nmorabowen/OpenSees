---
wp: PR-525
title: "525 -- 7 vanilla row(s)"
pr: "#525"
files: ["`SRC/material/nD/NDMaterial.{h,cpp}`", "`SRC/material/nD/ElasticIsotropicPlaneStrain2D.{h,cpp}`", "`SRC/material/nD/J2PlaneStrain.{h,cpp}`", "`SRC/material/nD/PlaneStrainMaterial.{h,cpp}`", "`SRC/material/nD/UWmaterials/DruckerPragerPlaneStrain.{h,cpp}`", "`SRC/element/fourNodeQuad/FourNodeQuad.cpp`", "`SRC/element/triangle/Tri31.cpp`"]
table: "main"
legacy_seq: [269, 270, 271, 272, 273, 274, 275]
---
| `SRC/material/nD/NDMaterial.{h,cpp}` | `// Ladruno` plane-strain σ_zz exposure (apeGmsh request): add `virtual double getStressZZ(void)` defaulting to **quiet NaN** ("not available" — a material that doesn't implement it never reports a misleading 0). Strictly additive base virtual; `+#include <limits>`. | [#525](https://github.com/nmorabowen/OpenSees/pull/525) |
| `SRC/material/nD/ElasticIsotropicPlaneStrain2D.{h,cpp}` | `// Ladruno` σ_zz: override `getStressZZ()` = λ(ε₀+ε₁), λ = νE/((1+ν)(1−2ν)) — the exact elastic out-of-plane normal stress under ε_zz = 0. | [#525](https://github.com/nmorabowen/OpenSees/pull/525) |
| `SRC/material/nD/J2PlaneStrain.{h,cpp}` | `// Ladruno` σ_zz: override `getStressZZ()` = internal `stress(2,2)` — the 3D return mapping computes the full tensor; `getStress()` drops σ_zz (≠ ν(σxx+σyy) once plastic). | [#525](https://github.com/nmorabowen/OpenSees/pull/525) |
| `SRC/material/nD/PlaneStrainMaterial.{h,cpp}` | `// Ladruno` σ_zz: override `getStressZZ()` = wrapped 3D material's `getStress()(2)` — one override covers ANY 3D constitutive model adapted to plane strain (the highest-leverage site). | [#525](https://github.com/nmorabowen/OpenSees/pull/525) |
| `SRC/material/nD/UWmaterials/DruckerPragerPlaneStrain.{h,cpp}` | `// Ladruno` σ_zz: override `getStressZZ()` = `mSigma(2)` (base DruckerPrager keeps the full 6-vector; the plane-strain view drops σ_zz). | [#525](https://github.com/nmorabowen/OpenSees/pull/525) |
| `SRC/element/fourNodeQuad/FourNodeQuad.cpp` | `// Ladruno` σ_zz: new element response `"stressesPlaneStrain"`/`"stressPlaneStrain"` (responseID 21, Vector(16)) = `[σxx, σyy, σxy, σzz]` per GP, 4th slot from `getStressZZ()`, ResponseType sigma33. Existing `"stresses"` (ID 3) untouched. | [#525](https://github.com/nmorabowen/OpenSees/pull/525) |
| `SRC/element/triangle/Tri31.cpp` | `// Ladruno` σ_zz: same new `"stressesPlaneStrain"` response (ID 21, Vector(4*numgp)). Existing `"stresses"` untouched. | [#525](https://github.com/nmorabowen/OpenSees/pull/525) |
