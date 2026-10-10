---
wp: PR-222
title: "222 -- 1 vanilla row(s)"
pr: "#222"
files: ["`SRC/material/nD/InitStrainNDMaterial.{h,cpp}`"]
table: "main"
legacy_seq: [173]
---
| `SRC/material/nD/InitStrainNDMaterial.{h,cpp}` | `// Ladruno`: make Petracca's fixed-prestrain `InitStrain` wrapper **dimension-general**. Was 3D-only — `getOrder()` hardcoded 6 and `getCopy(type)` only answered `"ThreeDimensional"`, so the base `NDMaterial::getCopy()` returned `0` for `PlaneStrain`/`AxiSymmetric` (the **default** view for `LadrunoQuad`/`LadrunoCST`) → null material → element construction failed. Fix: add an adopting ctor that preserves an already type-reduced inner view; `getCopy(type)` builds `PlaneStrain`/`AxiSymmetric` natively (inner = that view of the stored 3D material; `eps0` reduced via the LadrunoJ2 `vmap` {0,1,3} / {0,1,2,3}); `getOrder()` delegates to the inner; `setTrialStrain(Incr)` sized to the inner order; `getCopy()` clones preserving the current view. Keeps a canonical size-6 `epsInit3D` so the `eps0_ij` setParameter API has stable indices in every view (updateParameter edits 3D + re-projects; out-of-view comps are no-ops, no OOB). `sendSelf/recvSelf` carry `epsInit3D`+active order. `PlaneStress`/`PlateFiber`/`BeamFiber`/3D unchanged (legacy base-wrapper path). Additive — 3D behavior byte-identical. | [#222](https://github.com/nmorabowen/OpenSees/pull/222) |
