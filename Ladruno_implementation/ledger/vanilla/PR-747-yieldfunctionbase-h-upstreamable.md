---
wp: PR-747
title: "747 -- upstreamable-table row(s)"
pr: "#747"
files: ["`SRC/material/nD/ASDPlasticMaterial3D/YieldFunctionBase.h`", "`SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h`"]
table: "upstreamable"
legacy_seq: [468, 469]
---
| `SRC/material/nD/ASDPlasticMaterial3D/YieldFunctionBase.h` | `// Ladruno (ADR-84 P3)`: the `SPECIAL_RETURN` macro gains a trailing `int& return_quality` out-param, and two `SR_QUALITY_*` codes are defined next to it. `stiffness_return` is now contractually the **RAW active-set (Koiter) tangent** — the YF must NOT apply a tangent-operator policy. Only `MohrCoulombTensionCutoff_YF` implements the macro, so the signature change has exactly one implementor and one call site. | [#747](https://github.com/nmorabowen/OpenSees/pull/747) |
| `SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h` | `// Ladruno (ADR-84 P3)`: at the `special_return` call site in `Backward_Euler`, (a) apply the material's configured `INT_OPT_tangent_operator_type` to the raw tangent the hook returns — `Elastic → Eelastic`, `Continuum`/`Algorithmic` → raw, `Secant` (**the default**) → `(D+E)/2`, `Numerical_Algorithmic_*` → raw (they differentiate `compute_local_stress()`, which never calls the hook, so they would be strictly worse); and (b) when the hook reports `SR_QUALITY_FALLBACK` (the terminal vertex projection, which discards the deviatoric state) and `strict_convergence` is on, print an `opserr` line and return -1 instead of committing it. **P0 assigned `Stiffness = stiff_sr` unconditionally and the YF had already secant-blended it, so `tangent_type` was silently overridden on every hook-resolved Gauss point** — the measured cause of the Cerro Lindo M5 confined-shear stall (87% tangent error at the corner; ADR-84 §9). Default `Secant` keeps this byte-identical. | [#747](https://github.com/nmorabowen/OpenSees/pull/747) |
