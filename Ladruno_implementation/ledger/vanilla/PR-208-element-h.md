---
wp: PR-208
title: "208 -- 1 vanilla row(s)"
pr: "#208"
files: ["`SRC/element/Element.{h,cpp}`"]
table: "main"
legacy_seq: [151]
---
| `SRC/element/Element.{h,cpp}` | `// Ladruno` (ADR 23 §3, Phase 2 UR): add `virtual int getInterpolationGradients(const Vector& xi, Matrix& dNdx)` to the `Element` base — default returns −1 (not implemented). The gradient companion of `getInterpolationWeights`: fills `dNdx(i,j)=∂N_i/∂x_j` so `LadrunoEmbeddedNode`'s rotation (UR) tie can read the host continuum rotation `θ=½curl(u)=skew(∇u)` at the embedded point (weights alone can't — the rotation needs `∂N/∂x`). Overridden by fork hosts `LadrunoBrick` (via `shp3d`) + `BezierTet10` (via `computeJacobian`). Additive base-class virtual (vtable change ⇒ recompile-all, no existing behavior touched). | [#208](https://github.com/nmorabowen/OpenSees/pull/208) |
