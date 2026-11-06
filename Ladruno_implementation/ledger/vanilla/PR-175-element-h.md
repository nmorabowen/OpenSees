---
wp: PR-175
title: "175 -- 1 vanilla row(s)"
pr: "#175"
files: ["`SRC/element/Element.{h,cpp}`"]
table: "main"
legacy_seq: [135]
---
| `SRC/element/Element.{h,cpp}` | `// Ladruno` (ADR 20 §9): add `virtual int getInterpolationWeights(const Vector& xi, Vector& N)` to the `Element` base — default returns −1 (not implemented), host elements override. Single source of truth for nodal shape-function weights so `LadrunoEmbeddedRebar` can embed a rebar node via `-host eleTag -xi …` instead of re-supplied `-shape`. Additive base-class virtual (vtable change ⇒ recompile-all, but no existing behavior touched). Overridden by fork hosts `LadrunoBrick` (trilinear) + `BezierTet10` (Bernstein). | [#175](https://github.com/nmorabowen/OpenSees/pull/175) |
