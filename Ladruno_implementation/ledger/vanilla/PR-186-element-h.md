---
wp: PR-186
title: "186 -- 2 vanilla row(s)"
pr: "#186"
files: ["`SRC/element/Element.{h,cpp}`", "`SRC/analysis/integrator/CriticalTimeStep.cpp`"]
table: "main"
legacy_seq: [150, 152]
---
| `SRC/element/Element.{h,cpp}` | `// Ladruno` (ADR 20 §10.6.1): add `virtual double getExplicitCriticalTimeStep()` to the `Element` base — default returns −1 ("no opinion"). Lets an element self-report an explicit critical step its per-element `K v = λ M v` pencil can't express; overridden by `LadrunoEmbeddedRebar` to surface its bipenalty bound `2√(m_p/k_eff)` (invisible to the eigensolve, which sees the massless-host coupling as `λ_max=0`). Additive base-class virtual (vtable change ⇒ recompile-all, no existing behavior touched). | [#186](https://github.com/nmorabowen/OpenSees/pull/186) |
| `SRC/analysis/integrator/CriticalTimeStep.cpp` | `// Ladruno` (ADR 20 §10.6.1): at the top of the per-element loop in `computeCriticalTimeStep`, query `ele->getExplicitCriticalTimeStep()`; a non-negative value is folded into the running damped/undamped minima and the per-element eigensolve is skipped for that element. Default −1 ⇒ every existing element takes the unchanged eigensolve path. Makes `ops.criticalTimeStep`/`-cflAbort` honor a self-reporting element's bound. | [#186](https://github.com/nmorabowen/OpenSees/pull/186) |
