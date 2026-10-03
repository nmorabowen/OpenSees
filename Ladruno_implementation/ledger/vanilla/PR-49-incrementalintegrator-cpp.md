---
wp: PR-49
title: "49 -- 2 vanilla row(s)"
pr: "#49"
files: ["`SRC/analysis/integrator/IncrementalIntegrator.cpp`", "`SRC/analysis/integrator/TransientIntegrator.cpp`"]
table: "main"
legacy_seq: [118, 119]
---
| `SRC/analysis/integrator/IncrementalIntegrator.cpp` | `// Ladruno P3`: deep per-element-type timing — `OPS_PROFILE_SCOPE_DEEP_NAMED` around the `formTangent` and `formElementResidual` FE_Element loops + `OPS_PROFILE_FE_ELEM_SCOPE` per element (folds each `getTangent`/`getResidual` wall into a per-classTag `elem_by_type` bucket). Loops re-braced to scope the per-element timer; `getElement()`/`getClassTag()` touched ONLY when the deep gate is on (`tmr.engaged()` short-circuit) so an unprofiled run pays nothing. +`#include <Element.h>`, `<profiler/ProfilerMacros.h>`. Additive/behavior-preserving. | [#49](https://github.com/nmorabowen/OpenSees/pull/49) |
| `SRC/analysis/integrator/TransientIntegrator.cpp` | `// Ladruno P3`: same deep per-element-type timing on the FE_Element tangent loop of `TransientIntegrator::formTangent` (the transient/explicit assembly path; residual path reuses `IncrementalIntegrator::formElementResidual`). The DOF_Group mass/damping loop is left untimed (not element-keyed). +`#include <Element.h>`, `<profiler/ProfilerMacros.h>`. Additive. | [#49](https://github.com/nmorabowen/OpenSees/pull/49) |
