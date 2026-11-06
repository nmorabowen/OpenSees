---
wp: PR-650
title: "650 -- 2 vanilla row(s)"
pr: "#650"
files: ["`SRC/analysis/integrator/TransientIntegrator.cpp`", "`SRC/analysis/integrator/IncrementalIntegrator.cpp`"]
table: "main"
legacy_seq: [120, 121]
---
| `SRC/analysis/integrator/TransientIntegrator.cpp` | `// Ladruno ADR-77 T0`: **supersedes the "DOF_Group loop is left untimed" decision in the #49 row above.** Adds a coarse `OPS_PROFILE_SCOPE("dof.tangent")` around the DOF_Group `addA` loop in `formTangent`, re-braced to scope it. This loop carries the **nodal** mass/damping into A and exists ONLY on the transient path (`IncrementalIntegrator::formTangent` has no DOF_Group loop at all) — leaving it untimed silently folded nodal-mass assembly into "the rest of `formTangent`", which is precisely the five-loop gap ADR-75 P1f flagged. Coarse, not element-keyed (DOF_Groups have no classTag). Profiling-only; inert unless the profiler is enabled. | [#650](https://github.com/nmorabowen/OpenSees/pull/650) |
| `SRC/analysis/integrator/IncrementalIntegrator.cpp` | `// Ladruno ADR-77 T0`: adds a coarse `OPS_PROFILE_SCOPE("dof.residual")` in `formNodalUnbalance` (the nodal inertia/damping/load side of B), mirroring `dof.tangent`. Shared by the static and transient paths. Profiling-only; inert unless enabled. Measured cost on the ADR-77 Lane-B deck: 0.06-0.08% of step (the deck carries mass in the material, so this loop is nearly empty there — a lumped-nodal-mass deck would show more). | [#650](https://github.com/nmorabowen/OpenSees/pull/650) |
