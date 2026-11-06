---
wp: LEGACY
title: "Explicit zero-mass rotational DOF on an ndf=6 node ⇒ dt_cr silently overestimated"
legacy_seq: 73
---
### Explicit zero-mass rotational DOF on an ndf=6 node ⇒ `dt_cr` silently overestimated
- **Why it bites:** lumped element mass (beam lumped, `ASDShellQ4` rotational mass is
  EXPLICITLY omitted, `ASDShellQ4.cpp:1152`) and translational-only nodal `-mass` leave ZERO
  mass on the rotational dofs of an ndf=6 node ⇒ singular `M`. `CriticalTimeStep` does NOT
  warn — its DGGEV path FILTERS near-massless eigenpairs via a relative beta threshold
  (`betaTol = 1e-12*max|beta|`, `SRC/analysis/integrator/CriticalTimeStep.cpp:165-170`), so
  the zero-mass mode is dropped from omega_max and `dt_cr = 2/omega_max` comes back
  UNCONSERVATIVELY LARGE ⇒ the explicit run can go unstable with no diagnostic. In a mixed-ndf
  explicit model the binding constraint is MASS not stiffness: give every ACTIVE dof (incl.
  rotations on ndf=6 nodes) nonzero mass, use consistent mass, or restrain the massless dofs.
  Relevant to [[central_difference_ladruno_guide]]. See [[ndf_and_mixed_models_guide]] §7.
  2026-06-07.
