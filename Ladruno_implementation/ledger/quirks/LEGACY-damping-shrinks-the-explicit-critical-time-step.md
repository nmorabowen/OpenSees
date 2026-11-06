---
wp: LEGACY
title: "Damping shrinks the explicit critical time step — and βK is a trap"
legacy_seq: 24
---
### Damping shrinks the explicit critical time step — and βK is a trap
- **Bites:** adding damping *reduces* `dt_cr` in both coupled and explicit modes;
  stiffness-proportional damping (`βK`) is especially punishing under explicit
  integration.
- **Why:** βK inflates the effective high-frequency content that sets `dt_cr`.
- **Status:** use mass-proportional (`αM`) damping in explicit; captured in the
  CentralDifferenceLadruno plan (`project_robust_central_difference`).
