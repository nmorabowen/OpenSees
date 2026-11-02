---
wp: WP-220
title: "Path-control integrators built phat with the prior residual inside it"
date: 2026-11-02
---
### Load factor lags the element force after re-issuing `integrator DisplacementControl`, or adding a pattern at t != 0; residual stays at the drive DOF under `NormDispIncr`

- **Bites:** `ArcLength`, `ArcLength1`, `DisplacementControl`, `DistributedDisplacementControl`, `EQPath`,
  `HSConstraint`, `MinUnbalDispNorm`, `LadrunoArcLength`, `LadrunoIndirectControl`, when `domainChanged`
  runs with a non-zero unbalance. `StaticAnalysis::setIntegrator` resets `domainStamp`, so every re-issued
  integrator rebuilds phat.
- **Why:** phat = R(lambda+1) = q + R_prev. At a drive DOF r_new = r R_prev / (q + R_prev), and dU = 0
  makes `NormDispIncr` pass with r still present. With q = 1 N it fed on itself: 76 kN at the interior
  node of two trusses in series.
- **Fix (WP-220):** phat = R(lambda+1) - R(lambda), `LadrunoReferenceLoad.h`. A pattern on a `Constant` series
  has no such difference and keeps R(lambda+1), so the defect remains there. Identical when the unbalance is zero. Runs that
  re-issue the integrator get a different load factor, so anything read from it needs re-baselining.
- **`enforceSPs`:** the handler now returns -1 when an element fails to update. `updateDomain(t, dT)` honours it,
  so a transient step fails at newStep. `applyLoadDomain` still ignores it: a static step fails later in
  `Domain::update` with the same -3 as before, and tests pin that code.
- **Gate:** `tests/test_wp220_reference_load_and_sp_failure.py`. `EQPath` has no Python command and
  `DistributedDisplacementControl` is not exercised there.
