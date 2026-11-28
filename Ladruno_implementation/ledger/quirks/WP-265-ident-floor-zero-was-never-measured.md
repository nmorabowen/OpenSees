---
wp: WP-265
title: "A 0.0 mutation-gate IDENT floor meant 'never measured', not 'the suite cannot see the tangent'"
pr: "#1179"
date: 2026-11-28
status: "open"
---
### Mutation gate IDENT score / floor is 0 for every family: does nothing test the tangent?

- **Bites:** reading `"IDENT": 0.0` in `Ladruno_scripts/mutation_gate.py` as
  "no test pins the tangent". Also reading a high IDENT kill rate as
  "the tangent is verified".
- **Why:** before WP-265 no IDENT mutant was ever scored for CONTINUUM,
  SHELLMOD or CPML, because `ladruno_mutation.yml` builds ZERO and SCALE only.
  The floors were declared from the header's claim that Newton "usually still
  converges", and that claim is wrong for stiff static models. With K = I,
  dU = R diverges whenever the stiffness eigenvalues leave (0, 2). A `Linear`
  step lands at u = F. Implicit dynamics survives only where a0*M dominates.
  So IDENT kills most static tests by divergence. That proves only that a
  tangent exists. The tests that prove the tangent is RIGHT compare K with a
  finite difference of the force at a nonlinear state (`tests/*tangent_gate*`).
  IDENT also moves residuals built from K: Rayleigh betaK/betaK0,
  `criticalTimeStep` (reads `getInitialStiff`), shell drilling stiffness from
  the section initial tangent (ShellMITC4, ASDShellQ4), and ASDShellQ4's EAS
  residual. Kills on those paths are residual changes.
- **Fix (WP-265):** each gated family's IDENT suite (`paths_by_mode`) is its
  suite plus a K-vs-FD gate file. IDENT floors are > 0 and provisional until
  measured. `tests/test_mutation_tangent_gate_registry.py` fails if a gated
  family's floor returns to 0, or if a gate file leaks into ZERO/SCALE (where K
  and the FD are mutated consistently, so it can only survive).
