---
wp: WP-196
title: "ASDShellQ4 with AGQI (EAS): the kind of the LAST evaluation before the next step is state"
date: 2026-10-14
---
### ASDShellQ4 with AGQI (EAS, the default): the last evaluation's kind persists into the next step

- **What bites.** Every `calculateAll()` re-forms the AGQI condensation state
  `m_eas->{KQU, KUQ, KQQ_inv, Q_residual}`; an LHS-only call (`getTangentStiff`, and
  `getDamp` / `getRayleighDampingForces` through it) leaves `Q_residual` ZERO, an RHS call
  fills it, and the next `update()` (`AGQIupdate`) reads it. So any extra post-commit
  evaluation (a recorder, an energy sweep) that ENDS on a different kind than before changes
  the trajectory — not just its own output.
- **Where it showed.** WP-196 first let EAS shells skip the energy kernel's WP-184 check
  (which ends on a tangent). The sweep then ended on a residual, and KE_nod and the nodal
  displacements moved in the last bits from step 1 (`test_wp196_energy_commit_phase.py`,
  EAS arm). With the legacy sweep the EnergyBalance recorder is neutral (same displacements
  with and without it).
- **Rule.** Treat "evaluate it fewer (or more) times" as a state change for EAS ASDShellQ4.
  `ASDShellQ4::ladrunoEnergyRayleighOnly()` answers false with EAS; `-noeas` shells (the
  Paper 1 decks) are evaluation-count neutral.
