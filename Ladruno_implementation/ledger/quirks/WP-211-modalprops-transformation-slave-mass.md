---
wp: WP-211
title: "modalProperties drops slave-node mass under constraints Transformation"
date: 2026-10-25
status: "fixed-by WP-211"
---
### `modalProperties` gives `totalFreeMass` = 0 and Γ = 0 (RSA returns zero) for a rigid diaphragm with mass on the slave nodes

- **Bites:** `rigidDiaphragm` / `equalDOF` / any MP constraint, `constraints
  Transformation`, nodal mass on the constrained (slave) nodes — e.g. a meshed
  floor slaved to a diaphragm master (apeGmsh's default handler for MP models).
  `eigen` periods and `nodeEigenvector` at the slaves are right, but
  `modalProperties -return` shows `totalMass` full, `totalFreeMass` = 0, every
  `partiFactorM*` / `partiMassM*` = 0, and `responseSpectrumAnalysis` gives zero
  displacements, reactions and base shear. Penalty, Lagrange, or the mass lumped
  at the master were all right.
- **Why:** `DomainModalProperties` treated every DOF with DOF_Group ID -1 as
  fixed; Transformation numbers MP-retained DOFs -1 too (they are recovered from
  the master), so their mass dropped out of M, the free mass, the centre of mass
  and the ratio denominators.
- **Fix/status:** WP-211 — fixed means ID -1 **and** an SP on that (node, dof).
  An SP on top of an MP-retained DOF under Transformation stays fixed (mass
  leaves the free mass); under Penalty/Lagrange that DOF has an equation and
  counts as free, as before. Guard downstream on `totalFreeMass`, not
  `totalMass`. Learned 2026-10-25 (RSA gap study, gap 6).
