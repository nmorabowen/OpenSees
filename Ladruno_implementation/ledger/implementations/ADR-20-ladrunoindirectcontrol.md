---
wp: ADR-20
title: "LadrunoIndirectControl"
pr: "#184"
status: "shipped — Zone-A IC battery (reduce-to-DisplacementControl + CMOD snap-back)"
section: "table"
legacy_seq: 58
---
| **LadrunoIndirectControl** — indirect / CMOD displacement-control static integrator (classTag 33006), **ADR-20 §8 follow-up #3**. Generalizes stock `DisplacementControl` (a single DOF) to a **weighted multi-DOF control quantity** `ζ = c·U = Σₖ coefₖ·U(nodeₖ,dofₖ)` (e.g. crack-mouth-opening `CMOD = u_A − u_B`), which stays **monotone through snap-back** — where the load factor AND the individual nodal DOFs both reverse — so it follows equilibrium branches that geometric arc-length and single-DOF `DisplacementControl` cannot (de Borst indirect control). Constraint per step `c·ΔU = Δζ`; mirrors `DisplacementControl`'s dUhat/dUbar bookkeeping with the scalar control component replaced by the dot product `c·dUhat`. Parser `integrator LadrunoIndirectControl $incr -dof $node $dof $coef <-dof …> <-iter $numIter $dmin $dmax>`; optional Ramm increment adaptation; no DDM sensitivity (out of scope, as LadrunoArcLength). | Integrator | 33006 | `SRC/analysis/integrator/LadrunoIndirectControl.{cpp,h}`, `tests/test_ladrunoIndirectControl_integrator.py` | shipped — Zone-A IC battery (reduce-to-DisplacementControl + CMOD snap-back) | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
