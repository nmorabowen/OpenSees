---
wp: LEGACY
title: "Mixed-DOF EQ rows (w←(w,θ)) are ALREADY legal end-to-end — but a mortar basis change needs kernel help"
legacy_seq: 144
---
### Mixed-DOF EQ rows (w←(w,θ)) are ALREADY legal end-to-end — but a mortar basis change needs kernel help
- **Bites:** anyone extending a tie/constraint transfer beyond same-DOF rows — the ADR-62 P3.1
  `-hermite` Hermite w–θ shell-edge transfer, and any future shell-to-solid (`θ=½∇×u`) coupling.
- **Why (the good news, verified in source):** `EQ_Constraint` stores per-retained `(node, dof, coef)`
  TRIPLES (`rCoef_i·rDOF_i(rNode_i)`, EQ_Constraint.h:37-45) and `LadrunoProjectionHandler::buildGroups`
  creates one union-find vertex per `(node, dof)` walking arbitrary retained-DOF lists — there is NO
  retained-dof==constrained-dof assumption anywhere. So a slave-translation row that references master
  ROTATION DOFs is an ordinary EQ row: no constraint-class, kernel, or handler change (P3.1 shipped as a
  pure emission-level transform on the collocation path).
- **The catch:** this dodge is COLLOCATION-only for basis changes. `LadrunoMortarKernel::integratePair`
  returns only ACCUMULATED `D/M/g̃` (`PairResult` has no per-GP hook), so putting a different master basis
  (e.g. Hermite functions of `(w, θ_t)`) inside the weak-form `M` integral requires a kernel extension —
  that is why `-hermite -mortar` is a named refusal (mortar-Hermite = deferred P3.1b).
- **Bonus oracle gotcha:** when testing a Hermite w-row against slope-inconsistent (Mindlin) data, do NOT
  sample slaves at edge midpoints — the shear error enters through `H2+H4 ∝ ξ−3ξ²+2ξ³`, which has a root
  at exactly ξ=½, silently hiding the Kirchhoff-assumption error (proto_p3_1_hermite_tie.py T6).
