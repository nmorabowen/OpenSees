---
title: ADR-97 P1 — mutation gate record
project: Ladruno
status: complete
owner: nmora
tags:
  - implementation
  - material
  - verification
---

# ADR-97 P1 — mutation gate (gate 5)

ADR-87 D2: a green suite is only evidence if deleting the physics turns it red.
This is the record for `wp/97b-cp-smooth` (PR
[#819](https://github.com/nmorabowen/OpenSees/pull/819)), baseline build
`ec6091c4f`.

## What was dropped

ONE term of the consistent tangent — the `dl * dm/dσ` term of the **algorithmic
elastic modulus**

```
Xi = ( E^-1 + dl * dm/dsigma )^-1
```

removed from the closest-point Jacobian's `J_ss` block in
`SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h::cp_assemble`:

```diff
-                Jm(i, j) = ((i == j) ? 1.0 : 0.0) + dl * Edmds(i, j);
+                Jm(i, j) = ((i == j) ? 1.0 : 0.0); // ADR-97 MUTATION: dl*dm/ds dropped
```

Applied and reverted by `Ladruno_implementation/adr97_scripts/mutate_p1.py`
(`apply` / `revert`) on a scratch build never pushed; only the `OpenSeesPy`
target was rebuilt (one translation unit + link) and the resulting
`OpenSeesPy.dll` copied over `dist/bin/opensees.pyd`.

**Why this term.** With `dl * dm/dσ` gone, `Xi` collapses back to `E` and the
consistent tangent degenerates into *exactly* the shipped `Continuum` operator
`E − (E:m)(n:E)/(n:E:m − H)`. That equivalence is the whole finding of ADR-94
M3: the 57 % `Continuum` error **is** the size of this one term at realistic
step sizes. So the mutation does not merely perturb the new tangent — it turns
it back into the old one, which is the sharpest available test that the new code
does something the old code did not.

**Why the residual is untouched.** The mutation changes only the JACOBIAN. An
inexact Newton still converges to the exact residual, so the committed stress is
unaffected by construction and gates 1 (stress), 3, 4 and 6 must stay GREEN. A
mutation that killed everything would prove nothing about *which* test covers
the tangent.

## Result — 4 killed

| test | baseline | mutated | verdict |
|---|---|---|---|
| `test_gate2_algorithmic_matches_a_finite_difference_of_the_committed_map` | rel_err **1.51e-11** | rel_err **0.573447** | **KILLED** |
| `test_gate2_algorithmic_on_a_sheared_twelve_dof_rig` | rel_err **2.15e-10** | rel_err **0.670133** | **KILLED** |
| `test_gate2_two_cube_newton_cost` | testIter 4,4,4,4 (16 total, 7.1× vs `Continuum`) | testIter 6,6,6,6 (24 total, 4.7×) | **KILLED** |
| `test_gate1_newton_converges_in_at_most_five_iterations` | AF **3** iterations | AF **13** iterations | **KILLED** |

The two finite-difference numbers are the point of the whole exercise: the
mutated `Algorithmic` lands on **0.573447** and **0.670133**, which are — to
every digit printed — the `Backward_Euler`/`Continuum` values the P0 oracle
measured for the same two rigs (`adr97_oracle/README.md`, element block). The
mutation does not make the tangent *wrong in some direction*; it makes it
`Continuum`, exactly as the algebra says it must.

The fourth kill was not predicted and is a bonus: dropping the term also
destroys the LOCAL Newton's quadratic convergence, so the Armstrong–Frederick
return goes 3 → 13 iterations. Perfect and linear-hardening von Mises stay at 1
iteration even mutated, because their return is exactly radial and `dl = 0` at
the elastic-predictor start, so the dropped term is zero on the only Jacobian
those cases ever form.

## Result — survivors, and why each one *should* survive

39 of the 43 cases still passed. Every one of them is a case the mutation
provably cannot reach:

* **gate 1, committed stress (13 cases).** Unchanged to the last digit
  (7.4e-14 … 2.4e-12; the AF rotating-normal case moves 9.17e-9 → 3.14e-9, both
  inside the yield tolerance). This is the design of the mutation, not a hole:
  the residual is exact, so the root is exact.
* **gate 3 (path error).** `Closest_Point` 5.6795e-3 vs `Backward_Euler`
  2.3127e-2, ratio 4.07× — identical to baseline. Gate 3 measures the MAP, and
  the map did not change.
* **gate 4 (10 cases).** Byte identity is about `Backward_Euler`, which the
  mutation does not touch; the CP≡BE agreement and CP≠BE divergence numbers are
  properties of the map.
* **gate 6 (13 cases).** Refusals and parser gates.

That split is the useful output: **gate 2 is the only gate that covers the
consistent tangent**, and it does so on two independent rigs plus an iteration
count. If gate 2 were ever deleted or weakened, nothing else in this suite would
notice a tangent regression.

## Reverted

`mutate_p1.py revert`, rebuild, re-measure: gate 2 back to **1.51e-11** /
**2.15e-10** and the full ADR-97 suite back to 43/43 — identical to the
pre-mutation baseline, which is also the proof that the mutation build and the
revert build differed in nothing else. Full battery (ADR-84 + ADR-94 + ASDP +
ADR-97, 18 files): **134 passed, 2 skipped, 0 failed**.
