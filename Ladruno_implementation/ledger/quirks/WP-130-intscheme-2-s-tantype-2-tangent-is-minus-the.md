---
wp: WP-130
title: "IntScheme 2's TanType 2 tangent is MINUS the derivative of its own return map -- the element gets a negative-definite stiffness (WP-130)"
legacy_seq: 514
---
### `IntScheme 2`'s `TanType 2` tangent is MINUS the derivative of its own return map -- the element gets a negative-definite stiffness (WP-130)
- **Bites:** `ManzariDafalias::NewtonSol` ends `Cep = -1.0 * CSigma;`. The condensed CPPM system
  is `DSigma d_sigma = d_eps` (R1 = eStrain - TrialElasticStrain + ..., so d R1/d eps = -I), and
  `CSigma = (aC DSigma)^-1 aC = DSigma^-1`: the algorithmic tangent is `+CSigma`. The local Newton
  never uses `Cep` (its update is `delSig = -CSigma SConstant`), so the RETURN is right and only the
  tangent handed to the element is flipped -- a negative-definite 6x6 under `TanType 2`. Measured
  (`Ladruno_files/testbed/hypo_bearing/wp130_f18c/q_tangent_fd.py`, 3D cube, a plastic step whose
  replay reproduces the analysis to 4e-15): `||-T - D_fd|| / ||D_fd|| = 1.24e-3` (the one-iterate
  staleness), `||T - D_fd|| / ||D_fd|| = 2.0`; the diagonal of T is -1.66e4 where D_fd has +1.66e4.
- **Why it hid:** WP-105 (F12 §5.1) READ the code and called it "a genuine algorithmic tangent";
  nobody compared it with a finite difference. Every zero-free-DOF material-point deck converges in
  one global iteration whatever the tangent (no free equations), and every free-DOF scheme-2 failure
  was blamed on the recursive-halving ladder, which it also has. On F12's bearing deck the global
  Newton with the vanilla tangent DIVERGES from its first iteration (unbalance x3-5 per iteration),
  and only KrylovNewton's relaxed rung commits steps -- that, not the ladder alone, is why F12 found
  scheme 2 "475x shallower". Checklist rule, again: verify a tangent against a finite difference
  with free equations (`ladruno-new-material`, "Verifying a tangent").
- **Workaround/status (WP-130, #868):** `-cppmTangent fixed` hands out `+CSigma`
  (`tests/test_ladruno_sanisand_cppm_newton.py` pins both: vanilla -T within 1e-2 of D_fd, fixed
  +T within 1e-2, 3D and -- by a whole-run FD, sigma_33 not being exposed -- the plane-strain
  wrapper). **Owner decision (WP-130): `fixed` is the LadrunoSANISAND DEFAULT**; `-cppmTangent
  vanilla` reproduces the old binary; vanilla `ManzariDafalias` keeps the wrong sign (upstream
  report deferred). F12's bearing deck (x10z8, `h1.0_e0.6944`, 1200 s budget, TanType 2, driver unchanged), the RECOMMENDED recipe (`fixed` default + `-cppmOnFail refuse -cppmHalvings 3 -cppmLineSearch on`, no `-cppmStart`; build 6726f5e24, `wp130_f18c/tables_recipe.md`), measured back to back with an IntScheme-1 control on the same box, which was at 100 % CPU (so compare these two with each other only): the recipe reaches s/B 0.00293 / 0.00421 / 0.00523 / 0.00626 at 300 / 600 / 900 / 1200 s against IntScheme 1's 0.00138 / 0.00250 / 0.00442 / 0.00698 -- ahead at 300, 600 and 900 s, BEHIND at 1200 s -- with 3.2 global iterations per committed step against 16.8 (481 of 527 steps on the plain Newton rung), load-settlement within 0.2-1.4 % of IntScheme 1, and its refusals had spent 76 of the driver's 80 pinned subdivisions when the wall stopped it. Vanilla IntScheme 2 reaches 0.00002 and `-cppmTangent fixed` alone 0.00378 (earlier, unloaded runs). The global Newton is NOT quadratic even with the fixed tangent: median observed order 1.14 on the last three residuals (9 % of committed calls >= 1.8); see the four tangent error sources. (Pre-round-1 note, superseded: an arm WITH `-cppmStart explicit` -- whose accuracy review round 1 measured and rejected -- reached 0.00876 in 1081 s on an unloaded box.).
