---
wp: LEGACY
title: "setStrain (testUniaxialMaterial) COMMITS — central FD tangent is invalid for plasticity"
legacy_seq: 57
---
## `setStrain` (testUniaxialMaterial) COMMITS — central FD tangent is invalid for plasticity

The Python `setStrain(eps)` command (`SRC/interpreter/OpenSeesCommandsPython.cpp`,
`ops_setStrain`) calls **both** `setTrialStrain(eps)` **and** `commitState()`. So
the `_testbed.fem_checks.uniaxial_tangent_fd` central difference — which probes
`strain0-d`, `strain0`, `strain0+d` expecting all three from one committed state —
straddles an **elastic unload** on the minus side for any path-dependent (plastic)
material, returning ~`(E + E_alg)/2` instead of the consistent tangent (caught on
LadrunoUniaxialJ2: analytic 176.7 vs the broken FD 597 ≈ (1000+176)/2). The helper
is only valid for nonlinear-elastic materials.

- **Rule:** to FD a *plasticity* consistent tangent, probe each point as an
  INDEPENDENT one-step return from a FRESH material (`wipe` → redefine →
  `testUniaxialMaterial` → `setStrain`), so all probes are one-step-from-zero; the
  central difference of that one-step map IS the algorithmic tangent (O(d²)). See
  `test_consistent_tangent_fd` (V6) in `tests/test_ladrunoUniaxialJ2_material.py`.
  Also note `testUniaxialMaterial(tag)` returns the SAME object (not a copy), so it
  does not reset committed state — you must rebuild the material to get a clean one.
  Learned 2026-06-02 (LadrunoUniaxialJ2 polish).
