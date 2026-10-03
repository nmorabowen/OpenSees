---
wp: LEGACY
title: "TanType 2 under IntScheme 2 is a genuine algorithmic tangent — but any CPPM fallback silently overwrites it with ModifiedEuler's chained tangent, with no diagn…"
legacy_seq: 464
---
### `TanType 2` under `IntScheme 2` is a genuine algorithmic tangent — but any CPPM fallback silently overwrites it with `ModifiedEuler`'s chained tangent, with no diagnostic telling you which one you got
- **Bites:** you ask for `-TanType 2` expecting the consistent (algorithmic) tangent of the CPPM's
  own return map on every step. On any step where the CPPM instead fell back to
  `ModifiedEuler` — which under a global Newton at the campaign increment is most of them, see the
  first entry above — you silently get a different object: `ModifiedEuler`'s substep-chained
  continuum tangent. Nothing distinguishes the two in any response or log.
- **Why (verified on `634824e1f`):** `LadrunoSANISAND3D::getTangent()`
  (`SRC/material/nD/LadrunoSANISAND3D.cpp:170-178`, shadowing the identical
  `ManzariDafalias3D.cpp:134-141`) returns `mCep_Consistent` for `mTangType == 2`. Under scheme 2
  that member is written at `ManzariDafalias.cpp:2609` (`Cep_Consistent = aCepConsistent;`) from
  `NewtonIter2(...)` (`:2457`), which fills it in `NewtonSol` as the condensation of the 19x19 CPPM
  Jacobian (`:3448-3472`, ending `Cep = -1.0 * CSigma;` at `:3472`) — a genuine algorithmic
  tangent, one iterate stale (`NewtonIter2`'s loop tests convergence before the final `NewtonSol`
  call, so the tangent is evaluated at the second-to-last iterate, not the converged one). But
  every path out of the CPPM that reaches `explicit_integrator` — the low-`p` branch
  (`:2431-2436`) or ladder exhaustion (`:2584-2588`) — OVERWRITES `aCepConsistent` with
  `ModifiedEuler`'s own chained product (`:1835`,
  `aCep_Consistent = aCep_thisStep * (aD * aCep_Consistent + T * mIImix)`). `TanType 2` is not
  scheme-2-only: `ModifiedEuler` maintains its own `aCep_Consistent` (`:1490`, `:1835`), so the
  option is meaningful on scheme 1 too (the deck `sanisand_tau0_band.py:349` already uses it) — the
  two objects are different animals, a return-map Jacobian on scheme 2 versus a product of
  continuum tangents over the substep chain on scheme 1, and scheme 2 silently degrades into the
  latter on fallback. **Under `-implex` both are inert regardless:** the material hands out
  `Ce(p_n)` and says so (`LadrunoSANISAND.cpp:2097-2099`, "TanType ... is INERT under -implex").
- **CORRECTION (WP-130, #868): the "genuine algorithmic tangent" has the WRONG SIGN.** See the
  WP-130 entry "IntScheme 2's TanType-2 tangent is MINUS". The reading above (a return-map
  Jacobian, one iterate stale) is right about the object and wrong about its sign.
- **Workaround/status (2026-09-16, WP-105 / F12, no code changed):** not fixed; recorded as a
  read-only finding. If you need to know which tangent a step actually used, cross-reference the
  `substeps` response's `mSubstepsTakenInME` (non-zero iff `ModifiedEuler` ran) alongside
  `TanType 2` output — there is no dedicated flag for it.
