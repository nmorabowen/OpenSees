---
wp: LEGACY
title: "ManzariDafalias::ForwardEuler (IntScheme 5) has a shadowed Vector r — r is identically ZERO, both (n:r) terms are silently dropped"
legacy_seq: 365
---
## `ManzariDafalias::ForwardEuler` (`IntScheme 5`) has a shadowed `Vector r` — `r` is identically ZERO, both `(n:r)` terms are silently dropped

**Found 2026-09-05, ADR-92 P0 source extraction, verified at source.**

`:1330-1332` reads `Vector r(6); if (p > small) Vector r = GetDevPart(CurStress) / p;` — the
inner declaration is block-scoped, dies at the brace, and the outer `r` used at `:1334`,
`:1337`, `:1342`, `:1350` stays zero. `Kp`, `temp4` and `NextDGamma` lose their volumetric
coupling. The P0 oracle, reproducing it verbatim, puts scheme 5's terminal `q/p` at **1.188
against 0.946** for scheme 1 on the same path — a silent 26 % on mobilised strength.

- **Scope, measured not assumed:** `ModifiedEuler` (`:1483`) and `GetElastoPlasticTangent`
  (`:4839`) both use the corrected two-statement form `r = GetDevPart(...); r /= p;` with the
  one-liner commented out directly above — the trap was hit and fixed twice and missed once.
  `explicit_integrator`'s `default:` is `ModifiedEuler` (`:1060`), so scheme 2's fallback is
  clean too. **Reaches schemes 5, 9 and 4 only. Not the deck default; not the TIMs campaign.**
- **Fix owed (vanilla, two lines, separate PR):** `r = GetDevPart(CurStress); r /= p;` inside
  the `if`. Until then the entry above on schemes 3/5 having no error control has a second
  reason not to use 5.
- **✅ FIXED — WP-158, [#901](https://github.com/nmorabowen/OpenSees/pull/901) (2026-10-01).**
  `r = GetDevPart(CurStress) / p;` (`:1586` at d63f49750 + the fix). Upstream master
  (316cb2dbc, 2024-05-01) still has the shadow at `:1233`. Three things the earlier note did
  not know, all measured on the d63f49750 build vs the WP-158 build
  (`tests/test_manzari_forward_euler_r.py`, its `__main__`):
  - **Two more defects in the same function, tangent only.** `temp2 = 2G n - (n:r) I` was
    missing `K` (the multiplier's numerator is `2G n:de - K de_v (n:r)`), hidden while
    `r == 0`; fixing `r` alone would have left the scheme-5 tangent inconsistent with its own
    stress update. `temp1 = 2G mIIdevMix + K mIIvol` put **2G, not G, on the shear diagonal**
    (the mixed-variant identity; the stress update answers an engineering shear strain with
    G) — the WP-110 F15 family, never probed for scheme 5. Both fixed (`temp1 = aC`). FD
    gate: max |Ct − Cfd| / max |Cfd| = 0.56 before, 0.64 with only `r` + `K`, 1.2e-11 after.
  - **It is NOT a model error that survives refinement.** On a drained triaxial (p_cell
    100 kPa, 2 % axial, one SSPbrick) the OLD scheme 5 converges onto scheme 1 too: gap in
    (q/p, eps_v) at 3200 steps is (6.5e-4, 5.3e-4). The wrong multiplier drifts the stress
    off the yield surface; when it lands inside, `explicit_integrator` re-intersects it on
    the next step, which enforces consistency geometrically. What the bug broke is the
    ORDER of the step: Delta f/p of ONE step from an on-surface state shrinks 3.75x per 4x
    shorter step (first order) before, 15.9x (second order) after. At practical steps the
    error is large and erratic: 800 steps (1.7e-2, 5.3e-2) before vs (3.3e-3, 1.8e-3) after;
    50 steps, Newton diverges at eps_a 0.32 % before. The P0 oracle's 26 % was one coarse
    path. **Test lesson:** a comparison against a reference at a fine step can be green
    WITH this bug; test the order of one step (the consistency drift) instead.
  - **Reach, corrected:** scheme 5; scheme 4 only on increments where `MaxEnergyInc`'s
    energy test does not fire; schemes 7, 8, 9 only on increments ≤ `maxStrainInc` (1e-5) —
    above it `MaxStrainInc` hands FE uninitialised moduli (next entry); and the opt-in
    WP-130 `-cppmStart` guess walk. Of the WP-129 byte-identity decks only `ls3d_s5` moved
    (re-pinned in #901); `ls3d_s7/8/9` are byte-identical across the fix.
