---
wp: WP-144
title: "The initial state of LadrunoNorSand can be PROJECTED by the floor, and F0 is then evaluated at the PROJECTED stress; -pi0_auto is the unified rule (S.53), the…"
date: 2026-11-06
---
### The initial state of LadrunoNorSand can be PROJECTED by the floor, and `F0` is then evaluated at the PROJECTED stress; `-pi0_auto` is the unified rule (S.53), the refusal hint is the LEGACY rule (WP-144 round 3, sheet 5.4 / 9.7)
- **Bites:** (1) `initialState` floors the inverted `sigma0` (`n_f_init = 1`): the point then sits at `sigma(eps^e_f)` (`p = -p_min`; under HAR `q` also changes), not at the deck's `sigma0`. Evaluating `F(sigma0, pi_i0)` at the deck stress would classify a start the material does not have. The shell evaluates `F0` at `stress(s0)` when `n_f_init != 0`. (2) The unified rule (the surface through `(p_init, max(eta_init, c2 M))`, `c2 := 0` without a cap) is NOT the surface through `sigma0` when `eta_init < c2 M`: it is a larger surface. The refusal hint "give `-pi0 <=` the surface through sigma0" must use the pre-round-3 rule (`PI0_LEGACY`), or it would print a value that is not on the surface. (3) A `-pi0` on the deck overrides the rule; `-pi0_auto` asks for it; giving neither is still a refusal (no silent on-surface start), giving both is refused (207).
- **Rule:** after anything that can move the initial state (the floor), recompute every quantity derived from it from the state, not from the deck inputs.
- **Workaround/status:** IMPLEMENTED (WP-144 round 3, `LadrunoNorSand::buildInitialState`; the resolved `pi_i0` replaces the NaN so `recvSelf` rebuilds the same state). Final gate pending.
