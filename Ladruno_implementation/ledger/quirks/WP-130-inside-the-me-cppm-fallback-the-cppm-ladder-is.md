---
wp: WP-130
title: "Inside the ME -> CPPM fallback the CPPM ladder is guarded by mScheme == INT_BackwardEuler and its explicit exits re-enter ModifiedEuler (WP-130)"
legacy_seq: 486
---
### Inside the ME -> CPPM fallback the CPPM ladder is guarded by `mScheme == INT_BackwardEuler` and its explicit exits re-enter ModifiedEuler (WP-130)
- **Bites:** calling `BackwardEuler_CPPM` from an IntScheme-1 material skips the whole retry ladder
  (vanilla returns an UNCONVERGED state with errFlag 0), and any explicit exit (low-p branch, ladder
  exhaustion) calls `explicit_integrator` -> the ModifiedEuler that just hit its `-maxSubsteps`
  cap, which then hits it again at once (the cap counts per `integrate()`).
- **Fixed (WP-130, #868):** `-meFallback cppm` sets `mLadrunoInMEFallback` for the call: the ladder
  runs (`|| mLadrunoInMEFallback`) up to `-cppmHalvings`, and every explicit exit REFUSES. On
  success the cap flag is cleared; on failure it stays set, so the update is refused exactly as
  without the fallback, and every existing reader of `mSubstepCapHitInME` (the IMPL-EX companion
  included) stays consistent. Not qualified with `-implex` (the parser refuses the combination).
