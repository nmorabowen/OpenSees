---
wp: WP-130
title: "-cppmStart explicit clobbered the ladder's errFlag, and a successful guess accepted a coarse one-step BE root (WP-130 review r1)"
legacy_seq: 517
---
### `-cppmStart explicit` clobbered the ladder's errFlag, and a successful guess accepted a coarse one-step BE root (WP-130 review r1)
- **Bites:** the guess's `NewtonIter2` wrote `errFlag` before `if (errFlag == -1) SchemeControl = 3`,
  so trial non-convergence (0: halve) followed by a singular guess (-1) skipped halving and went
  straight to explicit/refuse (17 of 228 guess tries in the review's set). And an accepted guess
  root is one backward-Euler step over an increment the ladder would have halved: errors up to 0.77
  relative vs an oracle (`Check()` cannot see it: its p < 0 test is commented out).
- **Fixed / status (WP-130, #868):** the trial's errFlag and state are restored when the guess is
  not accepted; `gZ` joins the NaN check; the root is accepted only if dGamma >= 0, p > 0 and it
  agrees with the 50-substep explicit walk to 2 % (`LADRUNO_GUESS_AGREE`). Measured on the same
  300-increment oracle set: the error is more than 2x the default ladder's on 105 of 171 guess
  increments (74 with an extra 0.01 absolute margin; recount by confirmation review round 2);
  worst oracle-converged increment 0.54 vs 0.51; largest single degradation 0.09 -> 0.40 --
  (median 0.021 vs 0.007) -- so `-cppmStart explicit` is REMOVED from the recommended recipe and
  documented as a speed-for-accuracy trade.
