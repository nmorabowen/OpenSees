---
wp: WP-130
title: "A CPPM refusal under a DISCARDING element was committed -- the plain commitState path never read the refusal flag (WP-130, from WP-129's adversarial review)"
legacy_seq: 485
---
### A CPPM refusal under a DISCARDING element was committed -- the plain `commitState` path never read the refusal flag (WP-130, from WP-129's adversarial review)
- **Bites:** SSPquad, stdBrick (= Brick), BbarBrick, SSPbrick, BrickUP and LadrunoSolidShell drop
  `setTrialStrain`'s return code (F7 roster), so a `-cppmOnFail refuse` (or failed ME->CPPM
  fallback) refusal never cuts their step; Newton can converge on the refused, UNINTEGRATED trial
  state and `LadrunoSANISAND::commitState`'s plain (non-IMPL-EX) path committed it -- the refusal
  was honoured by forwarding elements only. WP-129's review found the same hole for the ME cap.
- **Fixed (WP-130, #868):** the plain path checks `mLadrunoCPPMRefused` first: declare to
  `Domain::commit()` (`ladrunoNoteCommitRefusal`, element-independent abort), latch
  (`mImplexCommitRefusedLatch`, reported in `implexRefusals[4]`), restore the trial, return
  `LADRUNO_MATERIAL_REFUSED`. Pinned on SSPquad (`test_cppm_refusal_under_a_discarding_element_does_not_commit`:
  analyze < 0, strain = last committed, further steps refused). WP-129 adds `mSubstepCapHitInME`
  to the same check: whichever merges second ORs the flags. The trial-time latch warning still
  says "-implex companion"; its text predates this second writer.
