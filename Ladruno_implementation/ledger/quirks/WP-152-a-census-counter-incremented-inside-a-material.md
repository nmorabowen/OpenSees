---
wp: WP-152
title: "A census counter incremented INSIDE a material's update counts update CALLS, not steps — an element may update a point several times per step, each call starti…"
legacy_seq: 535
---
### A census counter incremented INSIDE a material's update counts update CALLS, not steps — an element may update a point several times per step, each call starting from the committed state (WP-152)
- **Bites:** a "transition" counter (an entry, an exit, a re-seat) bumped in `setTrialStrain`/`integrate()`.
  - An element can call the material update more than once in one step: in the trial phase, again before the residual, and on every Newton iterate.
  - Each call recomputes the trial from the committed state and meets the same transition again.
  - Measured: WP-152's first build counted one re-contact as 2 on a single stdBrick path (`sepExits`).
  - A step that is cut and retried also counts transitions that never happened.
- **Rule:** Count a STATE TRANSITION when it commits. The update records the trial's transition in a trial-only field; `commitState()` counts it and clears the field; `revertToLastCommit()` clears it too. Per-call counters (substeps, rejections) are fine, but document them as per-call (cf. WP-151's `reseatHeld`: decisions, not reversals).
- **Workaround/status:** ✅ WP-152 (`sepEvent`, counted in `LadrunoSANISAND::commitState`). Test: `tests/test_ladruno_sanisand_tension_cutoff.py` asserts +1 per committed transition. [[152_sanisand_tension_cutoff]].
