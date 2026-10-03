---
wp: WP-128
title: "test_ladruno_sanisand_flip_determinism.py's two push-step tests fail on ladruno for a reason unrelated to MKL or the machine — their pinned constants predate a…"
legacy_seq: 500
---
### `test_ladruno_sanisand_flip_determinism.py`'s two push-step tests fail on `ladruno` for a reason unrelated to MKL or the machine — their pinned constants predate a numerical change (WP-128 note)
- **Bites:** `test_first_ten_push_steps_bit_identical_across_mkl_threads` and `test_default_first_step_is_immune_to_the_hold_lottery` fail at the first push step with `-3`. They fail on the unmodified `c03a1bd4b` too, so they look like a regression of whatever you just changed. The step fails identically under `system FullGeneral`, so it is not Pardiso/MKL threading. At `NormDispIncr 1e-6` it converges to 9.731325 kN/m, where the docstring pins 9.659111 (`48c0e99bc`). The deck (TanType 2) moved after that build; candidates in the range are WP-110's tangent fix `dee04dbe3` and WP-112. Not bisected.
- **Workaround/status:** treat as pre-existing until someone re-derives the pins and the 1e-8/100-iteration convergence on the current tree (`Ladruno_files/testbed/sanisand_ring_trace/flipdet_probe.py` reproduces it in ~2 min). **→ Resolved by WP-136 (row below): WP-110's tangent fix `dee04dbe3` moved the Newton path; the push now runs `FixedNumIter`.**
