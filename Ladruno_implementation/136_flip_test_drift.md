# WP-136 — why `test_ladruno_sanisand_flip_determinism.py`'s push-step pins drifted

Status: IN PROGRESS (draft PR). Investigation of the two failing tests
`test_first_ten_push_steps_bit_identical_across_mkl_threads` and
`test_default_first_step_is_immune_to_the_hold_lottery`, handed over by WP-128
(`128_sanisand_ring_trace.md` §7 on `wp/128-sanisand-ring-trace`, draft #869).

Plan: reproduce on `ladruno`; bisect `48c0e99bc..ladruno` over the SANISAND /
ManzariDafalias / deck-path commits, suspect `dee04dbe3` (WP-110 tangent fix)
first; classify fix-consequence vs regression; extend the "stale pins" quirk row.
