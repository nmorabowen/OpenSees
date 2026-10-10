---
wp: LEGACY
title: "Stage 2 PR-2b (2D objectivity/invariance/robustness gates) SHIPPED (test-only): the deferred Stage-2 validation gates,…"
section: "history"
legacy_seq: 121
---
- **Stage 2 PR-2b (2D objectivity/invariance/robustness gates) SHIPPED (test-only):** the deferred Stage-2 validation gates, verified by probing first (the PR-2a hinge already passes them). `tests/test_ladrunoDispBeamColumn2d_hinge.py` (8→15): corotational large-rotation `Gf`-dissipation at 74° tip rotation (pinned invariant under finite rotation), orientation invariance (0° vs 90° → identical M–θ to round-off), integration-objectivity nIP sweep (Lobatto 2..6 → invariant, no residual nIP drift), solver robustness (Newton/ModifiedNewton/NewtonLineSearch/KrylovNewton all dissipate `Gf`). The non-Newton "stale-α" hole does NOT bite (residual always post-`update()`), so that hardening is unneeded. **2D Stage-2 is gate-complete** except the `ladruno_drive` collapse test (blocked on RESERVED dissipation arc-length) and `-nl`+hinge.
