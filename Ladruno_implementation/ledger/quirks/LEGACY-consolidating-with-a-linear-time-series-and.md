---
wp: LEGACY
title: "Consolidating with a Linear time series and then adding a deviatoric pattern applies the deviator AT FULL AMPLITUDE on the first increment"
legacy_seq: 258
---
### Consolidating with a `Linear` time series and then adding a deviatoric pattern applies the deviator AT FULL AMPLITUDE on the first increment
- **Bites:** every single-element stress-path probe built as "ramp the
  confinement to lambda = 1, then add a second pattern and keep stepping". Both
  patterns scale with the SAME load factor, so the first post-consolidation
  increment (lambda = 1.0025) evaluates the deviatoric pattern at 1.0025x its
  full amplitude, not at 0.0025x. The probe then reports failure at step 0 with
  ZERO measured deviator — `cone_probe.py` avoided this with `Path` series keyed
  to pseudo-time and it is easy to miss when porting the same probe to a
  different material.
- **Tell:** every path of a stress probe returns the consolidation state
  unchanged (`sqrt(J2) = 0`, `alpha = 0`) and "failure at load factor ~1.00".
- **Rule:** `ops.loadConst("-time", 0.0)` after the consolidation stage, so the
  held pattern stops scaling and the new pattern ramps from zero. Same rule as
  the push stage of any displacement-controlled runner.
  *2026-07-30 (ADR-79 collapse study).*
