---
wp: LEGACY
title: "Energy-balance check script energy_check.py was stale vs the chunked DATA layout (fixed)"
legacy_seq: 45
---
### Energy-balance check script `energy_check.py` was stale vs the chunked `DATA` layout (fixed)
- **Bites:** `energy_check.py` failed with `TypeError: Only 1D arrays allowed for
  fancy indexing` — its `_read_result` iterated `grp["DATA"]` as a per-step *group*,
  but the recorder writes `ON_DOMAIN/ON_REGIONS energyBalance DATA` as a chunked
  `[T×nrows×ncomp]` **dataset** (the standardized streaming layout; recorder output
  is correct).
- **Fix:** read via `lf.iter_step_slices(grp)` (the canonical slicer the parity
  checks already use; tolerates both chunked and legacy `DATA/STEP_k` layouts).
  After the fix the energy kernel matches the EnergyBalance text sidecar to ~5e-9.
  Learned 2026-05-31.
