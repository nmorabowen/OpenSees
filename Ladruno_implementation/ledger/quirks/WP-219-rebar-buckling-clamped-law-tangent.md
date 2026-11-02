---
wp: WP-219
title: "LadrunoRebarBuckling used the law's tangent on the bare stress when g_law was clamped"
pr: "#1056"
date: 2026-11-01
---
### `LadrunoRebarBuckling` tangent differs from the finite difference of its stress for a stocky bar; `reduction` readout reports huge values

- **Bites:** s/d <= 2 only. There `fStar` exceeds the backbone onset stress, so r > 1 and the law sits
  below bare. `g_law` is clamped at 0 (`setTrialStrain`), so sigma = sigma_bare, but `bucklingBranch`
  took `kV1`: 2152.9 against a finite difference of 2001.9, on 19 to 24 states per history. A scan at s/d
  3, 4, 6, 8 and 12 over five histories found no clamped state, so the s/d = 6 recipe never reached it.
- **Readout:** `reduction` (sigma/sigma_bare) divides by a stress crossing zero while a raise is held:
  25251 on a fine step. It enters no response. Below 1 % of fy it now reports 1.
- **Fix (WP-219):** a clamped law passes the bare stress and tangent to the branch.
- **Gate:** `tests/test_ladrunoRebarBuckling_clamp_tangent.py`.
