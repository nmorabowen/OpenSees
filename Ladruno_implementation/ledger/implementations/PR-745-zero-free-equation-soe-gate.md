---
wp: PR-745
title: "Zero-free-equation SOE gate"
pr: "#745"
status: "shipped"
section: "table"
legacy_seq: 138
---
| **Zero-free-equation SOE gate** — 13 pytest cases locking in the dense-SOE `size==0` fix (see [[LEDGER_vanilla_files]], 6 rows): six systems (`FullGeneral`, `BandGeneral`, `BandSPD`, `ProfileSPD`, `SProfileSPD`, plus `UmfPack` as the always-worked control) x two scenarios, plus `Diagonal` on the fresh scenario only — *fresh* (a single stdBrick unit cube with 1/8-symmetry fixes and every remaining DOF driven by an `sp` under `constraints Transformation`, i.e. zero equations at the first `setSize()`) and *shrink* (node 7 free first, 3 equations, then `sp`'d in a new pattern so the SAME SOE is resized 3 -> 0). Each case runs in a **subprocess**: the pre-fix failure is `exit(-1)` inside the SOE, so a regression kills the interpreter and would otherwise take the whole pytest session down instead of failing one case; the child's nonzero exit code is the assertion. **Verified as a real gate: 6 of 13 FAIL on the pre-fix binary** — and stated honestly, those 6 are ALL the *fresh* cases (the 7th, UmfPack, is the control). **The shrink half does NOT reproduce the crash and passes pre-fix**: a 3 -> 0 resize has `size != oldSize`, so the wrapper block runs and rebuilds them at zero length — the null-wrapper path is reachable only on the FIRST `setSize`. Shrink is kept as a non-crashing regression guard on the quieter second defect (the ProfileSPD `iDiagLoc[-1]` UB read), which has no other coverage. All 13 pass after. Upstream `stdBrick` + `ElasticIsotropic` only, no MKL and no fork objects, so Zone-A runs it. | Tests (pytest, zone_a) | — | `tests/test_soe_zero_free_equations.py` | shipped | [#745](https://github.com/nmorabowen/OpenSees/pull/745) |
