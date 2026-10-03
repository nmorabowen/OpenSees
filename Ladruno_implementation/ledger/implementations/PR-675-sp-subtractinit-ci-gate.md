---
wp: PR-675
title: "sp -subtractInit CI gate"
pr: "#675"
status: "shipped"
section: "table"
legacy_seq: 105
---
| **`sp -subtractInit` CI gate** — 5 pytest cases locking in the openseespy `-subtractInit` fix (see [[LEDGER_vanilla_files]]): stage-1 premise (guards against a vacuous pass if `u0` were 0), no-flag-is-absolute, flag-is-incremental (`DELTA + u0` per `TransformationDOF_Group::enforceSPs`'s `setTrialDisp(value + initial)`), **flag-actually-changes-the-answer** (stated as a DIFFERENCE, not an absolute, so it survives any future revisit of the subtraction's sign convention — whatever the flag means, it must not mean *nothing*), and a case documenting that the flag is **inert under `constraints Plain`** by design so that is never mistaken for the bug returning. **Verified as a real gate, not decoration: 2 of the 5 FAIL on the pre-fix binary** (`assert 0.0 > 1e-09`) and all 5 pass after. Pure upstream `Truss` + `Elastic`, no MKL and no fork elements, so **Zone-A runs it** — unlike the PARDISO battery, which is Windows/MKL-gated. | Tests (pytest, zone_a) | — | `tests/test_sp_subtract_init.py` | shipped | [#675](https://github.com/nmorabowen/OpenSees/pull/675) |
