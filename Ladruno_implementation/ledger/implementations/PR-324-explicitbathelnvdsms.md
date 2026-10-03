---
wp: PR-324
title: "ExplicitBatheLNVDSMS"
pr: "#324"
status: "shipped — Zone-A (supra-stable, nodal-mass injection); base LNVD regression byt…"
section: "table"
legacy_seq: 46
---
| **ExplicitBatheLNVDSMS** ⟨**W1-E2 #419: COLLAPSED → `ExplicitBathe -lnvd -sms`; `.cpp/.h` DELETED, tag 33011 a deprecated alias**⟩ ([[36_ladruno_selective_mass_scaling_adr]]) — LUMPED selective mass scaling on the Noh-Bathe + FLAC-LNVD `ExplicitBatheLNVD` (sibling of ExplicitBatheSMS): protected `ExplicitBatheLNVD` classTag ctor; `domainChanged` injects nodal mass via `Ladruno::buildMassScaling` + restore. No solve hook (mass on RHS). FLAC `alpha` separate from sizing. Command takes `$p $alpha` first | Integrator | 33011 | `SRC/analysis/integrator/ExplicitBatheLNVDSMS.{cpp,h}`, `LadrunoMassScaling.h`, `ExplicitBatheLNVD` protected ctor, `tests/test_explicitBatheLNVDSMS_integrator.py` | shipped — Zone-A (supra-stable, nodal-mass injection); base LNVD regression byte-identical | [#324](https://github.com/nmorabowen/OpenSees/pull/324) |
