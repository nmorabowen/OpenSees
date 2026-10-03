---
wp: ADR-36
title: "ExplicitBatheSMS"
pr: "#324"
status: "shipped — Zone-A (supra-stable necessity, nodal-mass injection); base ExplicitB…"
section: "table"
legacy_seq: 44
---
| **ExplicitBatheSMS** ⟨**W1-E2 #419: COLLAPSED → `ExplicitBathe -sms`; `.cpp/.h` DELETED, tag 33009 a deprecated alias**⟩ ([[36_ladruno_selective_mass_scaling_adr]]) — LUMPED selective mass scaling on the Noh-Bathe `ExplicitBathe` (sibling of CentralDifferenceSMS): subclass via a new protected `ExplicitBathe` classTag ctor; `domainChanged` injects additive nodal mass via the shared `Ladruno::buildMassScaling` + restore lifecycle. ExplicitBathe assembles only the mass on the RHS, so the nodal injection is seen with NO solve-path change (the ADR-36 "trivial follow-up"). Command takes the Noh-Bathe `$p` first | Integrator | 33009 | `SRC/analysis/integrator/ExplicitBatheSMS.{cpp,h}`, `LadrunoMassScaling.h`, `ExplicitBathe` protected ctor, `tests/test_explicitBatheSMS_integrator.py` | shipped — Zone-A (supra-stable necessity, nodal-mass injection); base ExplicitBathe regression byte-identical | [#324](https://github.com/nmorabowen/OpenSees/pull/324) |
