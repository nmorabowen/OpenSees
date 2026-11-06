---
wp: PR-324
title: "ExplicitBatheLNVDSMSConsistent"
pr: "#324"
status: "shipped — Zone-A (supra-stable + no-nodal-mass, f1-preservation FFT, reduce-to-…"
section: "table"
legacy_seq: 47
---
| **ExplicitBatheLNVDSMSConsistent** ⟨**W1-E2 #419: COLLAPSED → `ExplicitBathe -lnvd -sms -consistent`; `.cpp/.h` DELETED, tag 33012 a deprecated alias**⟩ ([[38_ladruno_consistent_mass_scaling_adr]]) — CONSISTENT (Olovsson) mass scaling on `ExplicitBatheLNVD` (sibling of ExplicitBatheSMSConsistent): centroidal `M̄` + matrix-free PCG at BOTH Noh-Bathe sub-step solves via a new no-op `refineAccel()` hook on `ExplicitBatheLNVD` (default no-op ⇒ LNVD + lumped LNVDSMS byte-identical). Reuses `buildMassScalingConsistent` (slave+master MP exclusion). f1 −0.17% vs lumped −53% | Integrator | 33012 | `SRC/analysis/integrator/ExplicitBatheLNVDSMSConsistent.{cpp,h}`, `LadrunoMassScaling.h`, `ExplicitBatheLNVD` `refineAccel` hook, `tests/test_explicitBatheLNVDSMS_integrator.py` | shipped — Zone-A (supra-stable + no-nodal-mass, f1-preservation FFT, reduce-to-base bit-identical, PCG≤30 both sub-steps) | [#324](https://github.com/nmorabowen/OpenSees/pull/324) |
