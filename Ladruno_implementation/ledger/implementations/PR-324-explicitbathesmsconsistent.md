---
wp: PR-324
title: "ExplicitBatheSMSConsistent"
pr: "#324"
status: "shipped — Zone-A (supra-stable + no-nodal-mass, f1-preservation vs lumped FFT,…"
section: "table"
legacy_seq: 45
---
| **ExplicitBatheSMSConsistent** ⟨**W1-E2 #419: COLLAPSED → `ExplicitBathe -sms -consistent`; `.cpp/.h` DELETED, tag 33010 a deprecated alias**⟩ ([[38_ladruno_consistent_mass_scaling_adr]]) — CONSISTENT (Olovsson) mass scaling on the Noh-Bathe `ExplicitBathe` (sibling of CentralDifferenceSMSConsistent): centroidal `M̄=β[diag(m)−mmᵀ/Mₑ]` + matrix-free PCG (`consistentPCG`) refined at BOTH Noh-Bathe sub-step solves via a new no-op `refineAccel()` hook on `ExplicitBathe` (default no-op ⇒ ExplicitBathe AND lumped ExplicitBatheSMS byte-identical). Reuses `buildMassScalingConsistent` (slave+master MP exclusion). f1 −0.17% vs lumped −53% | Integrator | 33010 | `SRC/analysis/integrator/ExplicitBatheSMSConsistent.{cpp,h}`, `LadrunoMassScaling.h`, `ExplicitBathe` `refineAccel` hook, `tests/test_explicitBatheSMS_integrator.py` | shipped — Zone-A (supra-stable + no-nodal-mass, f1-preservation vs lumped FFT, reduce-to-base bit-identical, PCG≤30 at both sub-steps) | [#324](https://github.com/nmorabowen/OpenSees/pull/324) |
