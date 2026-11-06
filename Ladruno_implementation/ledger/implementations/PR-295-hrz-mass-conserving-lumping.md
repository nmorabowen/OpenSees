---
wp: PR-295
title: "HRZ mass-conserving lumping"
pr: "#295"
status: "shipped — Zone-A 6/6 + standalone"
section: "table"
legacy_seq: 40
---
| **HRZ mass-conserving lumping** ([[35_ladruno_hrz_lumped_mass_adr]]) — `-lump hrz` on all 3 explicit integrators (CDL/ExplicitBathe/LNVD); per-direction `α_d=m_d/S_d` scale, rotational DOFs by mean-α, positive + mass-conserving + rotation-aware; reduces to row-sum on regular elements; u-p/shell non-positive DOFs pass through; g++-verified OpenSees-free kernel | Integrator util | — | `SRC/analysis/integrator/LadrunoMassLumping.h` (`Ladruno::hrzLump`), `CTSLumping::HRZ` branch in `CriticalTimeStep.{h,cpp}`, parsers in `CentralDifferenceLadruno.cpp`/`ExplicitBathe.cpp`/`ExplicitBatheLNVD.cpp`, `tests/test_hrz_lumped_mass.py`, `tests/_hrz_verify/hrz_standalone.cpp` | shipped — Zone-A 6/6 + standalone | [#295](https://github.com/nmorabowen/OpenSees/pull/295) |
