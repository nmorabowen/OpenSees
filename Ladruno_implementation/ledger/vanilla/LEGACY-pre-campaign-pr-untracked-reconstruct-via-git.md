---
wp: LEGACY
title: "(pre-campaign; PR untracked — reconstruct via git log -- SRC/material/nD/ASDPlasticMaterial3D) -- 1 vanilla row(s)"
files: ["`SRC/material/nD/ASDPlasticMaterial3D/` (10 files: `ASD_material_definitions.cpp`, `All*` registries, `HardeningFunction.h`, `OPS_AllASDPlasticMaterial3Ds.cpp`, `gen_ASD_material_definitions_CPP.py`, `CMakeLists.txt`, `ASDPlasticMaterial3D.h`)"]
table: "main"
legacy_seq: [311]
---
| `SRC/material/nD/ASDPlasticMaterial3D/` (10 files: `ASD_material_definitions.cpp`, `All*` registries, `HardeningFunction.h`, `OPS_AllASDPlasticMaterial3Ds.cpp`, `gen_ASD_material_definitions_CPP.py`, `CMakeLists.txt`, `ASDPlasticMaterial3D.h`) | markers added under ADR-94 (R0 debt, 2026-09-07): extend jaabell's ASDPlasticMaterial3D template framework with **Hoek–Brown** (rock) and **StiffSoil shear+cap** model components — new YF/PF/EL/hardening headers (listed in [[LEDGER_implementations]] territory: `HoekBrown_YF/PF/Utils/ParameterTypes`, `StiffSoilShear_YF/PF`, `StiffSoilCap_YF/PF`, `StiffSoil_EL`, `StiffSoil_HardeningFunctions`, `test_HoekBrown.cpp`) + regenerated registry/definitions files (+2,944/−94 over 22 files). Also `ASDPlasticMaterial3D.h` `setResponse` XML labels + `Vector(-1)` guard (marked). | (pre-campaign; PR untracked — reconstruct via `git log -- SRC/material/nD/ASDPlasticMaterial3D`) |
