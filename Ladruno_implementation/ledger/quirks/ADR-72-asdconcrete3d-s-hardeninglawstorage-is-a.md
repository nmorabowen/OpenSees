---
wp: ADR-72
title: "ASDConcrete3D's HardeningLawStorage is a process-global store-if-absent registry keyed by MATERIAL TAG — it survives ops.wipe() (ADR 72 P1)"
legacy_seq: 179
---
## ASDConcrete3D's HardeningLawStorage is a process-global store-if-absent registry keyed by MATERIAL TAG — it survives `ops.wipe()` (ADR 72 P1)

`ASDConcrete3DMaterial.cpp::HardeningLawStorage::store` (static singleton;
`if (item == nullptr) item = make_shared(hl)`) latches the FIRST hardening law
ever constructed for a given material tag, for the life of the process;
`recover(tag, type)` hands that original law to every later material with the
same tag. `ops.wipe()` does not clear it. Consequence in a shared pytest
process (Zone-A runs the whole battery in ONE process): a test file that
builds `ASDConcrete3D` with a common tag (e.g. 1) and *different parameters*
than another file silently poisons the latter — the ADR-72 P1 battery's
advisory test (tag 1, toy Gf) drove `test_ladrunoBrick_asdconcrete.py`'s
mesh-objectivity dissipation ratio from ≈4 to 7.03 with zero diagnostic,
alphabetical test order deciding the victim. **Rule: every test file gives its
ASDConcrete3D materials a file-unique tag** (ADR-72 uses 337218). Same class
of static-registry risk: `CrackPlanesStorage` (same file). Upstream ASDEA
code — fix-on-touch only.
