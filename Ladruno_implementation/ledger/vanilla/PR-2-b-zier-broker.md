---
wp: PR-2
title: "2, #6, #8, #65, #97, #108 (Bézier broker) -- 1 vanilla row(s)"
pr: "#2, #6, #8, #65, #97, #108"
files: ["`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`"]
table: "main"
legacy_seq: [14]
---
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | Broker `getNewElement/Integrator/Recorder/NDMaterial` entries for the new class tags above (incl. `ELE_TAG_LadrunoBrick` → `new LadrunoBrick()`; `ND_TAG_LogStrainNDMaterial` → `new LogStrainNDMaterial()` so parallel/database `recvSelf` can reconstruct it — #70 shipped the material but missed the broker entry; `ELE_TAG_BezierTri6`/`ELE_TAG_BezierTet10` → `new BezierTri6()`/`new BezierTet10()` — both shipped elements were MISSING from the broker, so database/MPI restore failed with "no Element type exists for class tag 33000/33001"; surfaced by the BezierTet10 corot serialization round-trip test) | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#8](https://github.com/nmorabowen/OpenSees/pull/8), [#65](https://github.com/nmorabowen/OpenSees/pull/65), [#97](https://github.com/nmorabowen/OpenSees/pull/97), [#108](https://github.com/nmorabowen/OpenSees/pull/108) (Bézier broker) |
