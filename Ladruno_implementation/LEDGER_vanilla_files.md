---
title: Ledger — vanilla files we touched
project: Ladruno
tags:
  - ledger
  - provenance
  - upstream
---

# Ledger — vanilla OpenSees files we touched

Every upstream ("vanilla") OpenSees file the Ladruno fork modifies, **why** we
touched it, and the **PR** it landed in. This is the provenance record: if we
ever rebase onto a newer upstream, this table is the list of edits to re-apply
and re-verify.

## Conventions

- **Vanilla = pre-existing upstream file.** Brand-new files we author live in
  [[LEDGER_implementations]] instead — do not list them here.
- Keep edits minimal and marked. Every Ladruno edit in a vanilla file carries a
  `// Ladruno ...` comment so `grep -rn "Ladruno" SRC/` reconstructs this table.
- One row per (file, PR). If the same file is touched by several PRs, add a row
  per PR so the history stays linear.
- PR numbers are on the fork: `github.com/nmorabowen/OpenSees` (branch `ladruno`).
- When you touch a new vanilla file, **add the row in the same PR**.

## Ledger

| Vanilla file | Why touched | PR |
|---|---|---|
| `CMakeLists.txt` (root) | Eight build patches for SP/MP/Py toolchain (see [[../Ladruno_internal/01_compilation_journal]] / `project_opensees_cmake_gotchas`); add `openseesmp` target | build + [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/classTags.h` | Register Ladruno class tags: `ELE_TAG_BezierTri6`=272, `INTEGRATOR_TAGS_ExplicitBathe`=61, `…ExplicitBatheLNVD`=63, `RECORDER_TAGS_EnergyBalanceRecorder`=26, `RECORDER_TAGS_LadrunoRecorder`=33001 | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#8](https://github.com/nmorabowen/OpenSees/pull/8) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | Broker `getNewElement/Integrator/Recorder` entries for the new class tags above | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#8](https://github.com/nmorabowen/OpenSees/pull/8) |
| `SRC/interpreter/PythonModule.cpp` | openseespy base (Patch 9 multi-module split); Ladruno splash banner + auto-generated feature list (`FEATURES-START/END`) | banner patch + [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/interpreter/PythonMPIModule.cpp` | `openseesmp` MPI module entry point (Patch 9) — re-includes `PythonModule.cpp` with `OPS_PY_MODULE_NAME` redefined | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/interpreter/OpenSeesMiscCommands.cpp` | `partition` command METIS 5 API guard / error message (Patch 9) | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | Register `recorder` keywords for EnergyBalance and Ladruno (`recorder ladruno`) | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#8](https://github.com/nmorabowen/OpenSees/pull/8) |
| `SRC/interpreter/CMakeLists.txt` | Build wiring for `openseesmp` and new interpreter sources | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/recorder/CMakeLists.txt` | Add `EnergyBalanceRecorder` and `LadrunoRecorder` + `Ladruno_*` to the recorder target | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#8](https://github.com/nmorabowen/OpenSees/pull/8) |
| `SRC/element/elasticBeamColumn/ElasticBeam3d.cpp` | `// Ladruno`: add a `"localAxes"` element response (id 30) returning the 9 packed direction cosines from `theCoordTransf->getLocalAxes`, so the Ladruno recorder can write `MODEL/LOCAL_AXES` instead of a silent identity-quaternion fallback (apeGmsh beam-orientation gap). Additive — no existing response touched; pattern to replicate on the other beams. | [#38](https://github.com/nmorabowen/OpenSees/pull/38) |
| `SRC/element/elasticBeamColumn/ElasticBeam2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) → `Vector(9)` dir cosines from `theCoordTransf->getLocalAxes` (2D transf returns full 3D cosines, vz=(0,0,1)). Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/dispBeamColumn/DispBeamColumn2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/dispBeamColumn/DispBeamColumn3d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/forceBeamColumn/ForceBeamColumn2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/forceBeamColumn/ForceBeamColumn3d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/actor/channel/CMakeLists.txt` | MP build wiring | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/actor/address/CMakeLists.txt` | MP build wiring | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/tcl/tclMain.cpp` | Ladruno splash banner (`BANNER-START/END`) + auto-generated feature list (`FEATURES-START/END`) — generated by `Ladruno_scripts/patch_banner.py` | banner patch |
| `SRC/matrix/Matrix.cpp` | `// Ladruno P4`: profiler memory hooks — `OPS_PROFILE_COUNT_ALLOC` after each per-object `data` `new[]`, `OPS_PROFILE_COUNT_FREE` before each `delete[]` (under the `fromFree==0` ownership guard), tagged `ops_profiler::ALLOC_MATRIX`. Runtime-gated on `mem()`; the static `matrixWork`/`intWork` scratch is deliberately not counted (never-freed whitelist). Additive. | [#46](https://github.com/nmorabowen/OpenSees/pull/46) |
| `SRC/matrix/Vector.cpp` | `// Ladruno P4`: same profiler memory hooks on the per-object `theData` `new[]`/`delete[]` sites (ctor/copy/dtor/setData/resize/`operator[]`/`operator=`/move-`operator=`), tagged `ops_profiler::ALLOC_VECTOR`. Runtime-gated on `mem()`. Additive. | [#46](https://github.com/nmorabowen/OpenSees/pull/46) |
| `SRC/matrix/ID.cpp` | `// Ladruno P4`: same profiler memory hooks on the per-object `data` `new[]`/`delete[]` sites (ctor×2 + `malloc` branch / copy / dtor / setData / `unique` / `operator[]`-grow / `resize`-grow / `operator=` / `insert`-grow), tagged `ops_profiler::ALLOC_ID`; `arraySize` is the allocated capacity (exact frees). Bare-`if` deletes brace-converted; `operator=` captures the old `arraySize` in a temp before its overwrite. Additive. | [#47](https://github.com/nmorabowen/OpenSees/pull/47) |
| `SRC/actor/actor/MovableObject.cpp` | `// Ladruno P4`: TaggedObject live-component census — `OPS_PROFILE_CENSUS_BORN(classTag)` in both ctor bodies, `OPS_PROFILE_CENSUS_DIED(classTag)` in the dtor (+`#include <profiler/ProfilerMacros.h>`). Seam is `MovableObject` (not `TaggedObject`): `classTag` is a plain member valid through the dtor, whereas `TaggedObject::getClassTag()` is virtual and unsafe in a base ctor/dtor. Runtime-gated on `mem()`; off-path cost is one branch on this hot path. Counts every MovableObject by raw classTag (shared integer space → superset; viewer filters by classTag band). Additive. | [#48](https://github.com/nmorabowen/OpenSees/pull/48) |
| `SRC/analysis/integrator/IncrementalIntegrator.cpp` | `// Ladruno P3`: deep per-element-type timing — `OPS_PROFILE_SCOPE_DEEP_NAMED` around the `formTangent` and `formElementResidual` FE_Element loops + `OPS_PROFILE_FE_ELEM_SCOPE` per element (folds each `getTangent`/`getResidual` wall into a per-classTag `elem_by_type` bucket). Loops re-braced to scope the per-element timer; `getElement()`/`getClassTag()` touched ONLY when the deep gate is on (`tmr.engaged()` short-circuit) so an unprofiled run pays nothing. +`#include <Element.h>`, `<profiler/ProfilerMacros.h>`. Additive/behavior-preserving. | [#49](https://github.com/nmorabowen/OpenSees/pull/49) |
| `SRC/analysis/integrator/TransientIntegrator.cpp` | `// Ladruno P3`: same deep per-element-type timing on the FE_Element tangent loop of `TransientIntegrator::formTangent` (the transient/explicit assembly path; residual path reuses `IncrementalIntegrator::formElementResidual`). The DOF_Group mass/damping loop is left untimed (not element-keyed). +`#include <Element.h>`, `<profiler/ProfilerMacros.h>`. Additive. | [#49](https://github.com/nmorabowen/OpenSees/pull/49) |

> [!note] Upstreamable bugfixes
> Some PRs fix genuine upstream bugs (not fork-only features) and are candidates
> to send back to OpenSeesFramework. Track those in the table below so we know
> what to PR upstream.

| Vanilla file | Upstreamable fix | Fork PR |
|---|---|---|
| Lagrange quad/tri elements (`SRC/element/...`) | Fix `setResponse` NdMaterialOutput tags so MPCO sees material output | [#7](https://github.com/nmorabowen/OpenSees/pull/7) |
