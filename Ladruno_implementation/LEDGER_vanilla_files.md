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
| `SRC/classTags.h` | Register Ladruno class tags (private band ≥33000): `ELE_TAG_BezierTri6`=33000, `ELE_TAG_BezierTet10`=33001, `ELE_TAG_LadrunoBrick`=33002, `INTEGRATOR_TAGS_ExplicitBathe`=33000…`CentralDifferenceLadruno`=33003, `RECORDER_TAGS_EnergyBalanceRecorder`=33000, `RECORDER_TAGS_LadrunoRecorder`=33001 | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#8](https://github.com/nmorabowen/OpenSees/pull/8), [#21](https://github.com/nmorabowen/OpenSees/pull/21), [#25](https://github.com/nmorabowen/OpenSees/pull/25), [#65](https://github.com/nmorabowen/OpenSees/pull/65), _pending (BezierTet10)_ |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | Broker `getNewElement/Integrator/Recorder/NDMaterial` entries for the new class tags above (incl. `ELE_TAG_LadrunoBrick` → `new LadrunoBrick()`; `ND_TAG_LogStrainNDMaterial` → `new LogStrainNDMaterial()` so parallel/database `recvSelf` can reconstruct it — #70 shipped the material but missed the broker entry; `ELE_TAG_BezierTri6`/`ELE_TAG_BezierTet10` → `new BezierTri6()`/`new BezierTet10()` — both shipped elements were MISSING from the broker, so database/MPI restore failed with "no Element type exists for class tag 33000/33001"; surfaced by the BezierTet10 corot serialization round-trip test) | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#8](https://github.com/nmorabowen/OpenSees/pull/8), [#65](https://github.com/nmorabowen/OpenSees/pull/65), _pending (LogStrain broker)_, [#108](https://github.com/nmorabowen/OpenSees/pull/108) (Bézier broker) |
| `SRC/interpreter/PythonModule.cpp` | openseespy base (Patch 9 multi-module split); Ladruno splash banner + auto-generated feature list (`FEATURES-START/END`) | banner patch + [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/interpreter/PythonMPIModule.cpp` | `openseesmp` MPI module entry point (Patch 9) — re-includes `PythonModule.cpp` with `OPS_PY_MODULE_NAME` redefined | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/interpreter/OpenSeesMiscCommands.cpp` | `partition` command METIS 5 API guard / error message (Patch 9) | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | Register `recorder` keywords for EnergyBalance and Ladruno (`recorder ladruno`) — the shared map used by OpenSeesPy/openseesmp **and** the interpreter-based Tcl (`TclWrapper`) | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#8](https://github.com/nmorabowen/OpenSees/pull/8) |
| `SRC/recorder/TclRecorderCommands.cpp` | `// Ladruno`: classic-Tcl wiring — add the `recorder ladruno` branch (extern `OPS_LadrunoRecorder` + a `strcmp(argv[1],"ladruno")` block mirroring the `mpco` one), so `recorder ladruno` also works in the classic `OpenSees.exe` Tcl shell (`commands.cpp`→`TclAddRecorder`), not only the Python/interpreter path. Additive. | [#62](https://github.com/nmorabowen/OpenSees/pull/62) |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | Register `element` dispatch for `BezierTri6`, `BezierTet10`, and `LadrunoBrick`/`ladrunoBrick` (fwd-decl + `functionMap`) | [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#65](https://github.com/nmorabowen/OpenSees/pull/65), _pending (BezierTet10)_ |
| `SRC/element/CMakeLists.txt` | `add_subdirectory(bezierTriangle)` + `add_subdirectory(bezierTetrahedron)` + `add_subdirectory(solidTransformation)` + `add_subdirectory(ladrunoBrick)` | [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#65](https://github.com/nmorabowen/OpenSees/pull/65), _pending (BezierTet10, SolidTransformation)_ |
| `SRC/interpreter/CMakeLists.txt` | Build wiring for `openseesmp` and new interpreter sources | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/recorder/CMakeLists.txt` | Add `EnergyBalanceRecorder` and `LadrunoRecorder` + `Ladruno_*` to the recorder target | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#8](https://github.com/nmorabowen/OpenSees/pull/8) |
| `SRC/classTags.h` | Register `RECORDER_TAGS_LadrunoMonitorRecorder`=33002 (live analysis-monitor recorder) | [#64](https://github.com/nmorabowen/OpenSees/pull/64) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | Register `recorder Monitor` keyword → `OPS_LadrunoMonitorRecorder` (HDF5 block) | [#64](https://github.com/nmorabowen/OpenSees/pull/64) |
| `SRC/recorder/CMakeLists.txt` | Add `LadrunoMonitorRecorder` + `LadrunoMonitorSink` to the recorder target (HDF5≥1.12 block, SWMR) | [#64](https://github.com/nmorabowen/OpenSees/pull/64) |
| `SRC/element/elasticBeamColumn/ElasticBeam3d.cpp` | `// Ladruno`: add a `"localAxes"` element response (id 30) returning the 9 packed direction cosines from `theCoordTransf->getLocalAxes`, so the Ladruno recorder can write `MODEL/LOCAL_AXES` instead of a silent identity-quaternion fallback (apeGmsh beam-orientation gap). Additive — no existing response touched; pattern to replicate on the other beams. | [#38](https://github.com/nmorabowen/OpenSees/pull/38) |
| `SRC/element/elasticBeamColumn/ElasticBeam2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) → `Vector(9)` dir cosines from `theCoordTransf->getLocalAxes` (2D transf returns full 3D cosines, vz=(0,0,1)). Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/dispBeamColumn/DispBeamColumn2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/dispBeamColumn/DispBeamColumn3d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/forceBeamColumn/ForceBeamColumn2d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/element/forceBeamColumn/ForceBeamColumn3d.cpp` | `// Ladruno`: same `"localAxes"` response (id 30) from `crdTransf->getLocalAxes`. Additive. | [#39](https://github.com/nmorabowen/OpenSees/pull/39) |
| `SRC/actor/channel/CMakeLists.txt` | MP build wiring | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/actor/address/CMakeLists.txt` | MP build wiring | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/tcl/tclMain.cpp` | Ladruno splash banner (`BANNER-START/END`) + auto-generated feature list (`FEATURES-START/END`) — generated by `Ladruno_scripts/patch_banner.py` | banner patch |
| `SRC/classTags.h` | Register `ND_TAG_LadrunoJ2`=33011 (combined iso + Chaboche AF kinematic J2) | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoJ2.h"` + `case ND_TAG_LadrunoJ2: return new LadrunoJ2();` so parallel/database `recvSelf` can reconstruct it | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoJ2()` + `nDMaterialsMap["LadrunoJ2"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/material/nD/CMakeLists.txt` | Add `LadrunoJ2.cpp` to `OPS_Material` sources + `LadrunoJ2.h` to headers | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/interpreter/PythonModule.cpp` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — add LadrunoJ2 line | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/classTags.h` | Register `ND_TAG_LadrunoJ2Finite`=33012 (finite-strain-native combined-hardening J2 with co-rotating backstress) | finite-native-j2 |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoJ2Finite.h"` + `case ND_TAG_LadrunoJ2Finite: return new LadrunoJ2Finite();` so parallel/database `recvSelf` can reconstruct it | finite-native-j2 |
| `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoJ2Finite()` + `nDMaterialsMap["LadrunoJ2Finite"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | finite-native-j2 |
| `SRC/material/nD/CMakeLists.txt` | Add `LadrunoJ2Finite.cpp` to `OPS_Material` sources + `LadrunoJ2Finite.h` to headers | finite-native-j2 |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — add LadrunoJ2Finite line | finite-native-j2 |
| `SRC/matrix/Matrix.cpp` | `// Ladruno P4`: profiler memory hooks — `OPS_PROFILE_COUNT_ALLOC` after each per-object `data` `new[]`, `OPS_PROFILE_COUNT_FREE` before each `delete[]` (under the `fromFree==0` ownership guard), tagged `ops_profiler::ALLOC_MATRIX`. Runtime-gated on `mem()`; the static `matrixWork`/`intWork` scratch is deliberately not counted (never-freed whitelist). Additive. | [#46](https://github.com/nmorabowen/OpenSees/pull/46) |
| `SRC/matrix/Vector.cpp` | `// Ladruno P4`: same profiler memory hooks on the per-object `theData` `new[]`/`delete[]` sites (ctor/copy/dtor/setData/resize/`operator[]`/`operator=`/move-`operator=`), tagged `ops_profiler::ALLOC_VECTOR`. Runtime-gated on `mem()`. Additive. | [#46](https://github.com/nmorabowen/OpenSees/pull/46) |
| `SRC/matrix/ID.cpp` | `// Ladruno P4`: same profiler memory hooks on the per-object `data` `new[]`/`delete[]` sites (ctor×2 + `malloc` branch / copy / dtor / setData / `unique` / `operator[]`-grow / `resize`-grow / `operator=` / `insert`-grow), tagged `ops_profiler::ALLOC_ID`; `arraySize` is the allocated capacity (exact frees). Bare-`if` deletes brace-converted; `operator=` captures the old `arraySize` in a temp before its overwrite. Additive. | [#47](https://github.com/nmorabowen/OpenSees/pull/47) |
| `SRC/actor/actor/MovableObject.cpp` | `// Ladruno P4`: TaggedObject live-component census — `OPS_PROFILE_CENSUS_BORN(classTag)` in both ctor bodies, `OPS_PROFILE_CENSUS_DIED(classTag)` in the dtor (+`#include <profiler/ProfilerMacros.h>`). Seam is `MovableObject` (not `TaggedObject`): `classTag` is a plain member valid through the dtor, whereas `TaggedObject::getClassTag()` is virtual and unsafe in a base ctor/dtor. Runtime-gated on `mem()`; off-path cost is one branch on this hot path. Counts every MovableObject by raw classTag (shared integer space → superset; viewer filters by classTag band). Additive. | [#48](https://github.com/nmorabowen/OpenSees/pull/48) |
| `SRC/analysis/integrator/IncrementalIntegrator.cpp` | `// Ladruno P3`: deep per-element-type timing — `OPS_PROFILE_SCOPE_DEEP_NAMED` around the `formTangent` and `formElementResidual` FE_Element loops + `OPS_PROFILE_FE_ELEM_SCOPE` per element (folds each `getTangent`/`getResidual` wall into a per-classTag `elem_by_type` bucket). Loops re-braced to scope the per-element timer; `getElement()`/`getClassTag()` touched ONLY when the deep gate is on (`tmr.engaged()` short-circuit) so an unprofiled run pays nothing. +`#include <Element.h>`, `<profiler/ProfilerMacros.h>`. Additive/behavior-preserving. | [#49](https://github.com/nmorabowen/OpenSees/pull/49) |
| `SRC/analysis/integrator/TransientIntegrator.cpp` | `// Ladruno P3`: same deep per-element-type timing on the FE_Element tangent loop of `TransientIntegrator::formTangent` (the transient/explicit assembly path; residual path reuses `IncrementalIntegrator::formElementResidual`). The DOF_Group mass/damping loop is left untimed (not element-keyed). +`#include <Element.h>`, `<profiler/ProfilerMacros.h>`. Additive. | [#49](https://github.com/nmorabowen/OpenSees/pull/49) |
| `SRC/analysis/analysis/DirectIntegrationAnalysis.cpp` | Profiler phase seams (P2, `OPS_PROFILE_SCOPE` step/newStep/solveCurrentStep/commit in `analyzeStep`) + `// Ladruno P0#3` per-step series hook (`OPS_PROFILE_STEP(getCommitTag, getCurrentTime, dT, getNumIterations)` after commit). +`#include <profiler/ProfilerMacros.h>`. Additive/behavior-preserving. | [#35](https://github.com/nmorabowen/OpenSees/pull/35), [#52](https://github.com/nmorabowen/OpenSees/pull/52) |
| `SRC/analysis/analysis/StaticAnalysis.cpp` | Profiler phase seams (P2, `OPS_PROFILE_SCOPE` step/newStep/solveCurrentStep/commit in `analyze`) + `// Ladruno P0#3` per-step series hook (`OPS_PROFILE_STEP`, dt=0 for static load-stepping) at the end of each loop iteration. +`#include <profiler/ProfilerMacros.h>`. Additive/behavior-preserving. | [#35](https://github.com/nmorabowen/OpenSees/pull/35), [#52](https://github.com/nmorabowen/OpenSees/pull/52) |
| `SRC/classTags.h` | Register `MAT_TAG_LadrunoUniaxialJ2`=33000 (uniaxial combined iso + Chaboche AF kinematic J2; first Ladruno *uniaxial* tag, 33000-band is per-registry) | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoUniaxialJ2.h"` + `case MAT_TAG_LadrunoUniaxialJ2: return new LadrunoUniaxialJ2();` so parallel/database `recvSelf` can reconstruct it | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoUniaxialJ2()` + `uniaxialMaterialsMap["LadrunoUniaxialJ2"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp` | `// Ladruno`: extern `OPS_LadrunoUniaxialJ2()` + `strcmp(argv[1],"LadrunoUniaxialJ2")` dispatch block (classic-Tcl `OpenSees.exe` path) | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/material/uniaxial/CMakeLists.txt` | Add `LadrunoUniaxialJ2.cpp` to `OPS_Material` sources + `LadrunoUniaxialJ2.h` to headers | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/material/nD/CMakeLists.txt` | Add header-only `LadrunoHardening.h` (shared isotropic law, consumed by both `LadrunoJ2` and uniaxial `LadrunoUniaxialJ2` — the oracle contract) to `OPS_Material` PUBLIC headers | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/classTags.h` | Register `MAT_TAG_LadrunoRebarBuckling`=33001 (rebar-buckling wrapper, Dhakal-Maekawa; sibling of LadrunoUniaxialJ2=33000 in the uniaxial band) | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoRebarBuckling.h"` + `case MAT_TAG_LadrunoRebarBuckling: return new LadrunoRebarBuckling();` (nested-material `recvSelf` reconstruction) | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoRebarBuckling()` + `uniaxialMaterialsMap["LadrunoRebarBuckling"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp` | `// Ladruno`: extern `OPS_LadrunoRebarBuckling()` + `strcmp(argv[1],"LadrunoRebarBuckling")` dispatch block (classic-Tcl `OpenSees.exe` path) | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/material/uniaxial/CMakeLists.txt` | Add `LadrunoRebarBuckling.cpp` to `OPS_Material` sources + `LadrunoRebarBuckling.h` to headers | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/element/Element.{h,cpp}` | `// Ladruno` (ADR 20 §9): add `virtual int getInterpolationWeights(const Vector& xi, Vector& N)` to the `Element` base — default returns −1 (not implemented), host elements override. Single source of truth for nodal shape-function weights so `LadrunoEmbeddedRebar` can embed a rebar node via `-host eleTag -xi …` instead of re-supplied `-shape`. Additive base-class virtual (vtable change ⇒ recompile-all, but no existing behavior touched). Overridden by fork hosts `LadrunoBrick` (trilinear) + `BezierTet10` (Bernstein). | [#175](https://github.com/nmorabowen/OpenSees/pull/175) |

> [!note] Upstreamable bugfixes
> Some PRs fix genuine upstream bugs (not fork-only features) and are candidates
> to send back to OpenSeesFramework. Track those in the table below so we know
> what to PR upstream.

| Vanilla file | Upstreamable fix | Fork PR |
|---|---|---|
| Lagrange quad/tri elements (`SRC/element/...`) | Fix `setResponse` NdMaterialOutput tags so MPCO sees material output | [#7](https://github.com/nmorabowen/OpenSees/pull/7) |
