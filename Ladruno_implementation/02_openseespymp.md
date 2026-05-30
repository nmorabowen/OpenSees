---
title: openseesmp / openseessp — parallel Python modules
project: Ladruno
status: planned
priority: high
owner: nmora
tags:
  - implementation
  - python
  - parallel
  - solver
  - build-system
---

# openseesmp / openseessp — parallel Python modules

> [!note] Decisions locked (2026-05-16)
> - **Source-split: Option B** — three per-define static libs (`OPS_InterpPyCmds{,_MP,_SP}`). Not negotiable; it's the load-bearing fix and the explicit delta over jaabell's branch.
> - **Phasing: MP first, SP fast-follow.** Both in scope; they share ~90% of the scaffolding, so SP lands right after MP is green.
> - **Naming/coexistence: adopt jaabell's module-rename shim** — `import openseesmp` / `import openseessp` live beside `import opensees` in one venv. The old "separate dist trees" idea is dropped.
> - **Intel MPI runtime: ship it.** Bundle the pinned oneAPI-2025.1 MPI redistributable (loader DLLs + Hydra launcher) into `dist/bin` so `mpiexec`/`import openseesmp` work without an activated oneAPI shell; `build.bat` copies it.
> - **All decisions locked — plan is execution-ready.**

## What

Add two CMake targets producing **MPI-parallel Python** modules, the Python analogues of the Tcl `OpenSeesMP` / `OpenSeesSP` executables:

- **`OpenSeesPyMP`** → `import openseesmp`. `PythonModule.cpp` compiled with `_PARALLEL_INTERPRETERS`. Each MPI rank runs the full Python script with its own interpreter; the model differentiates via `ops.getPID()` / `ops.getNP()`. The *parametric / independent-interpreter* model (upstream `openseespymp`).
- **`OpenSeesPySP`** → `import openseessp`. `PythonModule.cpp` compiled with `_PARALLEL_PROCESSING`. One model partitioned across ranks; rank 0 orchestrates (upstream `openseespysp`). Narrower command subset than MP.

The module-rename shim (jaabell's, [adopted below](#prior-art-jaabellopenseespymp)) means all three import names coexist in a single directory/venv — **no separate dist trees**.

In scope: the CMake targets, the `OPS_InterpPyCmds` interpreter-source split that makes a parallel Python build *function at all* (see [§ The core blocker](#the-core-blocker)), the rename shim, `build.bat` wiring, smoke tests with analytical references.

**Not in scope** (follow-up): pip-wheel packaging mirroring upstream's three-package layout; an `mpi4py`-style high-level parallel API; an auto-selecting `import` shim.

## Why

We already build sequential `OpenSeesPy` and the Tcl `OpenSeesMP`/`OpenSeesSP`. The parallel command surface (`getPID`, `getNP`, `barrier`, `send`, `recv`, `Bcast`, the `mpi` MPI-diagonal-SOE flag, distributed SOE/numberer) **already exists** in the Python/DL command engine under `#ifdef _PARALLEL_INTERPRETERS` — see [`OpenSeesCommands.cpp:120-210, 1436-3854`](../SRC/interpreter/OpenSeesCommands.cpp). It is simply never compiled with the define in this build. Users running large soil-structure / DRM models in Python currently must drop to Tcl for domain decomposition; closing that gap is the motivation. The MUMPS LP64 / ScaLAPACK / Intel-MPI infrastructure is already proven by the Tcl MP/SP targets, so the remaining work is build-system surgery, not numerics.

## Prior art: jaabell/OpenSeesPyMP

[`jaabell/OpenSees@OpenSeesPyMP`](https://github.com/jaabell/OpenSees/tree/OpenSeesPyMP) (commits [`1de9d620e`](https://github.com/jaabell/OpenSees/commit/1de9d620e) "add OpenSeesPyMP", `997c84f9b` Conan/MUMPS). His baseline is structurally identical to ours: `OPS_INTERPRETER` is one OBJECT lib ([CMakeLists.txt:561](https://github.com/jaabell/OpenSees/blob/OpenSeesPyMP/CMakeLists.txt#L561)) with `OpenSeesCommands.cpp` compiled in once ([interp CMakeLists:5](https://github.com/jaabell/OpenSees/blob/OpenSeesPyMP/SRC/interpreter/CMakeLists.txt#L5)), folded into `OpenSeesLIB`.

**Reuse (the genuinely good parts):**

1. **Module-rename shim.** [`PythonMPIModule.cpp`](https://github.com/jaabell/OpenSees/blob/OpenSeesPyMP/SRC/interpreter/PythonMPIModule.cpp) is 4 lines: `#define OPS_PY_MODULE_NAME "openseesmp"` / `#define OPS_PY_INIT_FUNC PyInit_openseesmp` then `#include "PythonModule.cpp"`. `PythonModule.cpp` is parameterized with `#ifndef`-guarded macros defaulting to `opensees`/`PyInit_opensees` (touching `moduledef`, the init func, and the `OpenSeesError` exception name). This is minimal, upstream-mergeable, and gives coexisting import names. **We adopt this verbatim**, plus a mirror `PythonSPModule.cpp` (`openseessp`).
2. **Per-subdir `target_sources(OpenSeesPyMP …)`** in `SRC/actor/machineBroker`, `…/channel`, `…/address`, `…/sparseGEN` — idiomatic, matches how the Tcl SP/MP targets already collect their sources. **We follow this placement** instead of a top-level source list.

**The gap (what he did NOT do — and why his branch is cosmetic):** he did *not* split `OPS_INTERPRETER`, and his `PythonModule.cpp` has zero MPI code (verified). So `OpenSeesCommands::OpenSeesCommands()` — the *only* place `new MPI_MachineBroker(...)`, channel setup, `setMPIDSOEFlag`, and the `getPID`/`getNP`/`barrier` commands live — stays compiled **without** `_PARALLEL_INTERPRETERS` in the shared lib. `target_compile_definitions(... PUBLIC _PARALLEL_INTERPRETERS)` on his SHARED target only recompiles *its own* sources, not the prebuilt `OpenSeesLIB` objects. Net: `theMachineBroker` is never constructed and the parallel commands are `#ifdef`'d out — `openseesmp` imports and runs, but as effectively-sequential OpenSees that links unused MPI infra under `mpiexec` (high confidence from source structure; verifiable with `mpiexec -n 2 python -c "import openseesmp; …; print(ops.getPID())"` against his branch — predicted to fail or print 0 on every rank). He also **comments MUMPS out** (a `-fPIC` Linux caveat that doesn't map to our Windows `/MT` static build) and links a **reduced source set** (no `PartitionedDomain`/`ActorSubdomain`/`ShadowSubdomain`, no `OPS_SysOfEqnPython`, upstream `${CONAN_LIBS}`/`${HDF5_LIBRARIES}` instead of our `HDF5::HDF5` targets).

**Conclusion:** his branch = the easy 80% (rename + link infra). The `OPS_InterpPyCmds` split *is* the missing 20% that makes parallelism real. The shim and the split are **orthogonal**: the shim fixes naming/coexistence; only the split actually compiles the parallel command engine. We take his shim, supply the split, restore full MUMPS + the Tcl-parity source set.

## Where

- **Modify**: [`OpenSees/CMakeLists.txt`](../CMakeLists.txt) — `OPS_FINAL_TARGET` STRINGS list; new `OPS_InterpPyCmds{,_MP,_SP}` static libs; new `OpenSeesPyMP`/`OpenSeesPySP` targets modelled on Tcl `OpenSeesMP` ([~845-908](../CMakeLists.txt)) / `OpenSeesSP` ([~776-838](../CMakeLists.txt)); remove the two parallel-sensitive sources from `OPS_INTERPRETER`.
- **Modify**: [`OpenSees/SRC/interpreter/CMakeLists.txt`](../SRC/interpreter/CMakeLists.txt) — drop `OpenSeesCommands.cpp` + `OpenSeesMiscCommands.cpp` from `OPS_INTERPRETER`; add the three `OPS_InterpPyCmds*` source sets; `target_sources(OpenSeesPyMP/SP …)` for the Python entry trio.
- **Modify** (per-subdir `target_sources`, jaabell's pattern): `SRC/actor/machineBroker`, `SRC/actor/channel`, `SRC/actor/address`, `SRC/system_of_eqn/linearSOE/sparseGEN`, `SRC/system_of_eqn/linearSOE/mumps`, `SRC/domain/{domain/partitioned,subdomain}` CMakeLists.
- **New files**: `SRC/interpreter/PythonMPIModule.cpp` (port jaabell's verbatim), `SRC/interpreter/PythonSPModule.cpp` (mirror, `openseessp`); parameterize `SRC/interpreter/PythonModule.cpp` with `OPS_PY_MODULE_NAME`/`OPS_PY_INIT_FUNC` (port jaabell's diff).
- **Reference (precedent)**: existing **Patch 8** `OPS_INTERP_TCL_PER_TARGET_SOURCES` in [`CMakeLists.txt`](../CMakeLists.txt) — same technique, same class of reason, Tcl side.
- **Modify**: [`Ladruno_scripts/build.bat`](../Ladruno_scripts/build.bat) — `TARGETS` (~66); add `OpenSeesPyMP.dll → openseesmp.pyd` / `OpenSeesPySP.dll → openseessp.pyd` copies into the existing `dist/bin` (~156-183); the deliberately-skipped Intel-MPI-runtime copy (~183) is now in question.
- **Build journal**: record the CMake change as "Patch 9" in [[../Ladruno_internal/01_compilation_journal]] and the `OpenSees CMakeLists.txt patches` memory.

## How

### The core blocker

`OpenSeesCommands.cpp` (the Python/DL command engine, analogue of Tcl `commands.cpp`) lives in the **`OPS_INTERPRETER` OBJECT library** ([`CMakeLists.txt:589`](../CMakeLists.txt), [`SRC/interpreter/CMakeLists.txt:5`](../SRC/interpreter/CMakeLists.txt)), compiled **exactly once** with no parallel define, folded into `OpenSeesLIB` ([`CMakeLists.txt:679`](../CMakeLists.txt)) which *every* target links.

1. All `#ifdef _PARALLEL_INTERPRETERS` machinery in it — `new MPI_MachineBroker(...)`, channel setup, `setMPIDSOEFlag`, the parallel SOE/numberer commands — is **compiled out**. A define on a *new* target does **not** recompile a prebuilt OBJECT lib. "Copy OpenSeesPy + add the define" is exactly jaabell's branch, and it does nothing (see [Prior art](#prior-art-jaabellopenseespymp)).
2. You cannot hold two ABIs of one OBJECT lib in a single configure. Compiling `OpenSeesCommands.cpp` both in the shared lib (seq) and per-target (parallel) into one link → duplicate symbols (LNK2005) — the failure Patch 8's notes record for the abandoned INTERFACE-library attempt.

Enabling fix: the parallel-sensitive interpreter sources leave the single OBJECT lib and are compiled **per consuming target** with that target's own define — the Python analogue of Patch 8's `OPS_INTERP_TCL_PER_TARGET_SOURCES`.

Parallel-sensitive sources (grep `_PARALLEL_INTERPRETERS|_PARALLEL_PROCESSING` over `SRC/interpreter`): **`OpenSeesCommands.cpp`**, **`OpenSeesMiscCommands.cpp`**. Both use `OPS_GetDomain()` (returns `Domain*` polymorphically) and `OpenSeesCommands::theDomain` is always `new Domain` regardless of define ([`OpenSeesCommands.cpp:188`](../SRC/interpreter/OpenSeesCommands.cpp)) — so unlike Tcl `commands.cpp`/`PartitionedDomain theDomain`, **no conditional-global-type mangled-name landmine**. The risk is purely conditional-compilation + ODR, which the per-target split resolves cleanly.

### Decided design: Option B + jaabell's shim (orthogonal)

**Split (Option B, locked):** three static libs each compiling `OpenSeesCommands.cpp` + `OpenSeesMiscCommands.cpp` with one define:

- `OPS_InterpPyCmds` — no parallel define (seq `OpenSeesPy` links this).
- `OPS_InterpPyCmds_MP` — `_PARALLEL_INTERPRETERS ${MUMPS_FLAG}`.
- `OPS_InterpPyCmds_SP` — `_PARALLEL_PROCESSING ${MUMPS_FLAG}`.

`OPS_INTERPRETER` loses these two sources; each Python target links exactly one `OPS_InterpPyCmds*`. Tcl exes link the seq one — behaviour unchanged (they drive parallelism via `commands.cpp`, never the Python command engine's parallel paths). Mirrors how the distributed MUMPS/SOE sources are already added per-target; keeps the seq OBJECT lib clean; no per-exe source edits. *(Option A — scattering `OPS_INTERP_PY_PER_TARGET_SOURCES` into every consumer of `OPS_INTERPRETER` — rejected: larger blast radius, edits every Tcl exe.)*

**Naming (jaabell's shim, adopted):** parameterize `PythonModule.cpp`; add `PythonMPIModule.cpp` (`openseesmp`) and `PythonSPModule.cpp` (`openseessp`), each a 4-line `#define` + `#include "PythonModule.cpp"`. Compiled per-target (in `target_sources`, not `OPS_INTERPRETER`), so they pick up the target's define and emit distinct `PyInit_*`.

These are **independent**: the shim alone (jaabell) gives a renamed but sequential module; the split alone gives a working `opensees`-named parallel module that collides with seq. We need both.

### New targets

`OpenSeesPyMP` mirrors the **Tcl `OpenSeesMP`** parallel source set (not jaabell's reduced one — we want functional parity with the proven Tcl parallel build):

- Python entry (per-target, via interp CMakeLists): `PythonMPIModule.cpp`, `PythonWrapper.cpp`, `PythonStream.cpp` (do *not* also pull `PythonModule.cpp` directly — the shim `#include`s it).
- `FEM_ObjectBrokerAllClasses.cpp` (per-target, parallel define).
- Parallel infra (per-subdir `target_sources`): `MPI_MachineBroker.cpp`, `PartitionedDomain*.cpp`, `ActorSubdomain.cpp`, `ShadowSubdomain.cpp`, `MumpsParallelSOE/Solver.cpp`, `MumpsSOE/Solver.cpp`, `MPIDiagonalSOE/Solver.cpp`, plus Channel/address sources jaabell added.
- Link `OPS_InterpPyCmds_MP`.
- **Drop** Tcl-only: `mpiParameterMain.cpp`, `tclMain.cpp`, `commands.cpp`, `OPS_INTERP_TCL_PER_TARGET_SOURCES`, `OPS_InterpTcl`.

`OpenSeesPySP` — same shape modelled on Tcl `OpenSeesSP`: adds `DistributedSuperLU.cpp`, `DistributedSparseGen{Col,Row}LinSOE.cpp`; `PythonSPModule.cpp`; links `OPS_InterpPyCmds_SP`; define `_PARALLEL_PROCESSING`.

```cmake
set_target_properties(OpenSeesPyMP PROPERTIES
    PREFIX "" OUTPUT_NAME "openseesmp" EXCLUDE_FROM_ALL ON)
target_include_directories(OpenSeesPyMP PUBLIC ${Python_INCLUDE_DIRS})
target_include_directories(OpenSeesPyMP PRIVATE
    $<IF:$<BOOL:${MUMPS_INCLUDE_DIR}>,${MUMPS_INCLUDE_DIR},${MUMPS_DIR}/_deps/mumps-src/include>
    ${MPI_CXX_INCLUDE_DIRS})
target_compile_definitions(OpenSeesPyMP PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG}
    $<$<BOOL:${OPENMPI}>:_OPENMPI>)
target_link_libraries(OpenSeesPyMP
    OPS_InterpPyCmds_MP OpenSeesLIB OPS_SysOfEqnPython OPS_Reliability OPS_Recorder
    OPS_Numerics METIS SUPERLU_DIST ${MUMPS_LIBRARIES}
    HDF5::HDF5 ZLIB::ZLIB Eigen3::Eigen ${Python_LIBRARIES}
    ${MPI_CXX_LIBRARIES} ${SCALAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
```

(Link set = seq `OpenSeesPy`'s libs ∪ Tcl `OpenSeesMP`'s parallel stack, swapping `OPS_InterpTcl`→`OPS_InterpPyCmds_MP`. We keep our fork's `HDF5::HDF5`/`Eigen3::Eigen` Conan-target convention, not jaabell's `${CONAN_LIBS}`. `OpenSeesPySP` identical but `_PARALLEL_PROCESSING`, `openseessp`, SP source/lib set.)

### build.bat + packaging

- Add `OpenSeesPyMP OpenSeesPySP` to `TARGETS` (~66). They depend on prebuilt MUMPS, which `build.bat` already ensures.
- Outputs `OpenSeesPyMP.dll` / `OpenSeesPySP.dll` (CMake `OUTPUT_NAME` makes the stem `openseesmp`/`openseessp`). Copy → `openseesmp.pyd` / `openseessp.pyd` into the **existing `dist/bin`** alongside `opensees.pyd`. No separate trees, no venv collision — the rename shim is what buys this.
- **Intel MPI runtime — ship it (decided).** Replace the deliberate skip at build.bat ~183 with a copy of the pinned oneAPI-2025.1 MPI redistributable into `dist/bin`. Two parts are needed because users will `mpiexec -n N python driver.py` from a plain venv with no oneAPI env:
  1. **Loader DLLs** (resolved at `import openseesmp`): `impi.dll` (and `impimt.dll` if the MT layer is used), `libfabric.dll` + its OFI provider DLLs. Source: `%ONEAPI_ROOT%\mpi\latest\bin\` and `…\mpi\latest\opt\mpi\libfabric\bin\`.
  2. **Hydra launcher** (to spawn ranks): `mpiexec.exe`, `hydra_service.exe`, `hydra_pmi_proxy.exe`, `hydra_bstrap_proxy.exe` from `…\mpi\latest\bin\`.
  Pin the exact file set against the oneAPI **2025.1** we build with (versions drift between oneAPI releases); enumerate via `dumpbin /dependents openseesmp.pyd` + a clean-machine (no-oneAPI) smoke run, not by guesswork. Confirm the Intel oneAPI redistribution terms cover bundling these (Intel MPI runtime is redistributable, but cite the license file in the dist). Adds ~tens of MB to `dist/bin` — acceptable for a self-contained Ladruno drop.

### Runtime model (document for users)

- **MP** (`import openseesmp`): `mpiexec -n N python driver.py`. Every rank runs `driver.py` independently; branch on `ops.getPID()`/`ops.getNP()`. Parametric sweeps / independent models.
- **SP** (`import openseessp`): `mpiexec -n N python driver.py`. One decomposed model; rank 0 orchestrates. Narrower supported command subset than MP — state this explicitly in dist docs.

### Testing

Port the proven Tcl smoke tests to Python, same analytical references:

- `test_sp_mumps.py` — 10-element cantilever, tip load, `mpiexec -n {1,2} python …` with `import openseessp`; expect `uy = -1.66667` (= `PL³/3EI`).
- `test_mp.py` — 4 independent parametric trusses, `mpiexec -n 4` with `import openseesmp`; each rank prints `ux = 10/k`. **`ops.getPID()` returning the true rank (not 0 everywhere) is the regression that catches the jaabell-style cosmetic failure.**
- Seq regression: confirm the `OPS_INTERPRETER`→`OPS_InterpPyCmds` split left sequential `import opensees` byte-identical on an existing model.

## Risks / open questions

> [!note] Resolved — no open questions remain
> **Intel MPI runtime: ship it.** Decided 2026-05-16. `build.bat` bundles the pinned oneAPI-2025.1 MPI redistributable (loader DLLs + Hydra launcher) into `dist/bin` — see [build.bat + packaging](#buildbat--packaging) for the file set and the pin/verify procedure. Residual *tasks* (not decisions): (a) enumerate the exact DLL/exe set on the build machine and pin it; (b) verify on a clean no-oneAPI box; (c) cite the Intel redistribution license in the dist. None block the CMake work.

- **ODR / LNK2005**: each Python target must link *exactly one* `OPS_InterpPyCmds*`, and `OPS_INTERPRETER` must no longer contain the two sources. This is the whole point of the split — verify with `dumpbin /symbols` if a duplicate-symbol link error appears.
- **jaabell-parity regression**: the failure mode of his branch is silent (imports fine, runs sequentially). The MP smoke test's `getPID()`-per-rank assertion is the guard; do not mark MP "done" without it.
- **`setMPIDSOEFlag` global**: defined in `OpenSeesCommands.cpp` under `_PARALLEL_INTERPRETERS` ([line 121](../SRC/interpreter/OpenSeesCommands.cpp)). Per-lib split → each `.pyd` gets its own copy (separate DLLs, fine). Ensure no parallel `OPS_InterpPyCmds*` leaks into a Tcl exe link.
- **MUMPS as static `.lib` in a SHARED `.pyd`**: jaabell deferred MUMPS citing `-fPIC`. On our Windows `/MT` static build the Tcl `OpenSeesMP.exe` already links static MUMPS successfully; a SHARED `.pyd` pulling the same `.lib` (+ MKL/MPI) should behave identically, but this is the highest-uncertainty link step — validate the MP build links and runs MUMPS before declaring success; if it fails, the SP/MP `.pyd` may need MUMPS as a DLL.
- **MPI integer size**: settled — LP64 throughout (MUMPS `intsize64=OFF`, ScaLAPACK `*_lp64`); Python modules inherit the same link, no new ABI exposure. Don't revisit unless those flip (MUMPS-build memory / Patch 7).
- **Python ABI**: `*.pyd` links `python312.dll`; `-DPython_EXECUTABLE=` must point at the 3.12 install (build-env memory).
- **Backwards compatibility**: sequential `import opensees` must be byte-for-byte unchanged after the split — most-used target, and the split touches its command engine's compilation path.

## Implementation log

**2026-05-16 — MP implemented & verified working.** Status: `OpenSeesPyMP` (`import openseesmp`) builds and runs real MPI parallelism. SP deferred (fast-follow). Recorded as "Patch 9" in the cmake-gotchas memory and [[../Ladruno_internal/01_compilation_journal]] (§9, full detail there).

Done exactly per the locked design (Option B + jaabell rename-shim, orthogonal):

- `PythonModule.cpp` parameterized with `#ifndef`-guarded `OPS_PY_MODULE_NAME`/`OPS_PY_INIT_FUNC` (compile-time string concat for the exception name — cleaner than jaabell's runtime `snprintf`); sequential output byte-identical.
- New `PythonMPIModule.cpp` 4-line shim (`openseesmp`).
- Patch 9 CMake: `OpenSeesCommands.cpp`+`OpenSeesMiscCommands.cpp` removed from `OPS_INTERPRETER`; three STATIC libs `OPS_InterpPyCmds{,_MP,_SP}`; G3/OpenSees/SP/MP/seq-Py rewired to the seq lib (kept out of `OpenSeesLIB`); new `OpenSeesPyMP` SHARED target (`OUTPUT_NAME openseesmp`, Tcl-`OpenSeesMP`-parity sources, links `OPS_InterpPyCmds_MP`).
- The `OpenSeesCommands.h` ABI-stability pivot was verified **before** writing CMake — the split is provably free of a struct-layout ODR landmine.

Four build iterations (root causes + fixes, all from latent never-compiled parallel code or build-script bugs — none from the Patch-9 design itself):

1. `OpenSeesMiscCommands.cpp` `<metis.h>` under `_PARALLEL_INTERPRETERS` → `OPS_InterpPyCmds_MP/_SP` `PRIVATE`-link `METIS SUPERLU_DIST` for their interface includes.
2. `OPS_partition()` uses the **METIS 5 API**; fork bundles **METIS 4.0.1**. Guarded behind `OPS_HAVE_METIS5` (default → honest runtime error; body preserved verbatim). Not on the MP-first critical path.
3. `LNK2019 MPI_Channel::MPI_Channel(int)` → added `MPI_Channel.cpp`/`MPI_ChannelAddress.cpp` to `OpenSeesPyMP` in `SRC/actor/{channel,address}/CMakeLists.txt` (Tcl-parity).
4. **build.bat latent bug** (it masked iteration 3): `if errorlevel 1` inside the parenthesized `for` build-loop didn't trip on a cmake/ninja link failure (false exit 0, ran Step 5 on a non-built target). Replaced with `… || (echo … & exit /b 1)`.

**Packaging deviation from the plan (deliberate, "do it best"):** `openseesmp.pyd` + its Intel MPI/MKL runtime ship in a **dedicated self-contained `dist\openseesmp\`**, NOT `dist\bin`. Reason: a `dist\bin` `impi.dll` shadows the oneAPI one and reintroduces the documented MUMPS "Instance Error 1" for the co-located Tcl `OpenSeesSP/MP.exe`. Same-oneAPI-version `mpiexec`+`impi.dll` are co-located so there is no version-skew. The rename shim still gives `import openseesmp` ≠ `import opensees`, so coexistence is preserved without sharing a directory.

**Verification — MP PASSED:** `dist\openseesmp\test_openseesmp.py` via the bundled `mpiexec -n 4` + the py3.12 the build targeted → ranks print **distinct** `[rank 0/4]…[rank 3/4]`, NP=4, correct independent `ux = P/k`, clean `Process N Terminating`. The per-rank-PID assertion is the regression guard that proves a real parallel build (vs jaabell's cosmetic one). Seq-regression (full-set rebuild + `import opensees` / `OpenSees.exe`): in progress at time of writing.

### 2026-05-17 — Installer & venv usage (both modules)

The Inno installer (`build_inno_installer.ps1` → `installer.iss`) bundles **all** of `dist/` (so `dist\bin` + `dist\openseesmp` both land under `{app}`) and wires the chosen venv to **both** via `wire_venv_pth.py <bin> <openseesmp>`, which drops `_ladruno_opensees_boot.py` + a one-line `ladruno_opensees.pth` in site-packages. The boot module, per dir: `sys.path.insert` + `os.add_dll_directory` + **prepend to the process-local `os.environ['PATH']`**, then a `PMI_RANK/PMI_SIZE`-guarded `openseespy`→`opensees` alias (skipped under `mpiexec`).

**Gotcha (load-bearing):** `os.add_dll_directory` is enough to load `openseesmp.pyd`+`impi.dll`, but at `MPI_Init` `impi.dll` does a bare `LoadLibrary("libfabric.dll")` that searches **PATH**, not the add_dll_directory set → omitting the process-local PATH prepend makes every rank abort `MPIDI_OFI_mpi_init_hook … PMPI_Init` (1090959). The PATH edit is in-memory only (this Python process + its mpiexec children); it does **not** touch user/system PATH, so it can't shadow the oneAPI `impi.dll` the Tcl `OpenSeesSP/MP.exe` rely on. `wire_pyenv.ps1` (dev fallback) now delegates to the same `wire_venv_pth.py` (one source of truth). Verified end-to-end in a throwaway venv: seq `import opensees` (no PATH set), MP `import openseesmp` via `dist\openseesmp\mpiexec.exe -n 3` (distinct ranks, clean exit), and `import openseespy.opensees`. Usage from a wired venv: normal scripts `import opensees`; MPI scripts launched as `…\openseesmp\mpiexec.exe -n N <venv>\Scripts\python.exe driver.py` with `import openseesmp`.

### 2026-05-17 — OpenSeesPySP: DEFINITIVE finding, then SP build trace REMOVED

**Got to the bottom of SP.** Root cause is architectural, not a defect: the Python/DL command engine (`OpenSeesCommands.cpp`) has **no `_PARALLEL_PROCESSING` domain-decomposition support whatsoever**. Verified by reading the source:

- `OpenSeesCommands` ctor does `theDomain = new Domain;` **unconditionally** — never a `PartitionedDomain`, even under `_PARALLEL_PROCESSING`.
- There is **zero** `PartitionedDomain` / `ShadowSubdomain` / `runActors` / `OPS_PARALLEL_PROCESSING` / partition-in-`analyze` logic anywhere in the file. Its `_PARALLEL_*` code is all `_PARALLEL_INTERPRETERS` (MP): `MPI_MachineBroker`, `getPID/getNP/barrier/send/recv`, `setMPIDSOEFlag`.
- SP is a **Tcl-only paradigm**: `OpenSeesSP.exe` gets it from `mpiMain.cpp` (master/worker: rank≠0 → `MachineBroker::runActors()`; rank 0 → `ShadowSubdomain`s + run interpreter) over a global `extern PartitionedDomain theDomain` in `commands.cpp`, with the Tcl partition machinery. **None of that exists in the Python engine.** Upstream OpenSees never built a Python SP.

⇒ A real `openseessp` is **not a fix** — it is porting that entire Tcl SP subsystem into the Python engine (scope ≥ the whole MP effort), and **low value**: parallel Python is served by the now-working **`openseesmp`**; true domain-decomposed SP by the Tcl **`OpenSeesSP.exe`**.

A prototype `OpenSeesPySP` was built to *prove* the above: it confirmed `np==1` runs only a sequential MUMPS solve (`uy=-1.66667`) while `np>1` produces no real decomposition, plus an unresolved exit-time `c0000005` (MUMPS/MPI finalize order in a path OpenSees never designed for Python).

**Decision (2026-05-17): SP build trace REMOVED.** A Python module that, at best, is "sequential MUMPS" and at worst silently mis-parallelises is not worth shipping or maintaining. So the `openseessp` build was fully backed out — *the finding above is retained as the rationale*:

- **Removed:** `OpenSeesPySP` CMake target; `OPS_InterpPyCmds_SP` lib; the `OpenSeesPySP` `target_sources` in `SRC/interpreter` + `SRC/actor/{channel,address}` CMakeLists; `OpenSeesPySP` from `build.bat` `TARGETS` + the `dist\openseessp\` packaging/echoes; `SRC/interpreter/PythonSPModule.cpp`; `dist\openseessp\`. `OpenSeesCommands.cpp` was `git checkout`-reverted (every edit there had been SP-only — MP compiled it unmodified, so this is byte-clean for MP/seq).
- **Retained:** all MP work (`OpenSeesPyMP`/`openseesmp`, `OPS_InterpPyCmds_MP`, `PythonMPIModule.cpp`, the `OPS_PY_MODULE_NAME` parameterization of `PythonModule.cpp`, the `OPS_partition` METIS-4-vs-5 stub in `OpenSeesMiscCommands.cpp` — all MP-required) and this definitive write-up.
- **Verified clean** post-removal: CMake reconfigures with no dangling `OpenSeesPySP`; `openseesmp -n 4` → distinct rank 0–3, correct `ux=P/k`; seq `import opensees` → `ux=5.0 OK`.

> [!todo]
> **If a real Python-SP need ever arises:** it is a from-scratch subsystem port (mpiMain master/worker + `PartitionedDomain`-backed `theDomain` + `ShadowSubdomain` auto-creation + partition-in-`analyze`), not a revert. This section is the complete root-cause + scoping reference to start from. Default recommendation: use `openseesmp` for parallel Python and the Tcl `OpenSeesSP.exe` for true SP instead.

> [!todo]
> **`OPS_HAVE_METIS5`** — add a METIS 5 (Conan `metis/5.2.1`) for the parallel Python libs to restore the `partition` command and enable domain-partitioned MP / full SP. `idx_t` width must match LP64 (`IDXTYPEWIDTH 32`). Genuine dependency decision — weigh before SP work that needs partitioning.

> [!todo]
> Intel MPI redist residuals (from the original plan): pin the exact `dist\openseesmp` DLL/exe set to oneAPI 2025.1; clean no-oneAPI-shell box test; cite the Intel redistribution license in the dist.
