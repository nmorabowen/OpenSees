---
title: OpenSees Windows compilation journal
project: Ladruno
date_started: 2026-04-30
status: working
targets:
  - OpenSees.exe
  - OpenSeesSP.exe
  - OpenSeesMP.exe
  - opensees.pyd (OpenSeesPy)
features:
  - HDF5 1.14
  - MPCO recorder
  - MUMPS 5.5.1 (LP64)
  - Intel MPI
tags:
  - opensees
  - build
  - windows
  - intel-oneapi
  - mumps
  - hdf5
  - mpco
---

# OpenSees Windows compilation journal

This file is the canonical record of how this clone of OpenSees was made to compile on Windows with all four front ends, HDF5/MPCO output, and the MUMPS sparse direct solver. Everything here is **descriptive** — the actual steps live in [setup_env.bat](../setup_env.bat) and [build.bat](../build.bat). Read this when something breaks and you need to remember *why* a patch is there.

## Goal

One repeatable script (`build.bat`) that produces from a fresh clone of upstream OpenSees:

| Target | What it is | Use case |
| --- | --- | --- |
| `OpenSees.exe` | Sequential Tcl interpreter | Single-process Tcl scripts |
| `OpenSeesSP.exe` | Partitioned-domain MPI Tcl interpreter | One big problem, MPI-distributed |
| `OpenSeesMP.exe` | Multi-interpreter MPI | Parameter sweeps, embarrassingly parallel |
| `opensees.pyd` | Python module | OpenSeesPy with HDF5 + MPCO recorder |

## Toolchain

| Tool | Version | Path | Role |
| --- | --- | --- | --- |
| Visual Studio | 2022 Community | `C:\Program Files\Microsoft Visual Studio\2022\Community\` | `cl.exe` (C/C++), `link.exe`, Windows SDK |
| Intel oneAPI | 2025.1 | `C:\Program Files (x86)\Intel\oneAPI\` | `ifx.exe` (Fortran), MKL (BLAS/LAPACK/ScaLAPACK), Intel MPI |
| CMake | 4.3.0 | `C:\Program Files\CMake\bin\` | Build generator |
| Ninja | latest (winget) | on PATH | Fast build executor |
| Python | 3.11.9 | `C:\Users\nmora\AppData\Local\Programs\Python\Python311\` | OpenSeesPy ABI target + Conan host |
| Conan | 2.28.1 | pip --user | Fetches HDF5, TCL, ZLIB, Eigen3 |
| MUMPS | 5.5.1 | local build → `mumps-install/` | Parallel sparse direct solver |

## Reproducible flow

All build scripts live in [../Ladruno_scripts/](../Ladruno_scripts/). All build outputs land in [../Ladruno_files/](../Ladruno_files/) and [../dist/](../dist/).

```
cd C:\Users\nmora\Github\OpenSees_Compile
call Ladruno_scripts\setup_env.bat   :: per-shell: load MSVC + Intel oneAPI + Conan
Ladruno_scripts\build.bat            :: full pipeline (Conan, MUMPS once, CMake, Ninja, dist)
```

`Ladruno_scripts\build.bat <target>` rebuilds a single target. `Ladruno_scripts\build.bat clean` wipes everything except `mumps-install/`.

After `build.bat` succeeds, all four artifacts plus their runtime DLLs and Tcl init are in `dist/bin/` and `dist/lib/tcl8.6/`. `dist/` is self-contained — no need to source `setup_env.bat` to run `OpenSees.exe` or `import opensees`. The MPI variants still need `setup_env.bat` because they're launched through `mpiexec`.

## Distributable installer

To package `dist/` into something shareable:

```
powershell -ExecutionPolicy Bypass -File Ladruno_scripts\make_installer.ps1
```

Outputs to `Ladruno_files/`:

- `Ladruno_OpenSees_<YYYYMMDD>.zip` — the binaries (zip of `dist/`)
- `install.ps1` — self-contained end-user installer with the LADRUNO ASCII banner embedded. Run with `powershell -ExecutionPolicy Bypass -File install.ps1`. Prompts for install location (default `%LOCALAPPDATA%\Ladruno\OpenSees`), expands the zip, optionally adds `bin/` to user PATH.

Distribute both files together. Receivers need Python 3.11 for OpenSeesPy and Intel oneAPI for the MPI variants — those aren't bundled.

## Python virtualenv (developer workflow)

```
powershell -ExecutionPolicy Bypass -File Ladruno_scripts\wire_pyenv.ps1
```

Creates `opensees_venv/` at the project root if missing, then writes a `.pth` file at `<venv>/Lib/site-packages/ladruno_opensees.pth` containing one line — the absolute path to `dist/bin/`. Python's `site` module reads the .pth at interpreter startup and adds that line to `sys.path`, so `import opensees` finds `dist/bin/opensees.pyd` automatically. The MKL/iomp DLLs co-located with the .pyd are picked up by Windows' DLL loader because they're in the same directory.

This is **not** a wheel install. It's the developer convenience pattern:

- Edit OpenSees source
- `Ladruno_scripts\build.bat OpenSeesPy`
- Re-run venv Python — new behavior without any re-install step

Pass `-Force` to wipe and recreate the venv. The script also runs a smoke test (one-element truss → `ux=0.01`) and writes a record to `Ladruno_files/wire_pyenv.log`.

> [!warning]
> The Python interpreter the venv inherits MUST match the ABI `opensees.pyd` was built for — currently CPython 3.11. The script warns if the selected interpreter is a different minor version.

## Banner workflow

Splash banner lives in [`Ladruno_scripts/banner_ASCII.txt`](../Ladruno_scripts/banner_ASCII.txt) (UTF-8). To change it:

1. Edit `banner_ASCII.txt`
2. `python Ladruno_scripts\patch_banner.py`
3. Rebuild: `Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP`

The script splices the new banner between `// BANNER-START` and `// BANNER-END` markers in [`OpenSees/SRC/tcl/tclMain.cpp`](../OpenSees/SRC/tcl/tclMain.cpp:114), backing up the previous cpp to `Ladruno_files/patch_banner_backups/` first.

> [!warning]
> The script uses **plain string slicing**, not `re.sub`. Python's `re.sub` interprets `\n` in the replacement string as a real newline character (because it processes regex backreference escapes), which produces broken multi-line C string literals. The bug is silent — the file looks fine in a normal editor — but the C++ compiler errors out on `multi-line string literal`.

## Source patches we apply to upstream OpenSees

These are local edits to `OpenSees/CMakeLists.txt` and a few subdirectory files. If you re-clone upstream, re-apply them. Each patch has explanatory inline comments at the patch site.

### 1. Dangling `mumps` target reference

**Where:** `CMakeLists.txt:805`

**Bug:** Original code calls `add_dependencies(OpenSeesMP mumps)` when `MUMPS_DIR` is undefined, but the `mumps` target only exists if `ExternalProject_Add(mumps ...)` runs — which is commented out in the upstream file.

**Fix:** Guard with `if (TARGET mumps)`.

### 2. MPI include paths must be 8.3 short form on Windows

**Where:** `CMakeLists.txt:405` (after the second `find_package(MPI)`)

**Bug:** CMake's Ninja generator emits broken response-file quoting for paths containing both spaces *and* parentheses when used as Fortran include flags. The MPI include path `C:\Program Files (x86)\Intel\oneAPI\mpi\2021.15\include` ends up rendered as `-IC:\Program" Files "(x86)\...` — quote in the wrong place — and ifx parses `Files` and `(x86)\...` as separate "files" to compile.

**Fix:** A small loop converts every entry of `MPI_*_INCLUDE_DIRS` / `MPI_*_INCLUDE_PATH` to its 8.3 short name (`C:\PROGRA~2\Intel\...`) on Windows. No spaces, no parens, no quoting bug.

### 3. Recursive include glob must be C/CXX-only

**Where:** `CMakeLists.txt:350`

**Bug:** Upstream uses `file(GLOB_RECURSE ... *.h)` then adds every directory containing a header to global `include_directories()`. On this codebase that produces ~200 directories. The Fortran preprocess rule wraps everything in `cmd.exe /C "..."`, which has a hard ~8 KB command-line limit. With 200 `-I` flags the line balloons to ~18 KB and ifx fails with `The command line is too long.`

**Fix:** Wrap each glob-derived include in `$<$<COMPILE_LANGUAGE:C,CXX>:...>` so Fortran targets don't see them. Drops the Fortran INCLUDES line from ~18 KB to ~4 KB. Fortran sources don't reference any C++ headers, so this is harmless.

### 4. `DistributedSuperLU.cpp` `stat` collides with MSVC's POSIX `stat()`

**Where:** `OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.cpp:69`

**Bug:** MSVC's `<sys/stat.h>` (transitively pulled in by SuperLU headers) declares `stat()` as a function. The OpenSees source declares a global variable also named `stat` of type `SuperLUStat_t`. C++ overload resolution prefers the function over the variable, so passing `&stat` to MUMPS routines yields a function pointer instead of `SuperLUStat_t*`.

**Fix:** Rename the global `stat` → `superlu_stat` and update both `&stat` references inside the file.

### 5. `MUMPS_INCLUDE_DIR` variable for clean install layouts

**Where:** `CMakeLists.txt:741, 801` (the SP and MP `target_include_directories` calls)

**Bug:** Upstream hardcodes `${MUMPS_DIR}/_deps/mumps-src/include` for MUMPS headers, which is the `ExternalProject_Add` in-tree layout. When MUMPS is installed via `cmake --install` to a separate prefix, headers live at `<prefix>/include/`.

**Fix:** Add a `MUMPS_INCLUDE_DIR` variable; if set, use it; else fall back to the in-tree path. `build.bat` sets it to `mumps-install/include`.

### 6. MUMPS library names need `.lib` suffix

**Where:** `CMakeLists.txt:451`

**Bug:** Upstream sets `MUMPS_LIBRARIES "${MUMPS_DIR}/dmumps;..."` without the `.lib` suffix. CMake's auto-suffix only kicks in for plain library *names* (no path component), not for full paths. Ninja then tries to find a *target* named exactly `dmumps` and fails with `missing and no known rule to make it`.

**Fix:** Append `.lib` explicitly.

### 7. ScaLAPACK ILP64 → LP64

**Where:** `CMakeLists.txt:512`

**Bug:** Upstream picks `mkl_scalapack_ilp64` / `mkl_intel_ilp64` / `mkl_blacs_intelmpi_ilp64` — the 64-bit-integer MKL variants. They only work if MUMPS was built with `-Dintsize64=ON` *and* Intel MPI is in ILP64 mode. We build MUMPS LP64 (default), so we link the LP64 ScaLAPACK variants. Mismatch surfaces as MUMPS "Instance Error 1 in DMUMPS_F77" with garbage `INSTANCE_NUMBER`.

**Fix:** Replace every `ilp64` with `lp64` in this line. Documented inline at the patch site so anyone flipping MUMPS back to ILP64 knows what to flip here.

### 8. Per-target compilation of `myCommands.cpp` and `Tcl_generateInterfacePoints.cpp`

**Where:** `CMakeLists.txt:598, 705, 736, 800` and removals from `SRC/modelbuilder/tcl/CMakeLists.txt` and `SRC/tcl/CMakeLists.txt`

**Bug:** Both files declare `extern PartitionedDomain theDomain;` if `_PARALLEL_PROCESSING` is defined, otherwise `extern Domain theDomain;`. The matching definition lives in `commands.cpp`, which is compiled per-exe with the right define. Upstream puts both files inside the static library `OPS_InterpTcl`, which is built **once** with no parallel define. SP/MP exes then link an `OPS_InterpTcl` whose `myCommands.cpp` expects `Domain theDomain` while their own `commands.cpp` defines `PartitionedDomain theDomain`. Different mangled names → unresolved external.

**Considered alternative:** Convert `OPS_InterpTcl` to an INTERFACE library so each consumer compiles all of its sources with its own defines. That broke `OpenSeesPy` because `elementAPI_TCL.cpp` (also in `OPS_InterpTcl`) defines API functions that conflict with `OpenSeesCommands.cpp` in `OpenSeesLIB.lib` — duplicate-symbol link errors (`LNK2005`).

**Fix:** Pull only the two parallel-sensitive files out of `OPS_InterpTcl` and add them to each consuming exe via a top-level CMake variable `OPS_INTERP_TCL_PER_TARGET_SOURCES`. The rest of `OPS_InterpTcl` stays a single static library shared by all three exes. `OpenSeesPy` doesn't need them.

## MUMPS build

`build.bat` builds MUMPS once into `mumps-install/` if `mumps-install/lib/dmumps.lib` is missing.

- Source archive: `mumps-archive/mumps_src.tar.gz` (auto-downloaded from `https://mumps-solver.org/MUMPS_5.5.1.tar.gz` if absent — the original `mumps.enseeiht.fr` URL embedded in scivision/mumps's `libraries.json` is dead).
- Build system: `mumps-src/` is a clone of `https://github.com/scivision/mumps.git` at tag `v5.5.1.5`. It uses `FetchContent_Populate` with the `local` option to consume the local archive instead of fetching.
- Configure: Ninja, Intel ifx + cl, `-Darith=d`, `BUILD_SHARED_LIBS=OFF`, `CMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded`, install prefix `mumps-install/`.
- **Critical: `-Dintsize64=OFF` passed explicitly.** See below.
- Output: `mumps-install/lib/{dmumps,mumps_common,pord}.lib` and `mumps-install/include/dmumps_c.h` etc.

### The `intsize64` cache-stickiness gotcha

scivision/mumps's `option(intsize64 ...)` defaults to OFF, but the value is sticky in CMake's cache **and** in the FetchContent-pulled subbuild's own cache. If a previous configure ran with `-Dintsize64=ON` and you `rm -rf mumps-build/` and reconfigure without the flag, CMake's FetchContent layer retains the old subbuild cache and the build still applies `-i8` to Fortran.

Symptom: `mumps_int_def.h` reports `MUMPS_INTSIZE64`, the install's `dmumps.lib` was compiled with 8-byte INTEGER, but `MUMPS_INT` in `mumps_c_types.h` is `int` (4 bytes) because the header generator believed the new flag. Result: ABI split inside the same `DMUMPS_STRUC_C` struct — C writes 4 bytes for `id.job`, Fortran reads 8 bytes for `INSTANCE_NUMBER`. Surfaces at runtime as **"Instance Error 1 in DMUMPS_F77"** with `INSTANCE_NUMBER` printed as a wildly large garbage integer (e.g. `4645744490609377280`).

**Fix:** `build.bat` always passes `-Dintsize64=OFF` explicitly to defeat this stickiness. Don't remove that flag.

### Why LP64 over ILP64

LP64 (32-bit Fortran INTEGER) matches Intel MPI's default LP64 layer (32-bit MPI handles) and MKL's `_lp64` ScaLAPACK / BLACS libraries. Everything is 32-bit-int-aligned and there is no ABI split.

ILP64 is for problems with more than 2^31 nonzeros. To use it, you'd flip `intsize64` to ON, switch every `lp64` to `ilp64` in [`CMakeLists.txt:512`](../OpenSees/CMakeLists.txt:512), and ensure Intel MPI is in ILP64 mode (`I_MPI_F77`/`I_MPI_F90` env or the ILP64 wrapper DLLs). Don't switch unless you actually need it.

## OpenSeesPy specifics

- Built as `OpenSeesPy.dll` by CMake. `build.bat`'s dist step renames it to `opensees.pyd` because Python expects the `.pyd` extension and the import name `opensees` (per the upstream install rule that renames to `opensees.so` on Linux).
- Python 3.11 is forced via `-DPython_EXECUTABLE=...`. Without that, CMake's `find_package(Python)` picks up Python 3.14 from `C:\Users\nmora\AppData\Local\Python\pythoncore-3.14-64\` and you get an `.pyd` that won't load in 3.11. (CMake's second `find_package(Python3)` call later still picks up 3.14 — that's harmless because it's used by an unrelated build helper, not by `OpenSeesPy` itself.)
- For self-contained imports without sourcing `setup_env.bat`, `dist/bin/` ships these Intel runtime DLLs alongside `opensees.pyd`:
  - `mkl_intel_thread.2.dll`, `mkl_core.2.dll`, `mkl_def.2.dll` (kernel runtime — MKL `dlopen`s `mkl_def.2.dll` lazily; without it you get "INTEL oneMKL FATAL ERROR: Cannot load mkl_def.2.dll")
  - `mkl_avx2.2.dll`, `mkl_avx512.2.dll`, `mkl_mc3.2.dll` (CPU-specific kernels MKL picks at runtime)
  - `mkl_scalapack_lp64.2.dll`, `mkl_blacs_intelmpi_lp64.2.dll` (for SP/MP)
  - `libiomp5md.dll` (OpenMP runtime)
- We do **not** ship `impi.dll` in `dist/` — when present alongside the SP/MP exe it gets loaded preferentially over the system Intel MPI runtime and subtly mismatches what `mpiexec` was built against. SP/MP launches always go through `mpiexec` after `setup_env.bat`, so Intel MPI's official PATH layout always wins.

## `setup_env.bat` — Windows oneAPI activation gotchas

Four Windows-specific quirks that broke our environment activation, all baked into `setup_env.bat`:

1. **`if (...)` blocks plus `(x86)` paths.** Inside a parenthesized `if` body in cmd.exe batch, parens in expanded paths prematurely close the body. Use goto-style error labels, not `if not exist X (echo err & exit /b 1)`.

2. **Intel's `setvars.bat` iteration is broken on this machine.** `call setvars.bat intel64 vs2022` runs but emits `'vars.bat' is not recognized` for every component. We sidestep by calling `compiler\latest\env\vars.bat`, `mkl\latest\env\vars.bat`, and `mpi\latest\env\vars.bat` directly.

3. **`vswhere.exe` must be on PATH.** It lives at `C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe` but isn't on PATH by default. Conan's per-package `conanvcvars.bat` and TCL 8.6.11's `makefile.vc` both need vswhere — without it, TCL build fails cryptically with `'nmakehlp' not recognized` and `rules.vc(1070): syntax error`.

4. **`NoDefaultCurrentDirectoryInExePath=1`** is set in user environment as a Windows-hardening default. cmd then refuses to auto-search `.` for executables. TCL 8.6.11's makefile builds `nmakehlp.exe` in cwd then invokes it as bare `nmakehlp` (no `.\` prefix) — fails with the variable set. We `set "NoDefaultCurrentDirectoryInExePath="` for the build session only; the user/machine setting is untouched.

## Splash banner

[`OpenSees/SRC/tcl/tclMain.cpp`](../OpenSees/SRC/tcl/tclMain.cpp) prints an ASCII-art Lardruno banner (221×21) before the standard version block. Gated on `_isatty(_fileno(stderr))` — when stderr is redirected (file/pipe/CI), the art is suppressed and only the original banner prints. On a narrow terminal it visually wraps; the user is expected to maximize.

`<windows.h>` couldn't be included in this file because Tcl's headers and `<cstdio>` fight over `vsnprintf` / `Offset` macros. Workaround: forward-declare `_isatty` and `_fileno` directly with `extern "C"`, no Windows header needed.

`OpenSeesPy` doesn't go through `tclMain.cpp`, so `import opensees` is silent (correct — the recorder/STKO banner from MPCO at run-time is enough).

## Verification

Smoke tests live at the project root.

| Test | Target | Expected |
| --- | --- | --- |
| `dist\bin\OpenSees.exe test_sp_mumps.tcl` | sequential | Tcl banner + cantilever solve via BandSPD |
| `python test_mpco.py` | `opensees.pyd` | `ux=0.01`, `.mpco` file written via HDF5, ASDEA banner |
| `mpiexec -n 1 dist\bin\OpenSeesSP.exe test_sp_mumps.tcl` | SP + MUMPS | `uy = -1.66667` (matches `PL^3/3EI`) |
| `mpiexec -n 2 dist\bin\OpenSeesSP.exe test_sp_mumps.tcl` | SP + MUMPS, 2 ranks | partitioner messages, same answer |
| `mpiexec -n 4 dist\bin\OpenSeesMP.exe test_mp.tcl` | MP, 4 ranks | each rank prints its own `ux = 10/k` for `k = 1000…1300` |

All five passed as of the build documented here.

## See also

- [setup_env.bat](../setup_env.bat) — toolchain activation
- [build.bat](../build.bat) — full pipeline driver
- [opensees-msvc-static.profile](../opensees-msvc-static.profile) — Conan profile (msvc 194, static MT runtime)
- [opensees_banner.txt](../opensees_banner.txt) — raw ASCII banner
