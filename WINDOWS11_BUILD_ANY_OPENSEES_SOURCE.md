# Windows 11 Guide: Compile Any OpenSees Source Tree

Last updated: 2026-03-27

## 1) Goal

Build any OpenSees fork/branch on Windows 11 using the automated 4-step workflow.

What happens:

1. Install the required Windows tools.
2. A script clones your OpenSees source, vcpkg, MUMPS, and this build-harness repo.
3. The harness files are copied into the source tree automatically.
4. The build compiles all four targets.
5. Optionally package the result as an installer or portable zip.

## 2) Automated Workflow (Recommended)

### Step 1: Install prerequisites (run once, needs admin)

```powershell
git clone https://github.com/nmorabowen/opensees_windows_compilation_files.git C:\work\harness
cd C:\work\harness
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\1_install_dependencies.ps1
```

Installs: Git, Visual Studio 2022 (C++ workload), Intel oneAPI Base + HPC, CMake, Ninja, Python 3.11. Verifies each tool after install.

Use `-DryRun` to preview without installing. Use `-SkipInnoSetup` to skip the optional Inno Setup installer.

### Step 2: Fetch source and build harness

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\2_fetch_source.ps1 `
  -OpenSeesRepo "https://github.com/OpenSees/OpenSees.git" `
  -WorkDir "C:\work"
```

This clones:

- Your OpenSees source into `C:\work\OpenSees-src`
- vcpkg into `C:\work\OpenSees-src\third_party\vcpkg`
- MUMPS into `C:\work\OpenSees-src\third_party\mumps` (pinned commit)
- This harness repo into `C:\work\opensees_windows_compilation_files`

Then copies the build harness files (CMakeLists.txt, vcpkg.json, scripts, cmake modules, source patches) into the source tree.

To build a specific fork/branch:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\2_fetch_source.ps1 `
  -OpenSeesRepo "https://github.com/YOUR_USER/OpenSees.git" `
  -OpenSeesBranch "my-feature" `
  -WorkDir "D:\builds"
```

### Step 3: Build

```powershell
cd C:\work\OpenSees-src
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\3_build.ps1
```

The script auto-initializes the VS + oneAPI environment (sets `VS2022INSTALLDIR`, calls `setvars.bat`), bootstraps vcpkg, builds MUMPS, configures and compiles OpenSees.

Common flags:

- `-SkipMumps` -- MUMPS already built from a previous run
- `-SkipTests` -- skip smoke tests after build
- `-Parallel 4` -- limit parallel jobs (default: CPU count)
- `-DryRun` -- preview without building

### Step 4: Package (optional)

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\4_package.ps1
```

Produces a portable zip and optionally a Windows installer. See `WINDOWS11_INSTALLER_EL_LADRUNO.md` for details.

## 3) Artifacts You Should Get

After step 3, in `build-win11\`:

- `OpenSees.exe` -- serial Tcl interpreter
- `OpenSeesSP.exe` -- parallel single-program (MPI)
- `OpenSeesMP.exe` -- parallel multi-program (MPI)
- `opensees.pyd` -- Python 3.11 module
- `OpenSees-Launch.cmd` / `OpenSeesSP-Launch.cmd` / `OpenSeesMP-Launch.cmd` -- launcher scripts
- `lib\tcl8.6\` -- Tcl runtime

## 4) How To Run

Serial:

```cmd
build-win11\OpenSees-Launch.cmd
```

MPI:

```cmd
build-win11\OpenSeesSP-Launch.cmd
build-win11\OpenSeesMP-Launch.cmd
```

Or manual MPI:

```cmd
cd EXAMPLES\ParallelModelMP
mpiexec -n 2 ..\..\build-win11\OpenSeesMP.exe exampleMP.tcl
```

## 5) What the Build Harness Copies

Step 2 copies these files from this repo into the OpenSees source tree:

| Source | Destination | Purpose |
|--------|------------|---------|
| `CMakeLists.txt` | root | Cross-platform CMake build system |
| `vcpkg.json` | root | Pinned dependency manifest |
| `cmake\cmake\OpenSeesDependencies.cmake` | `cmake\cmake\` | Unified dependency discovery |
| `cmake\cmake\OpenSeesDependenciesWin.cmake` | `cmake\cmake\` | Legacy Windows deps (deprecated) |
| `SCRIPTS\*.ps1`, `SCRIPTS\*.cmd` | `SCRIPTS\` | Build, install, package, MPI repair scripts |
| `SRC\...\Tcl_generateInterfacePoints.cpp` | `SRC\...\` | MSVC linker fix (ops_TheActiveDomain) |
| `SRC\...\myCommands.cpp` | `SRC\...\` | MSVC linker fix (ops_TheActiveDomain) |

The source patches fix `extern Domain theDomain` link errors on MSVC by using `ops_TheActiveDomain` from `OPS_Globals.h`.

## 6) Manual Workflow (Alternative)

If you prefer to do each step by hand instead of using the numbered scripts:

```powershell
# Clone source
mkdir C:\work && cd C:\work
git clone <YOUR_OPENSEES_REPO_URL> OpenSees-src
cd OpenSees-src
git clone https://github.com/microsoft/vcpkg third_party\vcpkg
git clone https://github.com/OpenSees/mumps.git third_party\mumps

# Clone harness and copy files
git clone https://github.com/nmorabowen/opensees_windows_compilation_files.git ..\harness
New-Item -ItemType Directory -Force -Path .\SCRIPTS, .\cmake\cmake | Out-Null
Copy-Item "..\harness\CMakeLists.txt" .\CMakeLists.txt -Force
Copy-Item "..\harness\vcpkg.json" .\vcpkg.json -Force
Copy-Item "..\harness\cmake\cmake\OpenSeesDependencies.cmake" .\cmake\cmake\ -Force
Copy-Item "..\harness\SCRIPTS\*" .\SCRIPTS\ -Force
Copy-Item "..\harness\SRC\element\UWelements\Tcl_generateInterfacePoints.cpp" .\SRC\element\UWelements\ -Force
Copy-Item "..\harness\SRC\modelbuilder\tcl\myCommands.cpp" .\SRC\modelbuilder\tcl\ -Force

# Build
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\build_windows11_full.ps1 `
  -VcpkgRoot .\third_party\vcpkg `
  -MumpsRoot .\third_party\mumps `
  -BuildDir .\build-win11
```

## 7) Fast Triage

| Symptom | Fix |
|---------|-----|
| `ifx` not found after `setvars.bat` | `VS2022INSTALLDIR` not set. Use `init_oneapi_windows11.cmd` or set it manually to your VS2022 install path |
| Configure fails on Tcl | vcpkg should provide it. Check `vcpkg.json` exists and triplet is `x64-windows-static` |
| `init.tcl` not found at runtime | Use `OpenSees-Launch.cmd` (sets `TCL_LIBRARY`) |
| SP/MP fails at `MPI_Init` | Run `mpiexec -n 2 hostname` first. If it fails, run `fix_intel_mpi_windows11.ps1` |
| MUMPS link errors | Remove `-SkipMumps` to let the build script compile MUMPS with `ifx` |

## 8) Intel MPI Repair (Admin)

When `mpiexec -n 2 hostname` fails on Windows client OS:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\fix_intel_mpi_windows11.ps1 -Ranks 2 -SetWinRmAutomatic
```

## 9) Requirements Summary

- Windows 11 (x64)
- Visual Studio 2022 with C++ workload
- Intel oneAPI Base Toolkit (MKL) + HPC Toolkit (ifx, Intel MPI)
- CMake >= 3.29, Ninja, Python 3.11, Git

All installed automatically by `1_install_dependencies.ps1`.

## 10) Notes

- Python 3.11 is pinned for `OpenSeesPy` ABI compatibility. Tcl executables do not need Python at runtime.
- MPI launchers assume Intel MPI on Windows.
- The `CMakeLists.txt` is cross-platform -- it also works on Linux/macOS (Conan or system packages).
- Tcl is provided by vcpkg. No system Tcl install is needed.
- MUMPS is pinned to commit `ec5f340` for reproducibility.
