# OpenSees Windows 11 Runbook (This Repository)

Last updated: 2026-03-27

## 1) Scope

This runbook is for building this repository on Windows 11 x64 with:

- `OpenSees`
- `OpenSeesPy` (`opensees.pyd`)
- `OpenSeesSP`
- `OpenSeesMP`

Build policy in this repo:

- MSVC (`cl`) for C/C++ + Intel `ifx` for Fortran
- vcpkg manifest mode (`x64-windows-static`) with pinned versions:
  - `tcl` `8.6.10-3`
  - `hdf5` `1.14.6`
  - `zlib` `1.3.1`
  - `libaec` `1.1.4`
  - `eigen3` `3.4.1#1`
- Intel MKL for LAPACK/ScaLAPACK
- Intel MPI for parallel targets
- no Conan (vcpkg replaces it for the Windows path)

If you are not building this exact repository, use `WINDOWS11_BUILD_ANY_OPENSEES_SOURCE.md` instead.

## 2) Prerequisites

Install everything at once with the automated script (elevated PowerShell):

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\1_install_dependencies.ps1
```

This installs via `winget`: Git, Visual Studio 2022 (C++ workload), Intel oneAPI Base + HPC, CMake, Ninja, Python 3.11, and optionally Inno Setup. It verifies each tool after install.

Manual install commands (if you prefer not to use the script):

```powershell
winget install --id Git.Git -e --accept-source-agreements --accept-package-agreements
winget install --id Microsoft.VisualStudio.2022.Community -e --accept-source-agreements --accept-package-agreements --override "--wait --quiet --add Microsoft.VisualStudio.Workload.NativeDesktop --includeRecommended"
winget install --id Intel.OneAPI.BaseToolkit -e --accept-source-agreements --accept-package-agreements
winget install --id Intel.OneAPI.HPCToolkit -e --accept-source-agreements --accept-package-agreements
winget install --id Kitware.CMake -e --accept-source-agreements --accept-package-agreements
winget install --id Ninja-build.Ninja -e --accept-source-agreements --accept-package-agreements
winget install --id Python.Python.3.11 -e --accept-source-agreements --accept-package-agreements
```

Verify the environment:

```cmd
call SCRIPTS\init_oneapi_windows11.cmd
where cl
where ifx
where mpiexec
cmake --version
ninja --version
py -3.11 --version
```

Notes:

- `init_oneapi_windows11.cmd` sets `VS2022INSTALLDIR` before calling `setvars.bat`. This is required on systems where Intel oneAPI does not auto-detect the Visual Studio install path.
- `Intel.OneAPI.HPCToolkit` provides `ifx` and Intel MPI; `Intel.OneAPI.BaseToolkit` provides MKL.
- Python 3.11 is pinned for `OpenSeesPy` ABI compatibility.

## 3) Repository Layout Expected

From repo root:

```
third_party\vcpkg\          (git clone https://github.com/microsoft/vcpkg)
third_party\mumps\          (git clone https://github.com/OpenSees/mumps.git)
SCRIPTS\                    (build scripts)
vcpkg.json                  (dependency manifest)
CMakeLists.txt              (cross-platform build system)
cmake\cmake\                (OpenSeesDependencies.cmake etc.)
```

## 4) Build Commands

### Quick start (uses the numbered wrapper script)

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\3_build.ps1
```

The wrapper auto-detects `third_party\vcpkg` and `third_party\mumps`, initializes the VS + oneAPI environment if needed, and forwards to `build_windows11_full.ps1`.

### Direct call with explicit options

Quick build + quick smoke:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\build_windows11_full.ps1 `
  -VcpkgRoot .\third_party\vcpkg `
  -MumpsRoot .\third_party\mumps `
  -BuildDir .\build-win11 `
  -SkipMumps `
  -SmokeMode quick `
  -SmokeTimeoutSec 600
```

Full smoke (no rebuild):

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\build_windows11_full.ps1 `
  -VcpkgRoot .\third_party\vcpkg `
  -MumpsRoot .\third_party\mumps `
  -BuildDir .\build-win11 `
  -SkipMumps `
  -SkipBuild `
  -SmokeMode full `
  -SmokeTimeoutSec 7200
```

Additional flags:

- `-Parallel N` -- parallel compilation jobs (default: CPU count)
- `-DryRun` -- print commands without executing
- `-SkipTests` -- skip smoke tests

## 5) Output Folder

Primary outputs in `build-win11`:

- `OpenSees.exe`
- `OpenSeesSP.exe`
- `OpenSeesMP.exe`
- `opensees.pyd`
- `OpenSees-Launch.cmd`
- `OpenSeesSP-Launch.cmd`
- `OpenSeesMP-Launch.cmd`
- `lib\tcl8.6\...` (Tcl runtime)
- `build_windows11_full.log`

## 6) How To Launch

Interactive OpenSees prompt (`OpenSees >`):

- Double-click `build-win11\OpenSees-Launch.cmd`

Parallel executables:

- Do not double-click `OpenSeesSP.exe` or `OpenSeesMP.exe` directly.
- Use `OpenSeesSP-Launch.cmd` and `OpenSeesMP-Launch.cmd` (they call `mpiexec`).

Optional rank override for SP/MP launchers:

- set `OPS_MPI_N` before launching (default is `2`)

## 7) MPI Health Check

Quick check:

```cmd
call SCRIPTS\init_oneapi_windows11.cmd
mpiexec -n 2 hostname
```

If this fails, `OpenSeesSP`/`OpenSeesMP` launchers will also fail.

## 8) MPI Repair (Admin)

If MPI health check fails, run this once in an elevated PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\fix_intel_mpi_windows11.ps1 -Ranks 2 -SetWinRmAutomatic
```

This configures/starts WinRM, installs the Intel MPI Hydra service, and validates with `mpiexec -n 2 hostname`.

## 9) Logs

- `build-win11\build_windows11_full.log`
- `build-win11\OpenSeesSP-launch.log`
- `build-win11\OpenSeesMP-launch.log`

## 10) Installer Packaging

Create distributables from the built `build-win11` folder:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\4_package.ps1
```

Or directly:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\create_el_ladruno_installer.ps1 `
  -BuildDir .\build-win11 `
  -OutputDir .\dist
```

Dedicated installer guide: `WINDOWS11_INSTALLER_EL_LADRUNO.md`

## 11) Fast Triage

| Symptom | Check |
|---------|-------|
| Configure fails | `vcpkg.json` pinning, triplet `x64-windows-static`, oneAPI env (`ifx`, `MKLROOT`) |
| `ifx` not found after `setvars.bat` | `VS2022INSTALLDIR` must be set first -- use `init_oneapi_windows11.cmd` |
| Serial launch warns about `init.tcl` | Use `OpenSees-Launch.cmd` (sets `TCL_LIBRARY`) |
| SP/MP fails at `MPI_Init` | Run `mpiexec -n 2 hostname` first; fix MPI before debugging OpenSees |

## 12) Sharing This Build (Reproducible File Pack)

Another user can reproduce this exact build by running the 4-step workflow:

```powershell
# Step 1: Install prerequisites
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\1_install_dependencies.ps1

# Step 2: Fetch source + harness (point to this repo or any fork)
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\2_fetch_source.ps1 `
  -OpenSeesRepo "https://github.com/YOUR_USER/OpenSees.git"

# Step 3: Build
cd C:\work\OpenSees-src
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\3_build.ps1

# Step 4: Package (optional)
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\4_package.ps1
```

The file-pack repo is: `https://github.com/nmorabowen/opensees_windows_compilation_files`
