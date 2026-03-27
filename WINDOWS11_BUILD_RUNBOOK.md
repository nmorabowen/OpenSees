# OpenSees Windows 11 Runbook (This Repository)

Last updated: 2026-03-09

## 1) Scope

This runbook is for building this repository on Windows 11 x64 with:

- `OpenSees`
- `OpenSeesPy` (`opensees.pyd`)
- `OpenSeesSP`
- `OpenSeesMP`

Build policy in this repo:

- MSVC + Intel oneAPI `ifx`
- no Conan
- vcpkg manifest mode with pinned versions:
  - `tcl` `8.6.10-3`
  - `hdf5` `1.14.6`
  - `zlib` `1.3.1`
  - `libaec` `1.1.4`
  - `eigen3` `3.4.1#1`

What this runbook assumes:

- you are already inside this repository
- this repository already contains the Windows build harness files
- after prerequisites are installed, you only need the third-party folders and the build command

If you are not building this exact repository, use `WINDOWS11_BUILD_ANY_OPENSEES_SOURCE.md` first. That guide explains how to copy the required Windows file-pack files into another OpenSees source tree before compiling.

## 2) Prerequisites

Install:

- Visual Studio 2022 with C++ tools
- Intel oneAPI Base + HPC (ifx, MKL, Intel MPI)
- CMake `>= 3.29`
- Ninja
- Python 3.11 x64

PowerShell bootstrap commands:

```powershell
winget install --id Microsoft.VisualStudio.2022.Community -e --accept-source-agreements --accept-package-agreements --override "--wait --quiet --add Microsoft.VisualStudio.Workload.NativeDesktop --includeRecommended"
winget install --id Intel.OneAPI.BaseToolkit -e --accept-source-agreements --accept-package-agreements
winget install --id Intel.OneAPI.HPCToolkit -e --accept-source-agreements --accept-package-agreements
winget install --id Kitware.CMake -e --accept-source-agreements --accept-package-agreements
winget install --id Ninja-build.Ninja -e --accept-source-agreements --accept-package-agreements
winget install --id Python.Python.3.11 -e --accept-source-agreements --accept-package-agreements
```

Recommended verification:

```powershell
cmd /c "\"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat\" intel64 && where ifx && where mpiexec"
cmake --version
ninja --version
py -3.11 --version
```

If `setvars.bat` warns that Visual Studio was not found, use the wrapper in this repo instead:

```cmd
call SCRIPTS\init_oneapi_windows11.cmd
where cl
where ifx
where mpiexec
```

Notes:

- The Visual Studio command installs Community 2022 with the native desktop C++ workload.
- `Intel.OneAPI.HPCToolkit` provides `ifx` and Intel MPI; `Intel.OneAPI.BaseToolkit` provides MKL and base tooling.
- If `winget` is blocked by policy on a machine, install the same products manually and keep the required versions compatible with this runbook.

## 3) Repository Layout Expected

From repo root:

- `third_party\vcpkg`
- `third_party\mumps`
- `SCRIPTS\build_windows11_full.ps1`
- `vcpkg.json`

## 4) Build Commands

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

Packaging only (regenerate launchers and copied runtimes):

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\build_windows11_full.ps1 `
  -VcpkgRoot .\third_party\vcpkg `
  -MumpsRoot .\third_party\mumps `
  -BuildDir .\build-win11 `
  -SkipMumps `
  -SkipBuild `
  -SkipTests
```

## 5) Output Folder

Primary outputs in `build-win11`:

- `OpenSees.exe`
- `OpenSeesSP.exe`
- `OpenSeesMP.exe`
- `opensees.pyd`
- `OpenSees-Launch.cmd`
- `OpenSeesSP-Launch.cmd`
- `OpenSeesMP-Launch.cmd`
- `lib\tcl8.6\...`
- `build_windows11_full.log`

## 6) How To Launch

Interactive OpenSees prompt (`OpenSees >`):

- Double-click `build-win11\OpenSees-Launch.cmd`

Parallel executables:

- Do not double-click `OpenSeesSP.exe` or `OpenSeesMP.exe` directly.
- Use `OpenSeesSP-Launch.cmd` and `OpenSeesMP-Launch.cmd` (they call `mpiexec`).

Optional rank override for SP/MP launchers:

- set `OPS_MPI_N` before launching (default is `2`)

## 7) Current MPI Status (This Machine)

Latest local check (2026-03-08):

- `mpiexec -n 2 hostname` passes
- quick smoke (`-SmokeMode quick`) passes end-to-end

Quick health check:

```cmd
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
mpiexec -n 2 hostname
```

To avoid the Visual Studio detection warning on this machine, prefer:

```cmd
call SCRIPTS\init_oneapi_windows11.cmd
mpiexec -n 2 hostname
```

If this fails, `OpenSeesSP`/`OpenSeesMP` launchers will also fail.

## 8) MPI Repair (Admin)

If MPI health check fails, run this once in an elevated PowerShell ("Run as Administrator"):

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\fix_intel_mpi_windows11.ps1 -Ranks 2 -SetWinRmAutomatic
```

This script:

- configures/starts WinRM
- installs and starts Intel MPI Hydra service
- validates with `mpiexec -n 2 hostname`

## 9) Logs For Troubleshooting

- `build-win11\build_windows11_full.log`
- `build-win11\OpenSeesSP-launch.log`
- `build-win11\OpenSeesMP-launch.log`

## 10) Installer Packaging (El Ladruno OpenSees)

Create distributables from the built `build-win11` folder:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\create_el_ladruno_installer.ps1 `
  -BuildDir .\build-win11 `
  -OutputDir .\dist
```

Outputs:

- `dist\ElLadrunoOpenSees_3.7.2_portable.zip` (portable package)
- `dist\ElLadrunoOpenSees.iss` (Inno Setup script)

If Inno Setup (`ISCC.exe`) is installed, rerun without `-SkipInnoCompile` to produce a `.exe` installer.

Dedicated installer guide:

- `WINDOWS11_INSTALLER_EL_LADRUNO.md`

## 11) Reproducible File Pack (What To Share)

If someone else will reproduce this exact Windows 11 flow, they can clone the published Windows file-pack repo and copy the files automatically:

```powershell
$FilePackRoot = "C:\work\opensees_windows_compilation_files"
git clone https://github.com/nmorabowen/opensees_windows_compilation_files.git $FilePackRoot

New-Item -ItemType Directory -Force -Path .\SCRIPTS, .\cmake\cmake | Out-Null

Copy-Item "$FilePackRoot\SCRIPTS\build_windows11_full.ps1" .\SCRIPTS\build_windows11_full.ps1 -Force
Copy-Item "$FilePackRoot\SCRIPTS\init_oneapi_windows11.cmd" .\SCRIPTS\init_oneapi_windows11.cmd -Force
Copy-Item "$FilePackRoot\SCRIPTS\fix_intel_mpi_windows11.ps1" .\SCRIPTS\fix_intel_mpi_windows11.ps1 -Force
Copy-Item "$FilePackRoot\SCRIPTS\create_el_ladruno_installer.ps1" .\SCRIPTS\create_el_ladruno_installer.ps1 -Force
Copy-Item "$FilePackRoot\vcpkg.json" .\vcpkg.json -Force
Copy-Item "$FilePackRoot\CMakeLists.txt" .\CMakeLists.txt -Force
Copy-Item "$FilePackRoot\cmake\cmake\OpenSeesDependenciesWin.cmake" .\cmake\cmake\OpenSeesDependenciesWin.cmake -Force
```

Those commands preserve the required destination paths:

- `SCRIPTS\build_windows11_full.ps1` -> target repo `SCRIPTS\`
- `SCRIPTS\init_oneapi_windows11.cmd` -> target repo `SCRIPTS\`
- `SCRIPTS\fix_intel_mpi_windows11.ps1` -> target repo `SCRIPTS\`
- `SCRIPTS\create_el_ladruno_installer.ps1` -> target repo `SCRIPTS\`
- `vcpkg.json` -> target repo root
- `CMakeLists.txt` -> target repo root
- `cmake\cmake\OpenSeesDependenciesWin.cmake` -> target repo `cmake\cmake\`

Recommended package name for sharing:

- `OpenSees_Win11_Repro_FilePack.zip`

In plain terms, the sharing workflow is:

1. Another user installs the prerequisites.
2. They clone their OpenSees source tree.
3. They clone `opensees_windows_compilation_files`.
4. They copy the file-pack files into that source tree.
5. They run the same build command from this runbook.
