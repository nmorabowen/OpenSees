# Windows 11 Guide: Compile Any OpenSees Source Tree

Last updated: 2026-03-10

## 1) Goal

Use this when you have a different OpenSees source tree (fork, branch, or fresh clone) and want a repeatable Windows 11 build.

What this guide does in practice:

1. Install the required Windows tools.
2. Download the OpenSees source tree you want to compile.
3. Download the Windows file-pack repo that contains the build harness used in this project.
4. Copy those Windows-specific files into the OpenSees source tree.
5. Run the build command.
6. Optionally package the result as an installer or portable zip.

Short version:

- this guide is not just "run one script"
- you are overlaying this project's Windows build files on top of another OpenSees source tree
- after that, you compile from the target repo root

## 2) Required Tools

- Visual Studio 2022 with C++ workload
- Intel oneAPI Base + HPC (`ifx`, MKL, Intel MPI)
- CMake `>= 3.29`
- Ninja
- Python 3.11 x64
- Git

## 3) Install Prerequisites (Copy/Paste)

Run these in an elevated PowerShell on a clean Windows 11 machine:

```powershell
winget install --id Git.Git -e --accept-source-agreements --accept-package-agreements
winget install --id Microsoft.VisualStudio.2022.Community -e --accept-source-agreements --accept-package-agreements --override "--wait --quiet --add Microsoft.VisualStudio.Workload.NativeDesktop --includeRecommended"
winget install --id Intel.OneAPI.BaseToolkit -e --accept-source-agreements --accept-package-agreements
winget install --id Intel.OneAPI.HPCToolkit -e --accept-source-agreements --accept-package-agreements
winget install --id Kitware.CMake -e --accept-source-agreements --accept-package-agreements
winget install --id Ninja-build.Ninja -e --accept-source-agreements --accept-package-agreements
winget install --id Python.Python.3.11 -e --accept-source-agreements --accept-package-agreements
```

Optional packaging tool:

```powershell
winget install --id JRSoftware.InnoSetup -e --accept-source-agreements --accept-package-agreements
```

Recommended verification after install:

```powershell
git --version
cmake --version
ninja --version
py -3.11 --version
cmd /c "\"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe\" -latest -products * -property installationPath"
cmd /c "set VS2022INSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools && call \"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat\" intel64 && where cl && where ifx && where mpiexec"
```

If your Visual Studio 2022 edition is not `BuildTools`, replace `VS2022INSTALLDIR` with the actual install path returned by `vswhere`.

If `winget` is blocked by policy, install the same products manually and then rerun the verification commands above.

## 4) Get Windows File Pack

Clone the published Windows file-pack repo first:

```powershell
$FilePackRoot = "C:\work\opensees_windows_compilation_files"
git clone https://github.com/nmorabowen/opensees_windows_compilation_files.git $FilePackRoot
```

This repo contains the Windows-specific files that make this workflow reproducible on another machine or another OpenSees branch.

## 5) Create Working Tree

Example:

```powershell
mkdir C:\work
cd C:\work
git clone <YOUR_OPENSEES_REPO_URL> OpenSees-src
cd OpenSees-src
mkdir third_party
git clone https://github.com/microsoft/vcpkg third_party\vcpkg
git clone https://github.com/OpenSees/mumps.git third_party\mumps
```

## 6) Copy Build Harness Files Into Target Tree

Copy the required files from the Windows file-pack repo into the target source tree.

```powershell
New-Item -ItemType Directory -Force -Path .\SCRIPTS, .\cmake\cmake | Out-Null

Copy-Item "$FilePackRoot\SCRIPTS\build_windows11_full.ps1" .\SCRIPTS\build_windows11_full.ps1 -Force
Copy-Item "$FilePackRoot\SCRIPTS\init_oneapi_windows11.cmd" .\SCRIPTS\init_oneapi_windows11.cmd -Force
Copy-Item "$FilePackRoot\SCRIPTS\fix_intel_mpi_windows11.ps1" .\SCRIPTS\fix_intel_mpi_windows11.ps1 -Force
Copy-Item "$FilePackRoot\SCRIPTS\create_el_ladruno_installer.ps1" .\SCRIPTS\create_el_ladruno_installer.ps1 -Force
Copy-Item "$FilePackRoot\vcpkg.json" .\vcpkg.json -Force
Copy-Item "$FilePackRoot\CMakeLists.txt" .\CMakeLists.txt -Force
Copy-Item "$FilePackRoot\cmake\cmake\OpenSeesDependenciesWin.cmake" .\cmake\cmake\OpenSeesDependenciesWin.cmake -Force
```

This step is the important part.

You are taking the Windows build harness from `opensees_windows_compilation_files` and placing it into the OpenSees source tree you want to compile.

That means:

- the target repo keeps its source code
- the Windows-specific build scripts come from the file-pack repo
- the copied `CMakeLists.txt` and `cmake\cmake\OpenSeesDependenciesWin.cmake` become the Windows porting layer used by this build flow

Optional: copy the notes themselves into the target tree too:

```powershell
Copy-Item "$FilePackRoot\WINDOWS11_BUILD_RUNBOOK.md" .\WINDOWS11_BUILD_RUNBOOK.md -Force
Copy-Item "$FilePackRoot\WINDOWS11_BUILD_ANY_OPENSEES_SOURCE.md" .\WINDOWS11_BUILD_ANY_OPENSEES_SOURCE.md -Force
Copy-Item "$FilePackRoot\WINDOWS11_INSTALLER_EL_LADRUNO.md" .\WINDOWS11_INSTALLER_EL_LADRUNO.md -Force
```

If the target source does not support the same CMake options yet, the copied `CMakeLists.txt` and `cmake\cmake\OpenSeesDependenciesWin.cmake` are the Windows-specific porting layer this workflow expects.

## 7) Build In Target Tree

From target repo root:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\build_windows11_full.ps1 `
  -VcpkgRoot .\third_party\vcpkg `
  -MumpsRoot .\third_party\mumps `
  -BuildDir .\build-win11 `
  -SkipMumps `
  -SmokeMode quick `
  -SmokeTimeoutSec 600
```

If MUMPS is not built yet, remove `-SkipMumps`.

If you only want the Tcl executables (`OpenSees.exe`, `OpenSeesSP.exe`, `OpenSeesMP.exe`), Python is not used at runtime.

If you also want `OpenSeesPy`, this workflow currently expects Python `3.11` during build because `opensees.pyd` is compiled against that Python ABI.

## 8) Artifacts You Should Get

- `build-win11\OpenSees.exe`
- `build-win11\OpenSeesSP.exe`
- `build-win11\OpenSeesMP.exe`
- `build-win11\opensees.pyd`
- launcher files in `build-win11\`

## 9) How To Run

Interactive serial prompt:

- `build-win11\OpenSees-Launch.cmd`

MPI runs:

- `build-win11\OpenSeesSP-Launch.cmd`
- `build-win11\OpenSeesMP-Launch.cmd`

Or manual MPI:

```powershell
cd EXAMPLES\ParallelModelMP
mpiexec -n 2 ..\..\build-win11\OpenSeesMP.exe exampleMP.tcl
```

If Intel oneAPI prints a Visual Studio detection warning on the target machine, initialize with:

```cmd
call SCRIPTS\init_oneapi_windows11.cmd
where cl
where ifx
where mpiexec
```

## 10) Fast Triage

If configure fails:

- check `vcpkg.json` pinning and installed triplet (`x64-windows-static`)
- check oneAPI initialized environment (`ifx`, `MKLROOT`, `I_MPI_ROOT`)

If serial launch warns about `init.tcl`:

- use `OpenSees-Launch.cmd` (sets `TCL_LIBRARY`)

If MP/SP fails at `MPI_Init` or `HYD_spawn`:

- run `mpiexec -n 2 hostname` first
- solve Intel MPI bootstrap/service issue before debugging OpenSees itself

## 11) Intel MPI Repair (Admin)

When `mpiexec -n 2 hostname` fails on Windows client OS, run this from an elevated PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\fix_intel_mpi_windows11.ps1 -Ranks 2 -SetWinRmAutomatic
```

If your target repo does not contain this script yet, copy it from this repository:

- `SCRIPTS\fix_intel_mpi_windows11.ps1`

## 12) Build Installer (El Ladruno OpenSees)

After successful build in the target tree, create a redistributable package:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\create_el_ladruno_installer.ps1 `
  -BuildDir .\build-win11 `
  -OutputDir .\dist
```

This generates:

- `dist\ElLadrunoOpenSees_<version>_portable.zip`
- `dist\ElLadrunoOpenSees.iss`

If Inno Setup is available (`ISCC.exe`), the script can also compile a Windows installer executable.

Dedicated installer guide:

- `WINDOWS11_INSTALLER_EL_LADRUNO.md`

## 13) File Pack Required In Any Target Source Tree

The bootstrap in sections `4)` and `6)` copies these files from the `opensees_windows_compilation_files` repo and preserves the required destination paths:

- `SCRIPTS\build_windows11_full.ps1` -> target repo `SCRIPTS\`
- `SCRIPTS\init_oneapi_windows11.cmd` -> target repo `SCRIPTS\`
- `SCRIPTS\fix_intel_mpi_windows11.ps1` -> target repo `SCRIPTS\`
- `SCRIPTS\create_el_ladruno_installer.ps1` -> target repo `SCRIPTS\`
- `vcpkg.json` -> target repo root
- `CMakeLists.txt` -> target repo root
- `cmake\cmake\OpenSeesDependenciesWin.cmake` -> target repo `cmake\cmake\`

If someone asks "what do I actually do?", the answer is:

1. Install prerequisites.
2. Clone the OpenSees repo you want to build.
3. Clone `https://github.com/nmorabowen/opensees_windows_compilation_files`.
4. Copy the file-pack files into the target repo.
5. Run `SCRIPTS\build_windows11_full.ps1`.
