# El Ladruno OpenSees Installer Note (Windows 11)

Last updated: 2026-03-10

## 1) Purpose

Create a redistributable installer from an already compiled OpenSees build (`build-win11`).

Branding used by default:

- Product name: `El Ladruno OpenSees`
- Publisher: `El Ladruno`

What this note is for:

- first you compile OpenSees successfully
- then you package that compiled result into a portable zip and optionally a Windows installer

This document does not replace the build guide.

Use:

- `WINDOWS11_BUILD_RUNBOOK.md` if you are building this repository
- `WINDOWS11_BUILD_ANY_OPENSEES_SOURCE.md` if you are building a different OpenSees source tree and need to copy the Windows file-pack files first

## 2) Inputs Required

You must already have a successful build in:

- `build-win11\OpenSees.exe`
- `build-win11\OpenSeesSP.exe`
- `build-win11\OpenSeesMP.exe`
- `build-win11\opensees.pyd`
- `build-win11\lib\...` (Tcl runtime)

If you do not have those files yet, stop here and run the build workflow first.

## 3) Main Command

From repository root:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\create_el_ladruno_installer.ps1 `
  -BuildDir .\build-win11 `
  -OutputDir .\dist
```

## 4) Outputs

Generated in `dist\`:

- `ElLadrunoOpenSeesSetup_<version>_x64.exe` (Windows installer, if Inno Setup compiler is available)
- `ElLadrunoOpenSees_<version>_portable.zip` (portable package)
- `ElLadrunoOpenSees.iss` (Inno Setup script)
- `ElLadrunoOpenSees-stage\` (staging directory used for packaging)

## 5) Optional Flags

Useful options:

- `-SkipInnoCompile`: only stage + zip + `.iss` (no installer compilation)
- `-SkipZip`: skip portable zip generation
- `-IncludeAllExamples`: include full `EXAMPLES` tree
- `-IncludeExamples:$false`: package without examples
- `-AppVersion "3.7.2"`: force installer version
- `-AppName "El Ladruno OpenSees"`: change product name
- `-Publisher "El Ladruno"`: change publisher
- `-InnoCompilerPath "C:\Path\To\ISCC.exe"`: explicit Inno compiler path

## 6) Inno Setup Requirement

If `ISCC.exe` is not found, the script still produces:

- staging directory
- portable zip
- `.iss` script

To install Inno Setup using winget:

```powershell
winget install --id JRSoftware.InnoSetup -e --accept-source-agreements --accept-package-agreements
```

## 7) Install On Another Machine

For end users:

1. Run `ElLadrunoOpenSeesSetup_<version>_x64.exe`.
1. Launch from Start Menu shortcut `El Ladruno OpenSees`.

Important runtime note:

- Serial OpenSees runs from packaged files.
- `OpenSeesSP` and `OpenSeesMP` require a working MPI environment (`mpiexec`) on the target machine.

That means the installer packages the OpenSees binaries and runtime files, but MPI execution still depends on Intel MPI being installed and working on the destination machine.

## 8) Recommended Validation After Install

On target machine:

```cmd
OpenSees-Launch.cmd
```

For MPI health:

```cmd
mpiexec -n 2 hostname
```

If MPI check fails, fix Intel MPI before testing SP/MP launchers.

## 9) Files This Chapter Depends On

To reproduce installer generation in another source tree, ensure these files exist at these paths:

- `SCRIPTS\create_el_ladruno_installer.ps1`
- `SCRIPTS\build_windows11_full.ps1` (provides launchers/runtime packaging)
- `SCRIPTS\init_oneapi_windows11.cmd` (avoids oneAPI Visual Studio detection issues on some machines)
- `SCRIPTS\fix_intel_mpi_windows11.ps1` (for MPI repair guidance in distributed docs)
- `vcpkg.json`
- `CMakeLists.txt`
- `cmake\cmake\OpenSeesDependenciesWin.cmake`

Why these files matter:

- the scripts in `SCRIPTS\` perform the build, oneAPI initialization, MPI repair, and packaging
- `vcpkg.json` and `CMakeLists.txt` define the Windows build configuration used by this workflow
- `cmake\cmake\OpenSeesDependenciesWin.cmake` provides the Windows dependency logic expected by the copied build harness

Automatic copy from the published Windows file-pack repo:

```powershell
$FilePackRoot = "C:\work\opensees_windows_compilation_files"
git clone https://github.com/nmorabowen/opensees_windows_compilation_files.git $FilePackRoot

New-Item -ItemType Directory -Force -Path .\SCRIPTS, .\cmake\cmake | Out-Null

Copy-Item "$FilePackRoot\SCRIPTS\create_el_ladruno_installer.ps1" .\SCRIPTS\create_el_ladruno_installer.ps1 -Force
Copy-Item "$FilePackRoot\SCRIPTS\build_windows11_full.ps1" .\SCRIPTS\build_windows11_full.ps1 -Force
Copy-Item "$FilePackRoot\SCRIPTS\init_oneapi_windows11.cmd" .\SCRIPTS\init_oneapi_windows11.cmd -Force
Copy-Item "$FilePackRoot\SCRIPTS\fix_intel_mpi_windows11.ps1" .\SCRIPTS\fix_intel_mpi_windows11.ps1 -Force
Copy-Item "$FilePackRoot\vcpkg.json" .\vcpkg.json -Force
Copy-Item "$FilePackRoot\CMakeLists.txt" .\CMakeLists.txt -Force
Copy-Item "$FilePackRoot\cmake\cmake\OpenSeesDependenciesWin.cmake" .\cmake\cmake\OpenSeesDependenciesWin.cmake -Force
```

If someone starts from a fresh OpenSees source download, the full process is:

1. Install prerequisites.
2. Clone the target OpenSees source tree.
3. Clone `https://github.com/nmorabowen/opensees_windows_compilation_files`.
4. Copy the file-pack files shown above into the target tree.
5. Run the build command from the build guide.
6. Confirm `build-win11` was created successfully.
7. Run the installer command from section `3)` of this note.
