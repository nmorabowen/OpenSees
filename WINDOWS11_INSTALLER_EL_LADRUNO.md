# El Ladruno OpenSees Installer Note (Windows 11)

Last updated: 2026-03-27

## 1) Purpose

Create a redistributable installer from an already compiled OpenSees build (`build-win11`).

Default branding:

- Product name: `El Ladruno OpenSees`
- Publisher: `El Ladruno`

This document covers packaging only. For building OpenSees first, see:

- `WINDOWS11_BUILD_RUNBOOK.md` -- this repo
- `WINDOWS11_BUILD_ANY_OPENSEES_SOURCE.md` -- any other OpenSees source tree

## 2) Prerequisites

A successful build must exist with at least:

- `build-win11\OpenSees.exe`
- `build-win11\OpenSeesSP.exe`
- `build-win11\OpenSeesMP.exe`
- `build-win11\opensees.pyd`
- `build-win11\lib\...` (Tcl runtime)

If you do not have these, run the build first.

## 3) Quick Command

Using the numbered wrapper:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\4_package.ps1
```

Or directly:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\create_el_ladruno_installer.ps1 `
  -BuildDir .\build-win11 `
  -OutputDir .\dist
```

## 4) Outputs

Generated in `dist\`:

- `ElLadrunoOpenSeesSetup_<version>_x64.exe` -- Windows installer (if Inno Setup is available)
- `ElLadrunoOpenSees_<version>_portable.zip` -- portable package
- `ElLadrunoOpenSees.iss` -- Inno Setup script
- `ElLadrunoOpenSees-stage\` -- staging directory

## 5) Customization

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\4_package.ps1 `
  -AppName "My Custom OpenSees" `
  -Publisher "My Lab" `
  -AppVersion "3.8.0"
```

Other flags:

| Flag | Effect |
|------|--------|
| `-SkipInnoCompile` | Only stage + zip + `.iss` (no installer exe) |
| `-SkipZip` | Skip portable zip |
| `-IncludeAllExamples` | Include full `EXAMPLES` tree |
| `-InnoCompilerPath "C:\...\ISCC.exe"` | Explicit Inno compiler path |

## 6) Inno Setup

If `ISCC.exe` is not found, the script still produces the staging directory, portable zip, and `.iss` script. To install Inno Setup:

```powershell
winget install --id JRSoftware.InnoSetup -e --accept-source-agreements --accept-package-agreements
```

## 7) End-User Installation

For end users who receive the installer:

1. Run `ElLadrunoOpenSeesSetup_<version>_x64.exe`.
2. Launch from Start Menu shortcut `El Ladruno OpenSees`.

Runtime notes:

- Serial `OpenSees` runs from packaged files (no additional dependencies).
- `OpenSeesSP` and `OpenSeesMP` require a working Intel MPI environment (`mpiexec`) on the target machine.

## 8) Validation After Install

```cmd
OpenSees-Launch.cmd
```

For MPI health:

```cmd
mpiexec -n 2 hostname
```

If MPI fails, install/repair Intel MPI before testing SP/MP launchers.
