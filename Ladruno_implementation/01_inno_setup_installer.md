---
title: Inno Setup .exe installer with venv wiring
project: Ladruno
status: planned
priority: high
owner: nmora
tags:
  - implementation
  - installer
  - distribution
---

# Inno Setup .exe installer with venv wiring

## What

Replace the PS2EXE-wrapped `install.exe` with a proper Inno Setup-based wizard installer that:

- Welcome page with the LADRUNO ASCII banner (or a banner image)
- Standard Inno wizard: license, install dir, components, ready, progress, finish
- A custom page that asks "Wire OpenSeesPy into a Python venv?" with:
  - radio: skip / create new venv / use existing venv
  - directory picker that defaults to `<install-dir>\opensees_venv`
  - python interpreter picker (defaults to "python on PATH")
- Optional component: "Add `bin\` to user PATH"
- Drops `.pth` file pointing at `<install-dir>\bin\` into the chosen venv's `Lib\site-packages\`

Out of scope (revisit later):
- Code signing
- Auto-update / version check
- Uninstaller that removes the .pth from the venv (Inno's default uninstaller only removes what the installer placed under `<install-dir>` — which is what we want anyway; the .pth-targeted venv may live elsewhere)

## Why

The current `install.ps1` works for power users; the PS2EXE-compiled `install.exe` wraps it but has flaky parameter-binding under Windows PowerShell 5.1 — `-VenvPath` silently skips. For non-developer recipients (the people we'd ship to) we want a real wizard with directory pickers and clear checkboxes, not a console with prompts.

Inno Setup is the standard Windows installer toolkit, free, scriptable in Pascal Script, well-documented. Single-file `setup.exe` output. Recipients are used to its UI from countless other installers.

## Where

- New file: `Ladruno_scripts/installer.iss` (Inno Setup script)
- New file: `Ladruno_scripts/build_inno_installer.ps1` (driver — runs `iscc` against the .iss, output to `Ladruno_files/`)
- Reference: existing [`Ladruno_scripts/make_installer.ps1`](../Ladruno_scripts/make_installer.ps1) — its end-user logic in the generated `install.ps1` is the reference for what the Pascal Script in the .iss should do at install time
- Reference: existing [`Ladruno_scripts/wire_pyenv.ps1`](../Ladruno_scripts/wire_pyenv.ps1) — the venv-creation + `.pth`-writing logic to port into the Inno script

Build: does NOT need a new dependency in `OpenSees/CMakeLists.txt`. Just needs Inno Setup compiler (`iscc.exe`) installed locally.

## How

### Tooling

```
winget install JRSoftware.InnoSetup
```

Installs to `C:\Program Files (x86)\Inno Setup 6\iscc.exe`. Add that dir to PATH in `Ladruno_scripts\build_inno_installer.ps1` or call by absolute path.

### `installer.iss` skeleton

```pascal
; Ladruno OpenSees installer
[Setup]
AppId={{8C8E2E87-LADRUNO-OPENSEES-XXXX}
AppName=Ladruno OpenSees
AppVersion={#GetEnv("LADRUNO_VERSION")}        ; passed via /D from build script
AppPublisher=Ladruno
DefaultDirName={localappdata}\Ladruno\OpenSees
DefaultGroupName=Ladruno OpenSees
OutputDir={#GetEnv("LADRUNO_OUT")}
OutputBaseFilename=Ladruno_OpenSees_{#GetEnv("LADRUNO_VERSION")}_setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog

[Files]
; Bundle every file under dist/ into the installer
Source: "{#GetEnv("LADRUNO_DIST")}\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Tasks]
Name: "addtopath"; Description: "Add {app}\bin to user PATH"; Flags: unchecked

[Code]
var
  VenvPage: TInputDirWizardPage;
  VenvOptionPage: TInputOptionWizardPage;
  PythonPathPage: TInputFileWizardPage;

procedure InitializeWizard;
begin
  VenvOptionPage := CreateInputOptionPage(
    wpSelectComponents,
    'Wire OpenSeesPy',
    'Make `import opensees` work in a Python virtualenv',
    'opensees.pyd is built for Python 3.11. Choose what to do:',
    True, False);
  VenvOptionPage.Add('Skip - I''ll wire it up later');
  VenvOptionPage.Add('Create a new venv at <install-dir>\opensees_venv');
  VenvOptionPage.Add('Use an existing venv (pick directory on next page)');
  VenvOptionPage.SelectedValueIndex := 1;

  VenvPage := CreateInputDirPage(
    VenvOptionPage.ID,
    'Existing virtualenv path',
    'Pick the venv you want to wire OpenSeesPy into',
    'Setup will write a .pth file into this venv''s Lib\site-packages.',
    False, '');
  VenvPage.Add('');

  PythonPathPage := CreateInputFilePage(
    VenvOptionPage.ID,
    'Base Python (for new venv)',
    'Pick the python.exe used to bootstrap the new venv',
    'Must be Python 3.11 - opensees.pyd is built against that exact ABI.');
  PythonPathPage.Add('python.exe', 'Python interpreter|python.exe', 'exe');
end;

function ShouldSkipPage(PageID: Integer): Boolean;
begin
  Result := False;
  if (PageID = VenvPage.ID) and (VenvOptionPage.SelectedValueIndex <> 2) then Result := True;
  if (PageID = PythonPathPage.ID) and (VenvOptionPage.SelectedValueIndex <> 1) then Result := True;
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  pthPath, sitePackages, venvPath, basePython, binDir: String;
  ResultCode: Integer;
begin
  if CurStep = ssPostInstall then
  begin
    binDir := ExpandConstant('{app}\bin');
    case VenvOptionPage.SelectedValueIndex of
      1: begin
           venvPath   := ExpandConstant('{app}\opensees_venv');
           basePython := PythonPathPage.Values[0];
           if basePython = '' then basePython := 'python';
           Exec(basePython, '-m venv "' + venvPath + '"', '', SW_HIDE, ewWaitUntilTerminated, ResultCode);
         end;
      2: begin
           venvPath := VenvPage.Values[0];
         end;
    end;

    if (venvPath <> '') and FileExists(venvPath + '\Scripts\python.exe') then
    begin
      // Ask the venv where its site-packages is, then write the .pth
      Exec(venvPath + '\Scripts\python.exe',
        '-c "import sysconfig,sys; open(sysconfig.get_paths()[''purelib'']+''/ladruno_opensees.pth'',''w'').write(r''' + binDir + ''')"',
        '', SW_HIDE, ewWaitUntilTerminated, ResultCode);
    end;
  end;
end;

[Run]
; Optional post-install: open the install dir
Filename: "{app}"; Description: "Open install folder"; Flags: postinstall shellexec skipifsilent
```

Skeleton above is illustrative — needs a clean pass to:
- Validate every Pascal-Script syntax detail (Inno's Pascal is finicky)
- Switch to a more robust `.pth`-writer (write through a small helper script that ships with the installer rather than constructing a Python one-liner with `''`-escaped paths)
- Embed the LADRUNO ASCII banner. Two reasonable options:
  - Use `WizardImageFile` + a rendered PNG of the ASCII art for a graphical splash (best UX)
  - Use `Setup.WelcomeLabel2` + a custom page with the raw text in a monospace `TNewMemo`

### `build_inno_installer.ps1` driver

```powershell
param([string]$Version = (Get-Date -Format "yyyyMMdd"))

$root = Split-Path -Parent $PSScriptRoot
$dist = Join-Path $root "dist"
$out  = Join-Path $root "Ladruno_files"
$iscc = "C:\Program Files (x86)\Inno Setup 6\iscc.exe"

if (-not (Test-Path $iscc)) {
    Write-Error "Inno Setup not found at $iscc. Install: winget install JRSoftware.InnoSetup"
    exit 1
}

$env:LADRUNO_VERSION = $Version
$env:LADRUNO_DIST    = $dist
$env:LADRUNO_OUT     = $out

& $iscc /Q (Join-Path $PSScriptRoot "installer.iss")
```

Output: `Ladruno_files\Ladruno_OpenSees_<version>_setup.exe` — a single self-contained installer file (no separate zip alongside).

## Risks / open questions

> [!question]
> Banner rendering — do we go with the Pascal-Script + monospace TNewMemo route (preserves the exact ASCII art chars), or render the art to a 410x300 PNG at build time and use it as `WizardImageFile`? PNG is more polished but loses the "made of ASCII chars" identity.

> [!question]
> Should the installer be per-user (`PrivilegesRequired=lowest`, default install to `%LOCALAPPDATA%`) or machine-wide (`admin`, `%PROGRAMFILES%`)? Per-user is safer and more typical for tool installers; machine-wide is required if multiple Windows users on the same box need it.

> [!question]
> Bundle Python? The installer currently relies on the user already having Python 3.11 installed. We could bundle the embeddable Python distribution (~10 MB) and use it as the base for the venv. Avoids one external dependency for non-developer recipients but adds size and complexity. Probably skip for v1.

- Inno Setup is licensed but free for any use including commercial. Confirm the license blurb in `[Setup]LicenseFile=` is what we want shown to recipients, or omit the license page entirely.
- `iscc` exit codes need to be checked in `build_inno_installer.ps1` — silent failures are easy.
- If recipients have Windows Defender SmartScreen on, an unsigned .exe shows the "Unrecognized app" warning. Long-term we want code signing; short-term, document that they should click "Run anyway".
- `Exec` calls in `[Code]` Pascal section run synchronously but capture no stdout. If the venv creation fails we silently skip — log to a file under `{app}\install.log` so the user can troubleshoot.

## Implementation log

*(to be filled in once the next session starts work)*

## Handoff state from previous session (2026-05-01)

The PS2EXE-based attempt is already in place but the .exe form's parameter binding is flaky. Don't spend time fixing it; replace with the Inno Setup path described above. Keep `Ladruno_scripts/make_installer.ps1` and the generated `install.ps1` as a fallback for power users who want a scriptable installer.

Existing artifacts to reuse:
- `Ladruno_scripts/banner_ASCII.txt` — single source of truth for the splash art
- `Ladruno_scripts/wire_pyenv.ps1` — already-tested venv-wiring logic; port the meaningful bits into `installer.iss`
- `dist/` — packaged binaries; the Inno installer's `[Files]` section just bundles everything under here
- `Ladruno_files/install.ps1` — keep as the scriptable / silent-mode fallback. Document in the README that `setup.exe` is for end users and `install.ps1 -InstallPath ... -VenvPath ... -Yes` is for automation
