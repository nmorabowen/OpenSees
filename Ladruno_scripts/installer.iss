; Ladruno OpenSees - Inno Setup script
;
; Build:
;   $env:LADRUNO_VERSION = '20260501'      ; or any string
;   $env:LADRUNO_DIST    = 'C:\...\dist'   ; abs path to the dist folder
;   $env:LADRUNO_OUT     = 'C:\...\Ladruno_files'
;   & 'C:\Program Files (x86)\Inno Setup 6\iscc.exe' installer.iss
;
; The build_inno_installer.ps1 driver does this for you.
;
; Install-time switches (Inno's own, plus two of ours):
;   /DIR="C:\path"          install location
;   /TASKS=addtopath        append {app}\bin to the user PATH
;   /SILENT | /VERYSILENT   no wizard          /LOG="file"   setup log
;   /VENV=skip              wire no venv  -- the DEFAULT for a silent install
;   /VENV=new               create {app}\opensees_venv and wire it
;   /VENV="C:\path\to\venv" wire an EXISTING venv (the ROOT folder, the one
;                           holding Scripts\ and Lib\ -- NOT python.exe)
;   /BASEPYTHON="C:\...\python.exe"   interpreter to bootstrap /VENV=new (3.12)
;
; The venv switches exist because Inno never shows custom wizard pages under
; /SILENT, so a scripted install used to be stuck with the page defaults
; ("create a venv from the literal string python.exe") -- which fails on any
; machine without such a python on the elevated PATH, and then blocked on a
; modal error box that no one was there to dismiss.
;
;   setup.exe /VERYSILENT /VENV="C:\Users\me\venv\opensees_env"

#define LadrunoVersion GetEnv("LADRUNO_VERSION")
#define LadrunoDist    GetEnv("LADRUNO_DIST")
#define LadrunoOut     GetEnv("LADRUNO_OUT")

#if LadrunoVersion == ""
  #error LADRUNO_VERSION env var must be set
#endif
#if LadrunoDist == ""
  #error LADRUNO_DIST env var must be set
#endif
#if LadrunoOut == ""
  #error LADRUNO_OUT env var must be set
#endif

[Setup]
AppId={{8C8E2E87-1A2B-4C3D-9E4F-0123456789AB}
AppName=Ladruno OpenSees
AppVersion={#LadrunoVersion}
AppPublisher=Ladruno
AppPublisherURL=https://github.com/nmorabowen
DefaultDirName={autopf}\Ladruno\OpenSees
; Inno defaults UsePreviousAppDir to yes: it pre-fills (and can silently skip
; showing, when combined with /DIR=) the destination page from whatever path
; a PRIOR run with this same AppId used, read from the uninstall registry key
; -- including a one-off /DIR= test/scratch run. Since every Ladruno OpenSees
; build shares one fixed AppId, that makes the destination "sticky" across
; unrelated builds/machines: a throwaway test install can silently redirect
; every subsequent real install until something else overwrites the registry
; value. Pin it off so DefaultDirName above is authoritative every time.
UsePreviousAppDir=no
DefaultGroupName=Ladruno OpenSees
DisableProgramGroupPage=yes
OutputDir={#LadrunoOut}
OutputBaseFilename=Ladruno_OpenSees_{#LadrunoVersion}_setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
; Default to a per-machine install (Program Files, requires admin/UAC).
; PrivilegesRequiredOverridesAllowed=dialog lets the user downgrade to a
; per-user install at runtime, in which case {autopf} -> %LOCALAPPDATA%\Programs.
PrivilegesRequired=admin
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesInstallIn64BitMode=x64compatible
ArchitecturesAllowed=x64compatible
UninstallDisplayName=Ladruno OpenSees {#LadrunoVersion}

[Files]
; The payload: everything under dist/ goes under {app}/, preserving subdirs.
Source: "{#LadrunoDist}\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion
; Helper resources extracted to {tmp} only at install time.
Source: "wire_venv_pth.py";        Flags: dontcopy
Source: "banner_ASCII_utf8bom.txt"; Flags: dontcopy

[Tasks]
Name: "addtopath"; Description: "Add {app}\bin to user PATH"; Flags: unchecked

[Registry]
; Append {app}\bin to user PATH when 'addtopath' task is selected.
; Inno broadcasts WM_SETTINGCHANGE after writing under HKCU\Environment.
Root: HKCU; Subkey: "Environment"; ValueType: expandsz; ValueName: "Path"; \
  ValueData: "{olddata};{app}\bin"; \
  Tasks: addtopath; Check: NeedsAddPath(ExpandConstant('{app}\bin'))

[Run]
Filename: "{app}"; Description: "Open install folder"; Flags: postinstall shellexec skipifsilent unchecked

[UninstallDelete]
Type: files; Name: "{app}\install.log"

[Code]
var
  BannerPage:     TWizardPage;
  VenvOptionPage: TInputOptionWizardPage;
  VenvDirPage:    TInputDirWizardPage;
  PythonPathPage: TInputFileWizardPage;

procedure LogToFile(const Msg: String);
var
  LogPath, Existing: String;
  Lst: TStringList;
begin
  LogPath := ExpandConstant('{app}\install.log');
  Lst := TStringList.Create;
  try
    if FileExists(LogPath) then
      Lst.LoadFromFile(LogPath);
    Lst.Add(Msg);
    Lst.SaveToFile(LogPath);
  finally
    Lst.Free;
  end;
end;

function NeedsAddPath(const Param: String): Boolean;
var
  OrigPath: String;
begin
  if not RegQueryStringValue(HKEY_CURRENT_USER, 'Environment', 'Path', OrigPath) then
  begin
    Result := True;
    Exit;
  end;
  Result := Pos(';' + Lowercase(Param) + ';', ';' + Lowercase(OrigPath) + ';') = 0;
end;

procedure InitializeWizard;
var
  Memo:       TNewMemo;
  Lst:        TStringList;
  BannerText: String;
  TmpFile:    String;
begin
  // ---- banner page (immediately after wpWelcome) -----------------------
  BannerPage := CreateCustomPage(wpWelcome,
    'Ladruno OpenSees',
    'OpenSees + HDF5/MPCO + MUMPS, packaged for Windows');

  ExtractTemporaryFile('banner_ASCII_utf8bom.txt');
  TmpFile := ExpandConstant('{tmp}\banner_ASCII_utf8bom.txt');
  Lst := TStringList.Create;
  try
    try
      Lst.LoadFromFile(TmpFile);
      BannerText := Lst.Text;
    except
      BannerText := '';
    end;
  finally
    Lst.Free;
  end;
  if BannerText = '' then
    BannerText := 'Ladruno OpenSees';

  Memo := TNewMemo.Create(BannerPage);
  Memo.Parent := BannerPage.Surface;
  Memo.Left := 0;
  Memo.Top := 0;
  Memo.Width := BannerPage.SurfaceWidth;
  Memo.Height := BannerPage.SurfaceHeight;
  Memo.Anchors := [akLeft, akTop, akRight, akBottom];
  Memo.ReadOnly := True;
  Memo.WordWrap := False;
  Memo.ScrollBars := ssBoth;
  Memo.Font.Name := 'Consolas';
  Memo.Font.Size := 8;
  Memo.Color := clBlack;
  Memo.Font.Color := clAqua;
  Memo.Text := BannerText + #13#10 + #13#10 +
    'Ladruno OpenSees ' + ExpandConstant('{#LadrunoVersion}') + #13#10 +
    'Includes: OpenSees.exe, OpenSeesSP.exe, OpenSeesMP.exe,' + #13#10 +
    '          opensees.pyd  (sequential Python, 3.12),' + #13#10 +
    '          openseesmp.pyd (MPI-parallel Python, self-contained' + #13#10 +
    '                          Intel MPI + mpiexec bundled),' + #13#10 +
    '          HDF5 + MPCO recorder, MUMPS sparse solver.' + #13#10 + #13#10 +
    'The Tcl OpenSeesSP / OpenSeesMP exes need Intel oneAPI''s' + #13#10 +
    'mpiexec on PATH (installed separately); openseesmp ships' + #13#10 +
    'its own matching mpiexec under <install>\openseesmp\.';

  // ---- venv option page (after wpSelectDir so {app} is final) ----------
  VenvOptionPage := CreateInputOptionPage(wpSelectDir,
    'Wire OpenSeesPy into a Python virtualenv',
    'Make ''import opensees'' AND ''import openseesmp'' work without sys.path fiddling',
    'Both modules are built for CPython 3.12. Choose what to do:',
    True, False);
  VenvOptionPage.Add('Skip - I''ll wire it up later');
  VenvOptionPage.Add('Create a new venv at <install-dir>\opensees_venv');
  VenvOptionPage.Add('Use an existing venv (choose path on the next page)');
  VenvOptionPage.SelectedValueIndex := 1;

  // ---- existing venv path page -----------------------------------------
  VenvDirPage := CreateInputDirPage(VenvOptionPage.ID,
    'Existing virtualenv path',
    'Choose the venv root folder (NOT python.exe)',
    'Point at the directory that contains "Scripts\python.exe" and "Lib\site-packages\" - e.g. C:\envs\myproject - not at python.exe itself. Setup writes a .pth file inside that venv so "import opensees" and "import openseesmp" find the binaries.',
    False, '');
  VenvDirPage.Add('Venv root folder');

  // ---- base python picker (only used when creating a new venv) ---------
  PythonPathPage := CreateInputFilePage(VenvOptionPage.ID,
    'Base Python (for new venv)',
    'Pick the python.exe used to bootstrap the new venv',
    'Must be Python 3.12 - opensees.pyd is built against that exact ABI.');
  PythonPathPage.Add('python.exe', 'Python interpreter|python.exe;py.exe', 'exe');
  PythonPathPage.Values[0] := 'python.exe';
end;

function ShouldSkipPage(PageID: Integer): Boolean;
begin
  Result := False;
  if (PageID = VenvDirPage.ID)    and (VenvOptionPage.SelectedValueIndex <> 2) then Result := True;
  if (PageID = PythonPathPage.ID) and (VenvOptionPage.SelectedValueIndex <> 1) then Result := True;
end;

function NextButtonClick(CurPageID: Integer): Boolean;
var
  VenvDir, VenvPython: String;
begin
  Result := True;
  if CurPageID = VenvDirPage.ID then
  begin
    VenvDir := Trim(VenvDirPage.Values[0]);
    if VenvDir = '' then
    begin
      MsgBox('Please pick a venv directory, or go back and choose Skip / Create.', mbError, MB_OK);
      Result := False;
      Exit;
    end;
    VenvPython := AddBackslash(VenvDir) + 'Scripts\python.exe';
    if not FileExists(VenvPython) then
    begin
      MsgBox('Could not find:' + #13#10 + '  ' + VenvPython + #13#10 + #13#10 +
             'That path does not look like a Python virtualenv. Pick another folder.',
             mbError, MB_OK);
      Result := False;
    end;
  end;
end;

// ---- /VENV= and /BASEPYTHON= : venv control for SCRIPTED installs -------
// Inno never shows custom wizard pages under /SILENT|/VERYSILENT, so without
// these a scripted install could only ever take the page DEFAULTS -- which
// meant "create a new venv, bootstrapped from the literal string python.exe".
// In an elevated installer there is usually no such python on PATH, so every
// silent install did the same thing: `python.exe -m venv exited 1`, then a
// MODAL ERROR BOX that blocks forever because nobody is there to click it.
// Measured 2026-08-11: a /VERYSILENT run sat on that box for ~3 minutes.
//
//   /VENV=skip              do not wire anything (the silent DEFAULT)
//   /VENV=new               create <install-dir>\opensees_venv
//   /VENV=<path-to-venv>    wire an EXISTING venv (dir containing Scripts\)
//   /BASEPYTHON=<py.exe>    interpreter used to bootstrap /VENV=new
//
// e.g.  setup.exe /VERYSILENT /VENV="C:\Users\me\venv\opensees_env"
function ParamOrEmpty(const Name: String): String;
begin
  Result := Trim(ExpandConstant('{param:' + Name + '|}'));
end;

// Resolve what to actually do, from (in order) the command line, silent-mode
// policy, then the wizard pages. Returns 0 = skip, 1 = create, 2 = existing.
function ResolveVenvChoice(var VenvDir: String; var BasePython: String): Integer;
var
  P: String;
begin
  VenvDir    := '';
  BasePython := ParamOrEmpty('basepython');
  P          := ParamOrEmpty('venv');

  if P <> '' then
  begin
    if CompareText(P, 'skip') = 0 then
    begin
      Result := 0;
      Exit;
    end;
    if CompareText(P, 'new') = 0 then
    begin
      Result  := 1;
      VenvDir := ExpandConstant('{app}\opensees_venv');
      Exit;
    end;
    Result  := 2;                       // anything else is a venv path
    VenvDir := RemoveBackslash(P);
    Exit;
  end;

  // No /VENV= given. A silent install must NOT guess: creating a venv from a
  // python that probably is not there is how the failure above happened, and a
  // half-done install that logs a warning nobody reads is worse than one that
  // plainly did nothing. Interactive runs still get the wizard's choice.
  if WizardSilent then
  begin
    Result := 0;
    Exit;
  end;

  Result := VenvOptionPage.SelectedValueIndex;
  if Result = 1 then
  begin
    VenvDir := ExpandConstant('{app}\opensees_venv');
    if BasePython = '' then BasePython := Trim(PythonPathPage.Values[0]);
  end;
  if Result = 2 then
    VenvDir := RemoveBackslash(Trim(VenvDirPage.Values[0]));
end;

// Never block a scripted install on a modal box; the log is the channel there.
procedure Complain(const Msg: String);
begin
  LogToFile('  ' + Msg);
  if not WizardSilent then
    MsgBox(Msg, mbError, MB_OK);
end;

procedure WirePth(const VenvPython, BinDir, MpDir: String);
var
  HelperPy: String;
  Params:   String;
  Code:     Integer;
begin
  ExtractTemporaryFile('wire_venv_pth.py');
  HelperPy := ExpandConstant('{tmp}\wire_venv_pth.py');
  Params := AddQuotes(HelperPy) + ' ' + AddQuotes(BinDir) + ' ' + AddQuotes(MpDir);
  if Exec(VenvPython, Params, '', SW_HIDE, ewWaitUntilTerminated, Code) then
  begin
    case Code of
      0: LogToFile('  .pth wired ok via ' + VenvPython);
      3: LogToFile('  WARNING: venv python is not 3.12 - opensees will fail to import');
    else
      LogToFile('  WARNING: wire_venv_pth.py exited ' + IntToStr(Code));
    end;
  end
  else
    LogToFile('  ERROR: failed to launch ' + VenvPython);
end;

procedure CreateAndWireVenv(const BasePython, VenvDir, BinDir, MpDir: String);
var
  Code: Integer;
  VenvPython: String;
begin
  LogToFile('Creating venv at ' + VenvDir + ' (base: ' + BasePython + ') ...');
  if not Exec(BasePython, '-m venv ' + AddQuotes(VenvDir), '', SW_HIDE, ewWaitUntilTerminated, Code) then
  begin
    Complain('ERROR: failed to launch ' + BasePython + '. Skipping venv wiring; ' +
             'run Ladruno_scripts\wire_venv_pth.py later, or re-run with ' +
             '/BASEPYTHON=<path to a 3.12 python.exe>.');
    Exit;
  end;
  if Code <> 0 then
  begin
    Complain('ERROR: ' + BasePython + ' -m venv exited ' + IntToStr(Code) +
             '. Skipping venv wiring; pass /BASEPYTHON=<path to a 3.12 ' +
             'python.exe>, or /VENV=<existing venv> to wire one you already have.');
    Exit;
  end;
  VenvPython := AddBackslash(VenvDir) + 'Scripts\python.exe';
  if FileExists(VenvPython) then
    WirePth(VenvPython, BinDir, MpDir)
  else
    LogToFile('  ERROR: ' + VenvPython + ' not produced after -m venv');
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  BinDir, MpDir, VenvDir, BasePython, VenvPython: String;
begin
  if CurStep = ssPostInstall then
  begin
    BinDir := ExpandConstant('{app}\bin');
    MpDir  := ExpandConstant('{app}\openseesmp');
    if not DirExists(MpDir) then
      MpDir := '';
    LogToFile('Ladruno OpenSees install log');
    LogToFile('  date       : ' + GetDateTimeString('yyyy-mm-dd hh:nn:ss', '-', ':'));
    LogToFile('  version    : ' + ExpandConstant('{#LadrunoVersion}'));
    LogToFile('  install dir: ' + ExpandConstant('{app}'));
    LogToFile('  bin dir    : ' + BinDir);
    LogToFile('  openseesmp : ' + MpDir);

    // /VENV= and /BASEPYTHON= win over the wizard; a silent run with neither
    // resolves to SKIP rather than guessing. See ResolveVenvChoice.
    case ResolveVenvChoice(VenvDir, BasePython) of
      0:
        if WizardSilent and (ParamOrEmpty('venv') = '') then
          LogToFile('  venv option: skip (silent install, no /VENV= given -- ' +
                    'pass /VENV=<path>, /VENV=new or /VENV=skip to be explicit)')
        else
          LogToFile('  venv option: skip');
      1:
        begin
          if BasePython = '' then BasePython := 'python';
          LogToFile('  venv option: create at ' + VenvDir);
          LogToFile('  base python: ' + BasePython);
          CreateAndWireVenv(BasePython, VenvDir, BinDir, MpDir);
        end;
      2:
        begin
          VenvPython := AddBackslash(VenvDir) + 'Scripts\python.exe';
          LogToFile('  venv option: existing at ' + VenvDir);
          if FileExists(VenvPython) then
            WirePth(VenvPython, BinDir, MpDir)
          else
            Complain('ERROR: ' + VenvPython + ' missing -- /VENV= wants the venv ' +
                     'ROOT (the folder containing Scripts\ and Lib\), not python.exe. ' +
                     'Nothing was wired.');
        end;
    end;

    if WizardIsTaskSelected('addtopath') then
      LogToFile('  PATH: appended ' + BinDir + ' to user PATH');
  end;
end;
