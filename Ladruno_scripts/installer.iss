; Ladruno OpenSees - Inno Setup script
;
; Build:
;   $env:LADRUNO_VERSION = '20260501'      ; or any string
;   $env:LADRUNO_DIST    = 'C:\...\dist'   ; abs path to the dist folder
;   $env:LADRUNO_OUT     = 'C:\...\Ladruno_files'
;   & 'C:\Program Files (x86)\Inno Setup 6\iscc.exe' installer.iss
;
; The build_inno_installer.ps1 driver does this for you.

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
    '          opensees.pyd (Python 3.12), HDF5 + MPCO recorder,' + #13#10 +
    '          MUMPS sparse solver.' + #13#10 + #13#10 +
    'OpenSeesSP / OpenSeesMP need Intel MPI''s mpiexec on PATH' + #13#10 +
    '(Intel oneAPI install required separately).';

  // ---- venv option page (after wpSelectDir so {app} is final) ----------
  VenvOptionPage := CreateInputOptionPage(wpSelectDir,
    'Wire OpenSeesPy into a Python virtualenv',
    'Make ''import opensees'' work without sys.path fiddling',
    'opensees.pyd is built for CPython 3.12. Choose what to do:',
    True, False);
  VenvOptionPage.Add('Skip - I''ll wire it up later');
  VenvOptionPage.Add('Create a new venv at <install-dir>\opensees_venv');
  VenvOptionPage.Add('Use an existing venv (choose path on the next page)');
  VenvOptionPage.SelectedValueIndex := 1;

  // ---- existing venv path page -----------------------------------------
  VenvDirPage := CreateInputDirPage(VenvOptionPage.ID,
    'Existing virtualenv path',
    'Choose the venv root folder (NOT python.exe)',
    'Point at the directory that contains "Scripts\python.exe" and "Lib\site-packages\" - e.g. C:\envs\myproject - not at python.exe itself. Setup writes a .pth file inside that venv so "import opensees" finds the binaries.',
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

procedure WirePth(const VenvPython, BinDir: String);
var
  HelperPy: String;
  Params:   String;
  Code:     Integer;
begin
  ExtractTemporaryFile('wire_venv_pth.py');
  HelperPy := ExpandConstant('{tmp}\wire_venv_pth.py');
  Params := AddQuotes(HelperPy) + ' ' + AddQuotes(BinDir);
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

procedure CreateAndWireVenv(const BasePython, VenvDir, BinDir: String);
var
  Code: Integer;
  VenvPython: String;
begin
  LogToFile('Creating venv at ' + VenvDir + ' (base: ' + BasePython + ') ...');
  if not Exec(BasePython, '-m venv ' + AddQuotes(VenvDir), '', SW_HIDE, ewWaitUntilTerminated, Code) then
  begin
    LogToFile('  ERROR: failed to launch ' + BasePython);
    MsgBox('Failed to launch ' + BasePython + '.' + #13#10 +
           'Skipping venv wiring; you can run Ladruno_scripts\wire_pyenv.ps1 later.',
           mbError, MB_OK);
    Exit;
  end;
  if Code <> 0 then
  begin
    LogToFile('  ERROR: ' + BasePython + ' -m venv exited ' + IntToStr(Code));
    MsgBox(BasePython + ' -m venv exited ' + IntToStr(Code) + '.' + #13#10 +
           'Skipping venv wiring.',
           mbError, MB_OK);
    Exit;
  end;
  VenvPython := AddBackslash(VenvDir) + 'Scripts\python.exe';
  if FileExists(VenvPython) then
    WirePth(VenvPython, BinDir)
  else
    LogToFile('  ERROR: ' + VenvPython + ' not produced after -m venv');
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  BinDir, VenvDir, BasePython, VenvPython: String;
begin
  if CurStep = ssPostInstall then
  begin
    BinDir := ExpandConstant('{app}\bin');
    LogToFile('Ladruno OpenSees install log');
    LogToFile('  date       : ' + GetDateTimeString('yyyy-mm-dd hh:nn:ss', '-', ':'));
    LogToFile('  version    : ' + ExpandConstant('{#LadrunoVersion}'));
    LogToFile('  install dir: ' + ExpandConstant('{app}'));
    LogToFile('  bin dir    : ' + BinDir);

    case VenvOptionPage.SelectedValueIndex of
      0:
        LogToFile('  venv option: skip');
      1:
        begin
          VenvDir    := ExpandConstant('{app}\opensees_venv');
          BasePython := Trim(PythonPathPage.Values[0]);
          if BasePython = '' then BasePython := 'python';
          LogToFile('  venv option: create at ' + VenvDir);
          LogToFile('  base python: ' + BasePython);
          CreateAndWireVenv(BasePython, VenvDir, BinDir);
        end;
      2:
        begin
          VenvDir    := VenvDirPage.Values[0];
          VenvPython := AddBackslash(VenvDir) + 'Scripts\python.exe';
          LogToFile('  venv option: existing at ' + VenvDir);
          if FileExists(VenvPython) then
            WirePth(VenvPython, BinDir)
          else
            LogToFile('  ERROR: ' + VenvPython + ' missing');
        end;
    end;

    if WizardIsTaskSelected('addtopath') then
      LogToFile('  PATH: appended ' + BinDir + ' to user PATH');
  end;
end;
