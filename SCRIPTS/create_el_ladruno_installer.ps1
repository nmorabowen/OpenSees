[CmdletBinding()]
param(
    [string]$BuildDir = "build-win11",
    [string]$OutputDir = "dist",
    [string]$AppName = "El Ladruno OpenSees",
    [string]$Publisher = "El Ladruno",
    [string]$AppVersion = "",
    [bool]$IncludeExamples = $true,
    [switch]$IncludeAllExamples,
    [switch]$SkipZip,
    [switch]$SkipInnoCompile,
    [string]$InnoCompilerPath = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Write-Log {
    param(
        [Parameter(Mandatory = $true)][string]$Message,
        [ValidateSet("INFO", "WARN", "ERROR")][string]$Level = "INFO"
    )
    $timestamp = (Get-Date).ToString("yyyy-MM-dd HH:mm:ss")
    $line = "[$timestamp] [$Level] $Message"
    switch ($Level) {
        "WARN" { Write-Host $line -ForegroundColor Yellow }
        "ERROR" { Write-Host $line -ForegroundColor Red }
        default { Write-Host $line }
    }
}

function Get-OpenSeesVersion {
    param(
        [Parameter(Mandatory = $true)][string]$OpenSeesExe,
        [Parameter(Mandatory = $true)][string]$RepoRoot
    )

    $cmakeLists = Join-Path $RepoRoot "CMakeLists.txt"
    if (Test-Path -Path $cmakeLists -PathType Leaf) {
        $text = Get-Content -Path $cmakeLists -Raw
        $m = [regex]::Match($text, 'project\s*\(\s*OpenSees\s+VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)', [System.Text.RegularExpressions.RegexOptions]::IgnoreCase)
        if ($m.Success) {
            return $m.Groups[1].Value
        }
    }

    $opsGlobals = Join-Path $RepoRoot "SRC\OPS_Globals.h"
    if (Test-Path -Path $opsGlobals -PathType Leaf) {
        $text = Get-Content -Path $opsGlobals -Raw
        $m = [regex]::Match($text, '#define\s+OPS_VERSION\s+"([0-9]+\.[0-9]+\.[0-9]+)"')
        if ($m.Success) {
            return $m.Groups[1].Value
        }
    }

    try {
        $output = cmd /c "echo exit | `"$OpenSeesExe`""
        $joined = ($output -join "`n")
        $m = [regex]::Match($joined, 'Version\s+([0-9]+\.[0-9]+\.[0-9]+)', [System.Text.RegularExpressions.RegexOptions]::IgnoreCase)
        if ($m.Success) {
            return $m.Groups[1].Value
        }
    } catch {
        Write-Log "Could not infer version from OpenSees executable: $($_.Exception.Message)" "WARN"
    }

    return (Get-Date).ToString("yyyy.MM.dd")
}

function Resolve-InnoCompiler {
    param([string]$Preferred)

    if ($Preferred -and (Test-Path -Path $Preferred -PathType Leaf)) {
        return (Resolve-Path $Preferred).Path
    }

    $known = @(
        (Join-Path $env:LOCALAPPDATA "Programs\Inno Setup 6\ISCC.exe"),
        "C:\Program Files (x86)\Inno Setup 6\ISCC.exe",
        "C:\Program Files\Inno Setup 6\ISCC.exe"
    )
    foreach ($candidate in $known) {
        if (Test-Path -Path $candidate -PathType Leaf) {
            return $candidate
        }
    }

    $command = Get-Command ISCC -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }

    return $null
}

function Copy-IfExists {
    param(
        [Parameter(Mandatory = $true)][string]$Source,
        [Parameter(Mandatory = $true)][string]$Destination
    )

    if (Test-Path -Path $Source) {
        Copy-Item -Path $Source -Destination $Destination -Force
    } else {
        Write-Log "Optional file not found and was skipped: $Source" "WARN"
    }
}

function Get-Vs2022InstallDir {
    $vswherePath = "C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe"
    if (Test-Path -Path $vswherePath -PathType Leaf) {
        $path = (& $vswherePath -latest -products * -property installationPath |
            Select-Object -First 1)
        if ($path) {
            return $path.Trim()
        }
    }

    $fallbacks = @(
        "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools",
        "C:\Program Files\Microsoft Visual Studio\2022\Community",
        "C:\Program Files\Microsoft Visual Studio\2022\Professional",
        "C:\Program Files\Microsoft Visual Studio\2022\Enterprise"
    )
    foreach ($candidate in $fallbacks) {
        if (Test-Path -Path $candidate -PathType Container) {
            return $candidate
        }
    }

    return $null
}

function New-InstalledLaunchers {
    param(
        [Parameter(Mandatory = $true)][string]$StageRoot
    )

    $tclDir = Get-ChildItem -Path (Join-Path $StageRoot "lib") -Directory -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -like "tcl*.*" } |
        Sort-Object Name |
        Select-Object -First 1

    $launcherHeader = @(
        "@echo off",
        "setlocal",
        "set `"OPS_ROOT=%~dp0`"",
        "set `"OPS_STATE=%LOCALAPPDATA%\ElLadrunoOpenSees`"",
        "if not exist `"%OPS_STATE%`" mkdir `"%OPS_STATE%`" >nul 2>&1"
    )
    $vsInstallDir = Get-Vs2022InstallDir
    if ($vsInstallDir) {
        $launcherHeader += "set `"VS2022INSTALLDIR=$vsInstallDir`""
    }
    $launcherHeader += "if exist `"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat`" call `"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat`" intel64 >nul 2>&1"
    if ($tclDir) {
        $launcherHeader += "set `"TCL_LIBRARY=%OPS_ROOT%lib\$($tclDir.Name)`""
    }

    $ladrunoHeader = @(
        "@echo off",
        "setlocal",
        "set `"OPS_ROOT=%~dp0`"",
        "set `"PATH=%OPS_ROOT%;%OPS_ROOT%opensees-bin;%OPS_ROOT%oneapi;%OPS_ROOT%plugins;%PATH%`""
    )
    if ($vsInstallDir) {
        $ladrunoHeader += "set `"VS2022INSTALLDIR=$vsInstallDir`""
    }
    $ladrunoHeader += @(
        "if exist `"%OPS_ROOT%oneapi\setvars.bat`" call `"%OPS_ROOT%oneapi\setvars.bat`" intel64 >nul 2>&1",
        "if exist `"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat`" call `"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat`" intel64 >nul 2>&1"
    )
    if ($tclDir) {
        $ladrunoHeader += "set `"TCL_LIBRARY=%OPS_ROOT%lib\$($tclDir.Name)`""
    }

    $serialLines = @($launcherHeader)
    $serialLines += @(
        "`"%OPS_ROOT%OpenSees.exe`" %*",
        "endlocal"
    )
    Set-Content -Path (Join-Path $StageRoot "OpenSees-Launch.cmd") -Value ($serialLines -join "`r`n") -Encoding Ascii

    $ladrunoSerialLines = @($ladrunoHeader)
    $ladrunoSerialLines += @(
        "set `"OPS_EXE=%OPS_ROOT%OpenSees.exe`"",
        "if not exist `"%OPS_EXE%`" set `"OPS_EXE=%OPS_ROOT%opensees-bin\OpenSees.exe`"",
        "if not exist `"%OPS_EXE%`" (",
        "  echo ERROR: OpenSees.exe not found under %OPS_ROOT%",
        "  endlocal & exit /b 1",
        ")",
        "if `"%~1`"==`"`" (",
        "  `"%OPS_EXE%`"",
        ") else (",
        "  pushd `"%~dp1`"",
        "  `"%OPS_EXE%`" `"%~nx1`"",
        "  set `"OPS_EXIT=%ERRORLEVEL%`"",
        "  popd",
        "  endlocal & exit /b %OPS_EXIT%",
        ")",
        "endlocal"
    )
    Set-Content -Path (Join-Path $StageRoot "opensees_ladruno.bat") -Value ($ladrunoSerialLines -join "`r`n") -Encoding Ascii

    $spLines = @($launcherHeader)
    $spLines += @(
        "if `"%OPS_MPI_N%`"==`"`" set `"OPS_MPI_N=2`"",
        "set `"OPS_LOG=%OPS_STATE%\OpenSeesSP-launch.log`"",
        "where mpiexec >nul 2>&1",
        "if errorlevel 1 (",
        "  echo ERROR: mpiexec not found. Intel MPI must be installed and initialized.",
        "  pause",
        "  endlocal & exit /b 1",
        ")",
        "if not `"%~1`"==`"`" goto run_args",
        "if not exist `"%OPS_ROOT%EXAMPLES\ParallelModelMP`" goto no_default",
        "pushd `"%OPS_ROOT%EXAMPLES\ParallelModelMP`"",
        "mpiexec -localonly -n %OPS_MPI_N% `"%OPS_ROOT%OpenSeesSP.exe`" exampleSP.tcl >> `"%OPS_LOG%`" 2>&1",
        "set `"OPS_EXIT=%ERRORLEVEL%`"",
        "popd",
        "goto finalize",
        ":run_args",
        "mpiexec -localonly -n %OPS_MPI_N% `"%OPS_ROOT%OpenSeesSP.exe`" %* >> `"%OPS_LOG%`" 2>&1",
        "set `"OPS_EXIT=%ERRORLEVEL%`"",
        "goto finalize",
        ":no_default",
        "echo ERROR: default SP example folder not found: %OPS_ROOT%EXAMPLES\ParallelModelMP",
        "echo Usage: OpenSeesSP-Launch.cmd your_parallel_script.tcl",
        "set `"OPS_EXIT=1`"",
        ":finalize",
        "if `"%OPS_EXIT%`"==`"0`" (",
        "  endlocal & exit /b 0",
        ")",
        "echo ERROR: OpenSeesSP MPI launch failed. See %OPS_LOG%",
        "echo Hint: open a oneAPI shell or run: call `"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat`" intel64",
        "pause",
        "endlocal & exit /b %OPS_EXIT%"
    )
    Set-Content -Path (Join-Path $StageRoot "OpenSeesSP-Launch.cmd") -Value ($spLines -join "`r`n") -Encoding Ascii

    $ladrunoSpLines = @($ladrunoHeader)
    $ladrunoSpLines += @(
        "set `"OPS_EXE=%OPS_ROOT%OpenSeesSP.exe`"",
        "if not exist `"%OPS_EXE%`" set `"OPS_EXE=%OPS_ROOT%opensees-bin\OpenSeesSP.exe`"",
        "if not exist `"%OPS_EXE%`" (",
        "  echo ERROR: OpenSeesSP.exe not found under %OPS_ROOT%",
        "  endlocal & exit /b 1",
        ")",
        "if `"%~2`"==`"`" (",
        "  echo RUNNING AS SEQUENTIAL",
        "  if `"%~1`"==`"`" (",
        "    `"%OPS_EXE%`"",
        "  ) else (",
        "    pushd `"%~dp1`"",
        "    `"%OPS_EXE%`" `"%~nx1`"",
        "    set `"OPS_EXIT=%ERRORLEVEL%`"",
        "    popd",
        "    endlocal & exit /b %OPS_EXIT%",
        "  )",
        "  endlocal & exit /b %ERRORLEVEL%",
        ")",
        "echo RUNNING AS PARALLEL",
        "where mpiexec >nul 2>&1",
        "if errorlevel 1 (",
        "  echo ERROR: mpiexec not found. Initialize Intel MPI first.",
        "  endlocal & exit /b 1",
        ")",
        "if `"%~1`"==`"`" (",
        "  mpiexec -n %2 `"%OPS_EXE%`"",
        ") else (",
        "  pushd `"%~dp1`"",
        "  mpiexec -n %2 `"%OPS_EXE%`" `"%~nx1`"",
        "  set `"OPS_EXIT=%ERRORLEVEL%`"",
        "  popd",
        "  endlocal & exit /b %OPS_EXIT%",
        ")",
        "endlocal"
    )
    Set-Content -Path (Join-Path $StageRoot "opensees_sp_ladruno.bat") -Value ($ladrunoSpLines -join "`r`n") -Encoding Ascii

    $mpLines = @($launcherHeader)
    $mpLines += @(
        "if `"%OPS_MPI_N%`"==`"`" set `"OPS_MPI_N=2`"",
        "set `"OPS_LOG=%OPS_STATE%\OpenSeesMP-launch.log`"",
        "where mpiexec >nul 2>&1",
        "if errorlevel 1 (",
        "  echo ERROR: mpiexec not found. Intel MPI must be installed and initialized.",
        "  pause",
        "  endlocal & exit /b 1",
        ")",
        "if not `"%~1`"==`"`" goto run_args",
        "if not exist `"%OPS_ROOT%EXAMPLES\ParallelModelMP`" goto no_default",
        "pushd `"%OPS_ROOT%EXAMPLES\ParallelModelMP`"",
        "mpiexec -localonly -n %OPS_MPI_N% `"%OPS_ROOT%OpenSeesMP.exe`" exampleMP.tcl >> `"%OPS_LOG%`" 2>&1",
        "set `"OPS_EXIT=%ERRORLEVEL%`"",
        "popd",
        "goto finalize",
        ":run_args",
        "mpiexec -localonly -n %OPS_MPI_N% `"%OPS_ROOT%OpenSeesMP.exe`" %* >> `"%OPS_LOG%`" 2>&1",
        "set `"OPS_EXIT=%ERRORLEVEL%`"",
        "goto finalize",
        ":no_default",
        "echo ERROR: default MP example folder not found: %OPS_ROOT%EXAMPLES\ParallelModelMP",
        "echo Usage: OpenSeesMP-Launch.cmd your_parallel_script.tcl",
        "set `"OPS_EXIT=1`"",
        ":finalize",
        "if `"%OPS_EXIT%`"==`"0`" (",
        "  endlocal & exit /b 0",
        ")",
        "echo ERROR: OpenSeesMP MPI launch failed. See %OPS_LOG%",
        "echo Hint: open a oneAPI shell or run: call `"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat`" intel64",
        "pause",
        "endlocal & exit /b %OPS_EXIT%"
    )
    Set-Content -Path (Join-Path $StageRoot "OpenSeesMP-Launch.cmd") -Value ($mpLines -join "`r`n") -Encoding Ascii

    $ladrunoMpLines = @($ladrunoHeader)
    $ladrunoMpLines += @(
        "set `"OPS_EXE=%OPS_ROOT%OpenSeesMP.exe`"",
        "if not exist `"%OPS_EXE%`" set `"OPS_EXE=%OPS_ROOT%opensees-bin\OpenSeesMP.exe`"",
        "if not exist `"%OPS_EXE%`" (",
        "  echo ERROR: OpenSeesMP.exe not found under %OPS_ROOT%",
        "  endlocal & exit /b 1",
        ")",
        "if `"%~2`"==`"`" (",
        "  echo RUNNING AS SEQUENTIAL",
        "  if `"%~1`"==`"`" (",
        "    `"%OPS_EXE%`"",
        "  ) else (",
        "    pushd `"%~dp1`"",
        "    `"%OPS_EXE%`" `"%~nx1`"",
        "    set `"OPS_EXIT=%ERRORLEVEL%`"",
        "    popd",
        "    endlocal & exit /b %OPS_EXIT%",
        "  )",
        "  endlocal & exit /b %ERRORLEVEL%",
        ")",
        "echo RUNNING AS PARALLEL",
        "where mpiexec >nul 2>&1",
        "if errorlevel 1 (",
        "  echo ERROR: mpiexec not found. Initialize Intel MPI first.",
        "  endlocal & exit /b 1",
        ")",
        "if `"%~1`"==`"`" (",
        "  mpiexec -n %2 `"%OPS_EXE%`"",
        ") else (",
        "  pushd `"%~dp1`"",
        "  mpiexec -n %2 `"%OPS_EXE%`" `"%~nx1`"",
        "  set `"OPS_EXIT=%ERRORLEVEL%`"",
        "  popd",
        "  endlocal & exit /b %OPS_EXIT%",
        ")",
        "endlocal"
    )
    Set-Content -Path (Join-Path $StageRoot "opensees_mp_ladruno.bat") -Value ($ladrunoMpLines -join "`r`n") -Encoding Ascii
}

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$buildRoot = if ([System.IO.Path]::IsPathRooted($BuildDir)) { $BuildDir } else { Join-Path $repoRoot $BuildDir }
$outputRoot = if ([System.IO.Path]::IsPathRooted($OutputDir)) { $OutputDir } else { Join-Path $repoRoot $OutputDir }

if (-not (Test-Path -Path $buildRoot -PathType Container)) {
    throw "Build directory does not exist: $buildRoot"
}

$openSeesExe = Join-Path $buildRoot "OpenSees.exe"
if (-not (Test-Path -Path $openSeesExe -PathType Leaf)) {
    throw "Required executable not found: $openSeesExe"
}

if (-not $AppVersion) {
    $AppVersion = Get-OpenSeesVersion -OpenSeesExe $openSeesExe -RepoRoot $repoRoot
}

$safeVersion = $AppVersion.Replace(" ", "_")
$stageRoot = Join-Path $outputRoot "ElLadrunoOpenSees-stage"
$issPath = Join-Path $outputRoot "ElLadrunoOpenSees.iss"
$portableZip = Join-Path $outputRoot ("ElLadrunoOpenSees_{0}_portable.zip" -f $safeVersion)

Write-Log "RepoRoot: $repoRoot"
Write-Log "BuildRoot: $buildRoot"
Write-Log "OutputRoot: $outputRoot"
Write-Log "AppName: $AppName"
Write-Log "AppVersion: $AppVersion"

New-Item -ItemType Directory -Path $outputRoot -Force | Out-Null
if (Test-Path -Path $stageRoot) {
    Remove-Item -Path $stageRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $stageRoot -Force | Out-Null

Write-Log "Copying core artifacts"
$requiredFiles = @(
    "OpenSees.exe",
    "OpenSeesSP.exe",
    "OpenSeesMP.exe",
    "OpenSeesPy.dll",
    "opensees.pyd",
    "OpenSees-Launch.cmd",
    "OpenSeesSP-Launch.cmd",
    "OpenSeesMP-Launch.cmd"
)

foreach ($file in $requiredFiles) {
    $src = Join-Path $buildRoot $file
    if (-not (Test-Path -Path $src -PathType Leaf)) {
        throw "Missing required artifact: $src"
    }
    Copy-Item -Path $src -Destination $stageRoot -Force
}

Write-Log "Copying runtime DLLs"
$dlls = Get-ChildItem -Path $buildRoot -Filter *.dll -File -ErrorAction SilentlyContinue
foreach ($dll in $dlls) {
    Copy-Item -Path $dll.FullName -Destination $stageRoot -Force
}

Write-Log "Copying Tcl runtime"
$libSource = Join-Path $buildRoot "lib"
if (-not (Test-Path -Path $libSource -PathType Container)) {
    throw "Missing Tcl runtime folder: $libSource"
}
Copy-Item -Path $libSource -Destination $stageRoot -Recurse -Force

Write-Log "Copying metadata files"
$versionTxt = @(
    "Product: $AppName",
    "Version: $AppVersion",
    "BuiltFrom: $repoRoot",
    "BuildDir: $buildRoot",
    "PackagedAt: $((Get-Date).ToString('yyyy-MM-dd HH:mm:ss'))"
) -join "`r`n"
Set-Content -Path (Join-Path $stageRoot "VERSION.txt") -Value $versionTxt -Encoding Ascii

$readmeTxt = @'
El Ladruno OpenSees
===================

Quick start:
1) Run OpenSees-Launch.cmd for interactive serial OpenSees.
2) Run OpenSeesSP-Launch.cmd or OpenSeesMP-Launch.cmd for MPI examples.

Notes:
- OpenSeesSP/OpenSeesMP require Intel MPI runtime on the target machine.
- If MPI fails, test with: mpiexec -n 2 hostname
'@
Set-Content -Path (Join-Path $stageRoot "README_ElLadruno.txt") -Value $readmeTxt -Encoding Ascii

if ($IncludeAllExamples -or $IncludeExamples) {
    $examplesRoot = Join-Path $repoRoot "EXAMPLES"
    if (Test-Path -Path $examplesRoot -PathType Container) {
        Write-Log "Copying examples"
        $destExamples = Join-Path $stageRoot "EXAMPLES"
        New-Item -ItemType Directory -Path $destExamples -Force | Out-Null
        if ($IncludeAllExamples) {
            Copy-Item -Path (Join-Path $examplesRoot "*") -Destination $destExamples -Recurse -Force
        } else {
            foreach ($name in @("ParallelModelMP", "SmallSP", "SmallMP", "ExamplePython")) {
                $src = Join-Path $examplesRoot $name
                if (Test-Path -Path $src -PathType Container) {
                    Copy-Item -Path $src -Destination $destExamples -Recurse -Force
                } else {
                    Write-Log "Example folder not found and skipped: $src" "WARN"
                }
            }
        }
    } else {
        Write-Log "EXAMPLES folder not found; skipping examples." "WARN"
    }
}

Write-Log "Rewriting packaged launchers for installed paths"
New-InstalledLaunchers -StageRoot $stageRoot

Write-Log "Generating Inno Setup script"
$myBuildRoot = $stageRoot
$issText = @"
#define MyAppName "$AppName"
#define MyAppVersion "$AppVersion"
#define MyAppPublisher "$Publisher"
#define MyAppExeName "OpenSees-Launch.cmd"
#define MyBuildRoot "$myBuildRoot"

[Setup]
AppId={{A1F0FA32-B138-49A5-9F34-39D23D69BB6F}}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={autopf}\El Ladruno OpenSees
DefaultGroupName=El Ladruno OpenSees
DisableProgramGroupPage=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir=$outputRoot
OutputBaseFilename=ElLadrunoOpenSeesSetup_{#MyAppVersion}_x64
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=admin
ChangesEnvironment=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional icons:"

[Files]
Source: "{#MyBuildRoot}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\El Ladruno OpenSees"; Filename: "{app}\OpenSees-Launch.cmd"
Name: "{group}\El Ladruno OpenSees SP"; Filename: "{app}\OpenSeesSP-Launch.cmd"
Name: "{group}\El Ladruno OpenSees MP"; Filename: "{app}\OpenSeesMP-Launch.cmd"
Name: "{autodesktop}\El Ladruno OpenSees"; Filename: "{app}\OpenSees-Launch.cmd"; Tasks: desktopicon

[Run]
Filename: "{app}\OpenSees-Launch.cmd"; Description: "Launch El Ladruno OpenSees"; Flags: nowait postinstall skipifsilent
"@
Set-Content -Path $issPath -Value $issText -Encoding Ascii

if (-not $SkipZip) {
    Write-Log "Creating portable ZIP: $portableZip"
    if (Test-Path -Path $portableZip) {
        Remove-Item -Path $portableZip -Force
    }
    Compress-Archive -Path (Join-Path $stageRoot "*") -DestinationPath $portableZip -CompressionLevel Optimal
}

if (-not $SkipInnoCompile) {
    $iscc = Resolve-InnoCompiler -Preferred $InnoCompilerPath
    if ($iscc) {
        Write-Log "Compiling installer with Inno Setup: $iscc"
        & $iscc $issPath
        if ($LASTEXITCODE -ne 0) {
            throw "Inno Setup compilation failed with exit code $LASTEXITCODE"
        }
        Write-Log "Installer compilation completed."
    } else {
        Write-Log "Inno Setup compiler (ISCC.exe) not found. Installer .iss generated only." "WARN"
    }
}

Write-Log "Staging folder: $stageRoot"
Write-Log "Inno script: $issPath"
if (-not $SkipZip) {
    Write-Log "Portable zip: $portableZip"
}
