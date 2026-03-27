[CmdletBinding()]
param(
    [string]$BuildDir = "build-win11",
    [string]$Triplet = "x64-windows-static",
    [string]$VcpkgRoot = "",
    [string]$MumpsRoot = "",
    [ValidateSet("quick", "full")][string]$SmokeMode = "quick",
    [int]$SmokeTimeoutSec = 600,
    [int]$Parallel = 0,
    [switch]$SkipMumps,
    [switch]$SkipBuild,
    [switch]$SkipTests,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$script:LogFile = $null
$script:DryRunMode = $DryRun.IsPresent

# Default parallel jobs to CPU count when not specified.
if ($Parallel -le 0) {
    $Parallel = [Environment]::ProcessorCount
}

# Pin MUMPS to a known-good commit for reproducibility.
# Pin MUMPS to the latest known-good commit from github.com/OpenSees/mumps.
# Update this hash when upgrading MUMPS.
$MumpsCommit = "ec5f340"

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

    if ($script:LogFile) {
        Add-Content -Path $script:LogFile -Value $line
    }
}

function Write-Step {
    param([string]$Message)
    Write-Host ""
    if ($script:LogFile) {
        Add-Content -Path $script:LogFile -Value ""
    }
    Write-Log "==> $Message"
}

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)][string]$FilePath,
        [Parameter()][string[]]$Arguments = @(),
        [Parameter()][int]$TimeoutSec = 0,
        [Parameter()][string]$WorkingDirectory = ""
    )

    $argText = ($Arguments -join " ")
    Write-Log ">> $FilePath $argText"

    if ($script:DryRunMode) {
        Write-Log "[DRY RUN] Would execute: $FilePath $argText"
        return
    }

    if ($TimeoutSec -gt 0) {
        $startParams = @{
            FilePath    = $FilePath
            ArgumentList = $Arguments
            NoNewWindow = $true
            PassThru    = $true
        }
        if ($WorkingDirectory) {
            $startParams["WorkingDirectory"] = $WorkingDirectory
        }

        $proc = Start-Process @startParams
        $timedOut = -not $proc.WaitForExit($TimeoutSec * 1000)
        if ($timedOut) {
            Write-Log "Command timed out after $TimeoutSec seconds: $FilePath" "ERROR"
            Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
            Start-Sleep -Milliseconds 200
            Stop-StaleOpenSeesProcesses
            throw "Command timed out after $TimeoutSec seconds: $FilePath"
        }

        $exitCode = $proc.ExitCode
        if ($null -eq $exitCode) {
            $proc.Refresh()
            $exitCode = $proc.ExitCode
        }
        if ($null -eq $exitCode) {
            $exitCode = 1
        }

        if ($exitCode -ne 0) {
            throw "Command failed with exit code ${exitCode}: $FilePath"
        }
        return
    }

    if ($WorkingDirectory) {
        Push-Location $WorkingDirectory
        try {
            & $FilePath @Arguments
            $exitCode = $LASTEXITCODE
        } finally {
            Pop-Location
        }
    } else {
        & $FilePath @Arguments
        $exitCode = $LASTEXITCODE
    }

    if ($exitCode -ne 0) {
        throw "Command failed with exit code ${exitCode}: $FilePath"
    }
}

function Invoke-MpiSmoke {
    param(
        [Parameter(Mandatory = $true)][string]$FilePath,
        [Parameter()][string[]]$Arguments = @(),
        [Parameter()][int]$TimeoutSec = 0,
        [Parameter()][string]$WorkingDirectory = ""
    )

    $argText = ($Arguments -join " ")
    Write-Log ">> $FilePath $argText"

    $outFile = Join-Path $env:TEMP ("opensees_mpi_out_" + [guid]::NewGuid().ToString("N") + ".log")
    $errFile = Join-Path $env:TEMP ("opensees_mpi_err_" + [guid]::NewGuid().ToString("N") + ".log")

    try {
        $startParams = @{
            FilePath = $FilePath
            ArgumentList = $Arguments
            NoNewWindow = $true
            PassThru = $true
            RedirectStandardOutput = $outFile
            RedirectStandardError = $errFile
        }
        if ($WorkingDirectory) {
            $startParams["WorkingDirectory"] = $WorkingDirectory
        }

        $proc = Start-Process @startParams

        if ($TimeoutSec -gt 0) {
            $timedOut = -not $proc.WaitForExit($TimeoutSec * 1000)
            if ($timedOut) {
                Write-Log "MPI smoke command timed out after $TimeoutSec seconds." "ERROR"
                Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
                Start-Sleep -Milliseconds 200
                Stop-StaleOpenSeesProcesses
                throw "MPI smoke command timed out after $TimeoutSec seconds."
            }
        } else {
            $proc.WaitForExit() | Out-Null
        }

        $combined = @()
        if (Test-Path -Path $outFile) {
            $combined += Get-Content -Path $outFile -ErrorAction SilentlyContinue
        }
        if (Test-Path -Path $errFile) {
            $combined += Get-Content -Path $errFile -ErrorAction SilentlyContinue
        }

        $combinedText = ($combined -join "`n")
        $fatalPatterns = @(
            "BAD TERMINATION",
            "application-specific initialization failed",
            "Fatal error in PMPI_Init",
            "Unknown error class",
            "Access violation",
            "The memory could not be read",
            "HYD_spawn",
            "unable to create stdout pipe",
            "Intel oneMKL ERROR",
            "ERROR --"
        )

        foreach ($pattern in $fatalPatterns) {
            if ($combinedText -match [regex]::Escape($pattern)) {
                $tail = ($combined | Select-Object -Last 80) -join "`n"
                if ($tail) {
                    Write-Log "MPI output tail:`n$tail" "ERROR"
                }
                throw "MPI smoke detected fatal pattern: '$pattern'"
            }
        }

        $exitCode = $proc.ExitCode
        if ($null -eq $exitCode) {
            $proc.Refresh()
            $exitCode = $proc.ExitCode
        }
        if ($null -eq $exitCode) {
            $exitCode = 1
        }

        if ($exitCode -ne 0) {
            $rank0Terminated = $combinedText -match "Process Terminating 0"
            $rank1Terminated = $combinedText -match "Process Terminating 1"
            if ($rank0Terminated -and $rank1Terminated) {
                Write-Log "MPI smoke returned non-zero exit code (${exitCode}) but both ranks reported clean termination. Continuing." "WARN"
                return
            }

            $tail = ($combined | Select-Object -Last 80) -join "`n"
            if ($tail) {
                Write-Log "MPI output tail:`n$tail" "ERROR"
            }
            throw "MPI smoke command failed with exit code ${exitCode}."
        }
    } finally {
        Remove-Item -Path $outFile, $errFile -Force -ErrorAction SilentlyContinue
    }
}

function Import-BatchEnvironment {
    param(
        [Parameter(Mandatory = $true)][string]$BatchPath,
        [string]$Arguments = ""
    )

    if (-not (Test-Path -Path $BatchPath -PathType Leaf)) {
        throw "Batch file not found: $BatchPath"
    }

    $cmdLine = if ([string]::IsNullOrWhiteSpace($Arguments)) {
        "`"$BatchPath`" && set"
    } else {
        "`"$BatchPath`" $Arguments && set"
    }

    $output = & cmd.exe /d /s /c $cmdLine
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to load environment from $BatchPath"
    }

    foreach ($line in $output) {
        $idx = $line.IndexOf("=")
        if ($idx -gt 0) {
            $name = $line.Substring(0, $idx)
            $value = $line.Substring($idx + 1)
            Set-Item -Path "Env:$name" -Value $value
        }
    }
}

function Get-LatestArtifact {
    param(
        [Parameter(Mandatory = $true)][string]$RootDir,
        [Parameter(Mandatory = $true)][string]$FileName
    )

    $item = Get-ChildItem -Path $RootDir -Recurse -File -Filter $FileName -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
    if (-not $item) {
        throw "Could not find '$FileName' under '$RootDir'"
    }
    return $item
}

function Get-DumpbinDependencies {
    param([Parameter(Mandatory = $true)][string]$BinaryPath)

    $dumpbin = Get-Command dumpbin.exe -ErrorAction SilentlyContinue
    if (-not $dumpbin) {
        return @()
    }

    $lines = & $dumpbin.Source /dependents $BinaryPath 2>$null
    if ($LASTEXITCODE -ne 0) {
        return @()
    }

    $deps = @()
    foreach ($line in $lines) {
        if ($line -match "^\s+([A-Za-z0-9_.-]+\.dll)\s*$") {
            $deps += $Matches[1]
        }
    }
    return $deps | Sort-Object -Unique
}

function Resolve-ToolPath {
    param(
        [Parameter(Mandatory = $true)][string]$Name,
        [string[]]$Candidates = @()
    )

    $cmd = Get-Command $Name -ErrorAction SilentlyContinue
    if ($cmd -and $cmd.Source) {
        return $cmd.Source
    }

    foreach ($candidate in $Candidates) {
        if ($candidate -and (Test-Path -Path $candidate -PathType Leaf)) {
            return (Resolve-Path $candidate).Path
        }
    }

    throw "Required tool '$Name' not found."
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

function Stop-StaleOpenSeesProcesses {
    $names = @("OpenSees", "OpenSeesSP", "OpenSeesMP", "mpiexec")
    foreach ($name in $names) {
        Get-Process -Name $name -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
    }
}

function New-LadrunoWrapperLines {
    param(
        [Parameter(Mandatory = $true)][ValidateSet("serial", "sp", "mp")][string]$Mode
    )

    $lines = @(
        "@echo off",
        "setlocal",
        "set `"OPS_ROOT=%~dp0`"",
        "set `"PATH=%OPS_ROOT%;%OPS_ROOT%opensees-bin;%OPS_ROOT%oneapi;%OPS_ROOT%plugins;%PATH%`"",
        "if defined VS2022INSTALLDIR set `"VS2022INSTALLDIR=%VS2022INSTALLDIR%`"",
        "if exist `"%OPS_ROOT%oneapi\setvars.bat`" call `"%OPS_ROOT%oneapi\setvars.bat`" intel64 >nul 2>&1",
        "if exist `"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat`" call `"%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat`" intel64 >nul 2>&1",
        "if exist `"%OPS_ROOT%lib\tcl8.6\init.tcl`" set `"TCL_LIBRARY=%OPS_ROOT%lib\tcl8.6`"",
        "if not defined TCL_LIBRARY if exist `"%OPS_ROOT%lib\tcl9.0\init.tcl`" set `"TCL_LIBRARY=%OPS_ROOT%lib\tcl9.0`""
    )

    switch ($Mode) {
        "serial" {
            $lines += @(
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
        }
        "sp" {
            $lines += @(
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
        }
        "mp" {
            $lines += @(
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
        }
    }

    return $lines
}

$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $RepoRoot

if (-not [System.IO.Path]::IsPathRooted($BuildDir)) {
    $BuildDir = Join-Path $RepoRoot $BuildDir
}
if ($SmokeTimeoutSec -le 0) {
    throw "SmokeTimeoutSec must be greater than zero."
}

if ([string]::IsNullOrWhiteSpace($VcpkgRoot)) {
    if ($env:VCPKG_ROOT) {
        $VcpkgRoot = $env:VCPKG_ROOT
    } else {
        $VcpkgRoot = Join-Path (Split-Path $RepoRoot -Parent) "vcpkg"
    }
}
if (-not [System.IO.Path]::IsPathRooted($VcpkgRoot)) {
    $VcpkgRoot = Join-Path $RepoRoot $VcpkgRoot
}
$env:VCPKG_ROOT = $VcpkgRoot

if ([string]::IsNullOrWhiteSpace($MumpsRoot)) {
    $MumpsRoot = Join-Path (Split-Path $RepoRoot -Parent) "mumps"
}
if (-not [System.IO.Path]::IsPathRooted($MumpsRoot)) {
    $MumpsRoot = Join-Path $RepoRoot $MumpsRoot
}

New-Item -ItemType Directory -Path $BuildDir -Force | Out-Null
$script:LogFile = Join-Path $BuildDir "build_windows11_full.log"
Set-Content -Path $script:LogFile -Value ("[{0}] [INFO] build_windows11_full.ps1 started" -f (Get-Date).ToString("yyyy-MM-dd HH:mm:ss"))
Write-Log "RepoRoot: $RepoRoot"
Write-Log "BuildDir: $BuildDir"
Write-Log "Options: Triplet=$Triplet SkipMumps=$SkipMumps SkipBuild=$SkipBuild SkipTests=$SkipTests SmokeMode=$SmokeMode SmokeTimeoutSec=$SmokeTimeoutSec"

$vcpkgInstallRoot = Join-Path $RepoRoot "vcpkg_installed"
$tclRuntimeName = "tcl8.6"
$tclRuntimeDir = $null
$tclPatchLevel = $null

try {
    Write-Step "Cleaning stale OpenSees/MPI processes"
    Stop-StaleOpenSeesProcesses

    Write-Step "Checking required command-line tools"
    $gitCmd = Resolve-ToolPath -Name "git" -Candidates @()
    $cmakeCmd = Resolve-ToolPath -Name "cmake" -Candidates @(
        "C:\Program Files\CMake\bin\cmake.exe",
        (Join-Path $env:ProgramFiles "CMake\bin\cmake.exe")
    )
    $python311Candidates = @(
        (Join-Path $env:LOCALAPPDATA "Programs\Python\Python311\python.exe"),
        "C:\Python311\python.exe"
    )
    $pythonCmd = $null
    foreach ($candidate in $python311Candidates) {
        if (Test-Path -Path $candidate -PathType Leaf) {
            $pythonCmd = (Resolve-Path $candidate).Path
            break
        }
    }
    if (-not $pythonCmd) {
        $pythonCmd = Resolve-ToolPath -Name "python" -Candidates @()
    }
    $pythonVersion = (& $pythonCmd -c "import sys; print(f'{sys.version_info[0]}.{sys.version_info[1]}')").Trim()
    if ($pythonVersion -ne "3.11") {
        throw "Python 3.11 is required for OpenSeesPy. Resolved: '$pythonCmd' (version $pythonVersion)."
    }

    if ($env:PATH -notlike "*$([IO.Path]::GetDirectoryName($cmakeCmd))*") {
        $env:PATH = "$([IO.Path]::GetDirectoryName($cmakeCmd));$env:PATH"
    }
    if ($env:PATH -notlike "*$([IO.Path]::GetDirectoryName($pythonCmd))*") {
        $env:PATH = "$([IO.Path]::GetDirectoryName($pythonCmd));$env:PATH"
    }
    if ((-not $SkipMumps) -and (-not $SkipBuild)) {
        $ninjaCmd = Resolve-ToolPath -Name "ninja" -Candidates @(
            (Join-Path $env:LOCALAPPDATA "Microsoft\WinGet\Packages\Ninja-build.Ninja_Microsoft.Winget.Source_8wekyb3d8bbwe\ninja.exe"),
            "C:\Program Files\Ninja\ninja.exe"
        )
        if ($env:PATH -notlike "*$([IO.Path]::GetDirectoryName($ninjaCmd))*") {
            $env:PATH = "$([IO.Path]::GetDirectoryName($ninjaCmd));$env:PATH"
        }
    }

    Write-Step "Loading Visual Studio Build Tools environment"
    $vsInstallPath = Get-Vs2022InstallDir
    if (-not $vsInstallPath) {
        $vsInstallPath = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
    }
    $vsInstallPath = $vsInstallPath.Trim()
    $vsDevCmdPath = Join-Path $vsInstallPath "Common7\Tools\VsDevCmd.bat"
    Import-BatchEnvironment -BatchPath $vsDevCmdPath -Arguments "-arch=x64 -host_arch=x64"

    Write-Step "Loading Intel oneAPI environment"
    # oneAPI relies on VS2022INSTALLDIR on some systems even when VsDevCmd was
    # already loaded in the current shell.
    $env:VS2022INSTALLDIR = $vsInstallPath
    $programFilesX86 = [Environment]::GetFolderPath("ProgramFilesX86")
    $setvarsPath = Join-Path $programFilesX86 "Intel\oneAPI\setvars.bat"
    Import-BatchEnvironment -BatchPath $setvarsPath -Arguments "intel64"

    $null = Get-Command ifx -ErrorAction Stop
    $null = Get-Command cl -ErrorAction Stop

    if (-not $SkipBuild) {
        Write-Step "Preparing vcpkg"
        $env:VCPKG_ROOT = $VcpkgRoot
        if (-not (Test-Path -Path $VcpkgRoot)) {
            Invoke-Checked -FilePath $gitCmd -Arguments @("clone", "https://github.com/microsoft/vcpkg.git", $VcpkgRoot)
        }

        $vcpkgExe = Join-Path $VcpkgRoot "vcpkg.exe"
        if (-not (Test-Path -Path $vcpkgExe -PathType Leaf)) {
            $bootstrap = Join-Path $VcpkgRoot "bootstrap-vcpkg.bat"
            Invoke-Checked -FilePath "cmd.exe" -Arguments @("/c", "`"$bootstrap`" -disableMetrics")
        }

        Invoke-Checked -FilePath $vcpkgExe -Arguments @(
            "install",
            "--triplet", $Triplet,
            "--x-manifest-root", $RepoRoot,
            "--x-install-root", $vcpkgInstallRoot
        )

        $tclCandidates = @(
            (Join-Path $vcpkgInstallRoot "$Triplet/tools/tcl/tcl8.6"),
            (Join-Path $vcpkgInstallRoot "$Triplet/tools/tcl/tcl9.0"),
            (Join-Path $vcpkgInstallRoot "$Triplet/tools/tcl/lib/tcl8.6"),
            (Join-Path $vcpkgInstallRoot "$Triplet/tools/tcl/lib/tcl9.0"),
            (Join-Path $VcpkgRoot "installed/$Triplet/tools/tcl/tcl8.6"),
            (Join-Path $VcpkgRoot "installed/$Triplet/tools/tcl/tcl9.0"),
            (Join-Path $VcpkgRoot "installed/$Triplet/tools/tcl/lib/tcl8.6"),
            (Join-Path $VcpkgRoot "installed/$Triplet/tools/tcl/lib/tcl9.0")
        )

        $tclIncludeDir = $null
        $tclIncludeCandidates = @(
            (Join-Path $vcpkgInstallRoot "$Triplet/include"),
            (Join-Path $VcpkgRoot "installed/$Triplet/include")
        )
        foreach ($candidate in $tclIncludeCandidates) {
            if (Test-Path -Path (Join-Path $candidate "tcl.h")) {
                $tclIncludeDir = (Resolve-Path $candidate).Path
                break
            }
        }
        if (-not $tclIncludeDir) {
            throw "Could not locate Tcl include directory with tcl.h."
        }
        $tclHeaderPath = Join-Path $tclIncludeDir "tcl.h"
        if (Test-Path -Path $tclHeaderPath -PathType Leaf) {
            $patchLine = Select-String -Path $tclHeaderPath -Pattern '^\s*#define\s+TCL_PATCH_LEVEL\s+"([0-9.]+)"' |
                Select-Object -First 1
            if ($patchLine -and $patchLine.Matches.Count -gt 0) {
                $tclPatchLevel = $patchLine.Matches[0].Groups[1].Value
            }
        }

        $tclLibrary = $null
        $tclLibSearchRoots = @(
            (Join-Path $vcpkgInstallRoot "$Triplet/lib"),
            (Join-Path $VcpkgRoot "installed/$Triplet/lib")
        )
        foreach ($root in $tclLibSearchRoots) {
            if (-not (Test-Path -Path $root)) {
                continue
            }
            $candidate = Get-ChildItem -Path $root -File -Filter "tcl*ts.lib" -ErrorAction SilentlyContinue |
                Where-Object { $_.Name -notlike "*g.lib" } |
                Sort-Object Name |
                Select-Object -First 1
            if ($candidate) {
                $tclLibrary = $candidate.FullName
                break
            }
        }
        if (-not $tclLibrary) {
            throw "Could not locate Tcl static library (expected pattern: tcl*ts.lib)."
        }

        $tclshExe = Get-ChildItem -Path (Join-Path $vcpkgInstallRoot "$Triplet/tools/tcl") -Recurse -File -Filter "tclsh*.exe" -ErrorAction SilentlyContinue |
            Sort-Object FullName |
            Select-Object -First 1

        if ($tclLibrary -match "tcl(\d)(\d)") {
            $tclRuntimeName = "tcl$($Matches[1]).$($Matches[2])"
        }

        $tclRuntimeCandidates = @(
            (Join-Path $vcpkgInstallRoot "$Triplet/tools/tcl/$tclRuntimeName"),
            (Join-Path $vcpkgInstallRoot "$Triplet/tools/tcl/lib/$tclRuntimeName"),
            (Join-Path $VcpkgRoot "installed/$Triplet/tools/tcl/$tclRuntimeName"),
            (Join-Path $VcpkgRoot "installed/$Triplet/tools/tcl/lib/$tclRuntimeName")
        )

        $tclBuildtreeSrcRoot = Join-Path $VcpkgRoot "buildtrees/tcl/src"
        if (Test-Path -Path $tclBuildtreeSrcRoot -PathType Container) {
            foreach ($srcDir in (Get-ChildItem -Path $tclBuildtreeSrcRoot -Directory -ErrorAction SilentlyContinue)) {
                $tclHeader = Join-Path $srcDir.FullName "generic/tcl.h"
                $tclLibDir = Join-Path $srcDir.FullName "library"
                if (-not (Test-Path -Path $tclHeader -PathType Leaf) -or -not (Test-Path -Path (Join-Path $tclLibDir "init.tcl") -PathType Leaf)) {
                    continue
                }

                $majorMatch = Select-String -Path $tclHeader -Pattern "^\s*#define\s+TCL_MAJOR_VERSION\s+([0-9]+)" |
                    Select-Object -First 1
                $minorMatch = Select-String -Path $tclHeader -Pattern "^\s*#define\s+TCL_MINOR_VERSION\s+([0-9]+)" |
                    Select-Object -First 1
                if ($majorMatch -and $minorMatch) {
                    $major = [regex]::Match($majorMatch.Line, "([0-9]+)$").Groups[1].Value
                    $minor = [regex]::Match($minorMatch.Line, "([0-9]+)$").Groups[1].Value
                    if ("tcl$major.$minor" -eq $tclRuntimeName) {
                        $tclRuntimeCandidates += $tclLibDir
                    }
                }
            }
        }

        $tclRuntimeCandidates += $tclCandidates
        foreach ($candidate in ($tclRuntimeCandidates | Where-Object { $_ } | Select-Object -Unique)) {
            if (Test-Path -Path (Join-Path $candidate "init.tcl")) {
                $tclRuntimeDir = (Resolve-Path $candidate).Path
                break
            }
        }
        if (-not $tclRuntimeDir) {
            throw "Tcl runtime directory with init.tcl was not found for runtime '$tclRuntimeName'."
        }

        if (-not $SkipMumps) {
            Write-Step "Preparing and building MUMPS"
            if (-not (Test-Path -Path $MumpsRoot)) {
                Invoke-Checked -FilePath $gitCmd -Arguments @("clone", "https://github.com/OpenSees/mumps.git", $MumpsRoot)
                if ($MumpsCommit) {
                    Invoke-Checked -FilePath $gitCmd -Arguments @("-C", $MumpsRoot, "checkout", $MumpsCommit)
                }
            }

            $mumpsBuildDir = Join-Path $MumpsRoot "build"
            if (Test-Path -Path $mumpsBuildDir -PathType Container) {
                $mumpsBuildDirResolved = [System.IO.Path]::GetFullPath($mumpsBuildDir)
                $mumpsRootResolved = [System.IO.Path]::GetFullPath($MumpsRoot)
                $resetMumpsCache = $false
                foreach ($cache in Get-ChildItem -Path $mumpsBuildDir -Recurse -File -Filter "CMakeCache.txt" -ErrorAction SilentlyContinue) {
                    foreach ($entry in Get-Content -Path $cache.FullName) {
                        if ($entry -like "CMAKE_CACHEFILE_DIR:INTERNAL=*") {
                            $cachedBuildDir = [System.IO.Path]::GetFullPath(($entry.Split("=", 2)[1]).Trim())
                            if (-not $cachedBuildDir.StartsWith($mumpsBuildDirResolved, [System.StringComparison]::OrdinalIgnoreCase)) {
                                $resetMumpsCache = $true
                                break
                            }
                        } elseif ($entry -like "CMAKE_HOME_DIRECTORY:INTERNAL=*") {
                            $cachedSourceDir = [System.IO.Path]::GetFullPath(($entry.Split("=", 2)[1]).Trim())
                            if (-not $cachedSourceDir.StartsWith($mumpsRootResolved, [System.StringComparison]::OrdinalIgnoreCase)) {
                                $resetMumpsCache = $true
                                break
                            }
                        }
                    }
                    if ($resetMumpsCache) {
                        break
                    }
                }

                if ($resetMumpsCache) {
                    Write-Log "Resetting stale MUMPS build tree after source/build directory move." "WARN"
                    Remove-Item -Path $mumpsBuildDir -Recurse -Force
                }
            }
            Invoke-Checked -FilePath $cmakeCmd -Arguments @(
                "-S", $MumpsRoot,
                "-B", $mumpsBuildDir,
                "-G", "Ninja",
                "-Darith=d",
                "-DCMAKE_BUILD_TYPE=Release",
                "-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded",
                "-DCMAKE_Fortran_COMPILER=ifx"
            )
            Invoke-Checked -FilePath $cmakeCmd -Arguments @(
                "--build", $mumpsBuildDir,
                "--config", "Release",
                "--parallel", "$Parallel"
            )
        } else {
            $mumpsBuildDir = Join-Path $MumpsRoot "build"
        }

        if (-not (Test-Path -Path $mumpsBuildDir)) {
            throw "MUMPS build directory does not exist: $mumpsBuildDir"
        }
        $mumpsBuildDir = (Resolve-Path $mumpsBuildDir).Path

        Write-Step "Preparing MKL library arguments"
        if (-not $env:MKLROOT) {
            throw "MKLROOT was not set by setvars.bat"
        }
        $mklInterfaceFull = "intel_lp64"
        $mklInterfaceSuffix = "lp64"

        $mklLibPath = Join-Path $env:MKLROOT "lib\intel64"
        if (-not (Test-Path -Path (Join-Path $mklLibPath "mkl_core.lib"))) {
            $mklLibPath = Join-Path $env:MKLROOT "lib"
        }
        if (-not (Test-Path -Path (Join-Path $mklLibPath "mkl_core.lib"))) {
            throw "Could not locate mkl_core.lib under MKLROOT: $($env:MKLROOT)"
        }

        $scalapackLibList = @(
            (Join-Path $mklLibPath "mkl_scalapack_$mklInterfaceSuffix.lib"),
            (Join-Path $mklLibPath "mkl_intel_$mklInterfaceSuffix.lib"),
            (Join-Path $mklLibPath "mkl_sequential.lib"),
            (Join-Path $mklLibPath "mkl_core.lib"),
            (Join-Path $mklLibPath "mkl_blacs_intelmpi_$mklInterfaceSuffix.lib")
        )
        $lapackLibList = @(
            (Join-Path $mklLibPath "mkl_intel_$mklInterfaceSuffix.lib"),
            (Join-Path $mklLibPath "mkl_sequential.lib"),
            (Join-Path $mklLibPath "mkl_core.lib")
        )
        foreach ($libPath in ($scalapackLibList + $lapackLibList | Select-Object -Unique)) {
            if (-not (Test-Path -Path $libPath)) {
                throw "Expected MKL library not found: $libPath"
            }
        }

        $scalapackLibraries = $scalapackLibList -join ";"
        $lapackLibraries = $lapackLibList -join ";"

        Write-Step "Configuring OpenSees"
        $toolchainFile = Join-Path $VcpkgRoot "scripts\buildsystems\vcpkg.cmake"
        if (-not (Test-Path -Path $toolchainFile -PathType Leaf)) {
            throw "Vcpkg toolchain file not found: $toolchainFile"
        }

        $cacheFile = Join-Path $BuildDir "CMakeCache.txt"
        if (Test-Path -Path $cacheFile -PathType Leaf) {
            $generatorLine = Get-Content -Path $cacheFile |
                Where-Object { $_ -like "CMAKE_GENERATOR:INTERNAL=*" } |
                Select-Object -First 1

            if ($generatorLine) {
                $cachedGenerator = ($generatorLine.Split("=", 2)[1]).Trim()
                if ($cachedGenerator -ne "Ninja") {
                    Write-Log "Existing build cache uses '$cachedGenerator'. Resetting CMake cache for Ninja." "WARN"
                    Remove-Item -Path $cacheFile -Force
                    $cacheDir = Join-Path $BuildDir "CMakeFiles"
                    if (Test-Path -Path $cacheDir) {
                        Remove-Item -Path $cacheDir -Recurse -Force
                    }
                }
            }
        }

        $cmakeConfigureArgs = @(
            "-S", $RepoRoot,
            "-B", $BuildDir,
            "-G", "Ninja",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_NINJA_FORCE_RESPONSE_FILE=ON",
            "-DCMAKE_C_USE_RESPONSE_FILE_FOR_INCLUDES=ON",
            "-DCMAKE_CXX_USE_RESPONSE_FILE_FOR_INCLUDES=ON",
            "-DCMAKE_Fortran_USE_RESPONSE_FILE_FOR_INCLUDES=ON",
            "-DCMAKE_C_COMPILER=cl",
            "-DCMAKE_CXX_COMPILER=cl",
            "-DCMAKE_Fortran_COMPILER=ifx",
            "-DCMAKE_TOOLCHAIN_FILE=$toolchainFile",
            "-DVCPKG_TARGET_TRIPLET=$Triplet",
            "-DVCPKG_MANIFEST_DIR=$RepoRoot",
            "-DVCPKG_INSTALLED_DIR=$vcpkgInstallRoot",
            "-DPython_EXECUTABLE=$pythonCmd",
            "-DMUMPS_DIR=$mumpsBuildDir",
            "-DOPS_TCL_RUNTIME_DIR=$tclRuntimeDir",
            "-DTCL_INCLUDE_PATH=$tclIncludeDir",
            "-DTCL_LIBRARY=$tclLibrary",
            "-DBLA_STATIC=ON",
            "-DMKL_LINK=static",
            "-DMKL_INTERFACE_FULL=$mklInterfaceFull",
            "-DSCALAPACK_LIBRARIES=$scalapackLibraries",
            "-DLAPACK_LIBRARIES=$lapackLibraries"
        )
        if ($tclshExe) {
            $cmakeConfigureArgs += "-DTCL_TCLSH=$($tclshExe.FullName)"
        }
        Invoke-Checked -FilePath $cmakeCmd -Arguments $cmakeConfigureArgs

        Write-Step "Building OpenSees targets"
        Stop-StaleOpenSeesProcesses
        Invoke-Checked -FilePath $cmakeCmd -Arguments @(
            "--build", $BuildDir,
            "--target", "OpenSees", "OpenSeesPy", "OpenSeesSP", "OpenSeesMP",
            "--parallel", "$Parallel"
        )
    } else {
        Write-Step "Skipping configure/build stages (-SkipBuild)"
        $existingTclCandidates = @(
            (Join-Path $BuildDir "lib\tcl8.6"),
            (Join-Path $BuildDir "lib\tcl9.0")
        )
        foreach ($candidate in $existingTclCandidates) {
            if (Test-Path -Path (Join-Path $candidate "init.tcl") -PathType Leaf) {
                $tclRuntimeDir = (Resolve-Path $candidate).Path
                $tclRuntimeName = Split-Path -Path $tclRuntimeDir -Leaf
                break
            }
        }
        if ($tclRuntimeDir) {
            Write-Log "Using existing Tcl runtime from build directory: $tclRuntimeDir"
        } else {
            Write-Log "No Tcl runtime found under build/lib. Tests may fail unless TCL_LIBRARY is already valid." "WARN"
        }
    }

    Write-Step "Locating build artifacts"
    $openSeesExe = Get-LatestArtifact -RootDir $BuildDir -FileName "OpenSees.exe"
    $openSeesSpExe = Get-LatestArtifact -RootDir $BuildDir -FileName "OpenSeesSP.exe"
    $openSeesMpExe = Get-LatestArtifact -RootDir $BuildDir -FileName "OpenSeesMP.exe"
    $openSeesPyDll = Get-LatestArtifact -RootDir $BuildDir -FileName "OpenSeesPy.dll"

    $artifactDir = $openSeesExe.DirectoryName
    $openSeesPyPyd = Join-Path $artifactDir "opensees.pyd"
    Copy-Item -Path $openSeesPyDll.FullName -Destination $openSeesPyPyd -Force

    if ($tclRuntimeDir) {
        Write-Step "Ensuring Tcl scripts are present under build/lib/$tclRuntimeName"
        $tclTargetDir = Join-Path $BuildDir ("lib\" + $tclRuntimeName)
        New-Item -ItemType Directory -Path $tclTargetDir -Force | Out-Null
        $resolvedRuntime = (Resolve-Path $tclRuntimeDir).Path
        $resolvedTarget = (Resolve-Path $tclTargetDir).Path
        if ($resolvedRuntime -ieq $resolvedTarget) {
            Write-Log "Tcl runtime already present at target location; skipping copy."
        } else {
            Copy-Item -Path (Join-Path $tclRuntimeDir "*") -Destination $tclTargetDir -Recurse -Force
        }
    } else {
        $tclTargetDir = $null
        foreach ($candidate in @(
            (Join-Path $BuildDir "lib\tcl8.6"),
            (Join-Path $BuildDir "lib\tcl9.0")
        )) {
            if (Test-Path -Path (Join-Path $candidate "init.tcl") -PathType Leaf) {
                $tclTargetDir = (Resolve-Path $candidate).Path
                break
            }
        }
    }

    if ($tclTargetDir) {
        if (-not $tclPatchLevel) {
            $tclHeaderFallback = Join-Path $vcpkgInstallRoot "$Triplet\include\tcl.h"
            if (Test-Path -Path $tclHeaderFallback -PathType Leaf) {
                $patchLine = Select-String -Path $tclHeaderFallback -Pattern '^\s*#define\s+TCL_PATCH_LEVEL\s+"([0-9.]+)"' |
                    Select-Object -First 1
                if ($patchLine -and $patchLine.Matches.Count -gt 0) {
                    $tclPatchLevel = $patchLine.Matches[0].Groups[1].Value
                }
            }
        }
        if (-not $tclPatchLevel) {
            $tclPatchLevel = "8.6.10"
        }

        $tclUserFallbackDir = Join-Path $env:USERPROFILE ("tcl" + $tclPatchLevel + "\library")
        try {
            New-Item -ItemType Directory -Path $tclUserFallbackDir -Force | Out-Null
            Copy-Item -Path (Join-Path $tclTargetDir "*") -Destination $tclUserFallbackDir -Recurse -Force
            Write-Log "Installed Tcl fallback runtime for desktop launch: $tclUserFallbackDir"
        } catch {
            Write-Log "Failed to refresh Tcl user fallback runtime at '$tclUserFallbackDir': $($_.Exception.Message)" "WARN"
        }
    }

    Write-Step "Creating OpenSees desktop launchers"
    $setvarsBat = $null
    if ($env:ONEAPI_ROOT) {
        $setvarsCandidate = Join-Path $env:ONEAPI_ROOT "setvars.bat"
        if (Test-Path -Path $setvarsCandidate -PathType Leaf) {
            $setvarsBat = (Resolve-Path $setvarsCandidate).Path
        }
    }

    $parallelExampleDir = Join-Path $RepoRoot "EXAMPLES\ParallelModelMP"
    if (-not (Test-Path -Path $parallelExampleDir -PathType Container)) {
        $parallelExampleDir = $null
    }

    $launcherVsInstallDir = Get-Vs2022InstallDir

    $baseLauncherLines = @(
        "@echo off",
        "setlocal"
    )
    if ($launcherVsInstallDir) {
        $baseLauncherLines += "set `"VS2022INSTALLDIR=$launcherVsInstallDir`""
    }
    if ($setvarsBat) {
        $baseLauncherLines += "if exist `"$setvarsBat`" call `"$setvarsBat`" intel64 >nul 2>&1"
    }
    if ($tclTargetDir) {
        $baseLauncherLines += "set `"TCL_LIBRARY=$tclTargetDir`""
    }

    $serialLauncherPath = Join-Path $artifactDir "OpenSees-Launch.cmd"
    $serialLauncherLines = @($baseLauncherLines)
    $serialLauncherLines += "`"$($openSeesExe.FullName)`" %*"
    $serialLauncherLines += "endlocal"
    Set-Content -Path $serialLauncherPath -Value ($serialLauncherLines -join "`r`n") -Encoding Ascii
    Write-Log "Created launcher: $serialLauncherPath"

    $ladrunoSerialPath = Join-Path $artifactDir "opensees_ladruno.bat"
    Set-Content -Path $ladrunoSerialPath -Value ((New-LadrunoWrapperLines -Mode "serial") -join "`r`n") -Encoding Ascii
    Write-Log "Created launcher: $ladrunoSerialPath"

    $spLauncherPath = Join-Path $artifactDir "OpenSeesSP-Launch.cmd"
    $spLauncherLines = @($baseLauncherLines)
    $spLauncherLines += @(
        "if `"%OPS_MPI_N%`"==`"`" set `"OPS_MPI_N=2`"",
        "set `"OPS_LOG=%~dp0OpenSeesSP-launch.log`"",
        "where mpiexec >nul 2>&1",
        "if errorlevel 1 (",
        "  echo ERROR: mpiexec not found. Install/initialize Intel MPI or run from oneAPI shell.",
        "  echo ERROR: mpiexec not found.>> `"%OPS_LOG%`"",
        "  pause",
        "  endlocal & exit /b 1",
        ")",
        "echo [%date% %time%] OpenSeesSP launcher start > `"%OPS_LOG%`"",
        "echo OPS_MPI_N=%OPS_MPI_N%>> `"%OPS_LOG%`""
    )
    if ($parallelExampleDir) {
        $spLauncherLines += @(
            "if not `"%~1`"==`"`" goto run_args",
            "if not exist `"$parallelExampleDir`" goto no_default",
            "pushd `"$parallelExampleDir`"",
            "mpiexec -localonly -n %OPS_MPI_N% `"$($openSeesSpExe.FullName)`" exampleSP.tcl >> `"%OPS_LOG%`" 2>&1",
            "set `"OPS_EXIT=%ERRORLEVEL%`"",
            "if not `"%OPS_EXIT%`"==`"0`" (",
            "  mpiexec -launcher service -localonly -n %OPS_MPI_N% `"$($openSeesSpExe.FullName)`" exampleSP.tcl >> `"%OPS_LOG%`" 2>&1",
            "  set `"OPS_EXIT=%ERRORLEVEL%`"",
            ")",
            "popd",
            "goto finalize",
            ":run_args",
            "mpiexec -localonly -n %OPS_MPI_N% `"$($openSeesSpExe.FullName)`" %* >> `"%OPS_LOG%`" 2>&1",
            "set `"OPS_EXIT=%ERRORLEVEL%`"",
            "if not `"%OPS_EXIT%`"==`"0`" (",
            "  mpiexec -launcher service -localonly -n %OPS_MPI_N% `"$($openSeesSpExe.FullName)`" %* >> `"%OPS_LOG%`" 2>&1",
            "  set `"OPS_EXIT=%ERRORLEVEL%`"",
            ")",
            "goto finalize",
            ":no_default",
            "echo ERROR: default SP example folder not found: $parallelExampleDir",
            "echo Usage: OpenSeesSP-Launch.cmd your_parallel_script.tcl",
            "set `"OPS_EXIT=1`"",
            "goto finalize",
            ":finalize",
            "if `"%OPS_EXIT%`"==`"0`" (",
            "  endlocal & exit /b 0",
            ")",
            "echo ERROR: OpenSeesSP MPI launch failed (exit %OPS_EXIT%).",
            "echo Full log: %OPS_LOG%",
            "echo Hint: run from Intel oneAPI command prompt and test: mpiexec -n 2 hostname",
            "echo Hint: PowerShell profile errors (for example missing 'mise') can break Intel MPI bootstrap.",
            "echo Hint: run as Administrator: powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\fix_intel_mpi_windows11.ps1 -Ranks 2 -SetWinRmAutomatic",
            "if exist `"%OPS_LOG%`" type `"%OPS_LOG%`"",
            "pause",
            "endlocal & exit /b %OPS_EXIT%"
        )
    } else {
        $spLauncherLines += @(
            "if not `"%~1`"==`"`" goto run_args",
            "echo Usage: OpenSeesSP-Launch.cmd your_parallel_script.tcl",
            "echo No default SP example path was detected.",
            "set `"OPS_EXIT=1`"",
            "goto finalize",
            ":run_args",
            "mpiexec -localonly -n %OPS_MPI_N% `"$($openSeesSpExe.FullName)`" %* >> `"%OPS_LOG%`" 2>&1",
            "set `"OPS_EXIT=%ERRORLEVEL%`"",
            "if not `"%OPS_EXIT%`"==`"0`" (",
            "  mpiexec -launcher service -localonly -n %OPS_MPI_N% `"$($openSeesSpExe.FullName)`" %* >> `"%OPS_LOG%`" 2>&1",
            "  set `"OPS_EXIT=%ERRORLEVEL%`"",
            ")",
            "goto finalize",
            ":finalize",
            "if `"%OPS_EXIT%`"==`"0`" (",
            "  endlocal & exit /b 0",
            ")",
            "echo ERROR: OpenSeesSP MPI launch failed (exit %OPS_EXIT%).",
            "echo Full log: %OPS_LOG%",
            "echo Hint: run from Intel oneAPI command prompt and test: mpiexec -n 2 hostname",
            "echo Hint: PowerShell profile errors (for example missing 'mise') can break Intel MPI bootstrap.",
            "echo Hint: run as Administrator: powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\fix_intel_mpi_windows11.ps1 -Ranks 2 -SetWinRmAutomatic",
            "if exist `"%OPS_LOG%`" type `"%OPS_LOG%`"",
            "pause",
            "endlocal & exit /b %OPS_EXIT%"
        )
    }
    Set-Content -Path $spLauncherPath -Value ($spLauncherLines -join "`r`n") -Encoding Ascii
    Write-Log "Created launcher: $spLauncherPath"

    $ladrunoSpPath = Join-Path $artifactDir "opensees_sp_ladruno.bat"
    Set-Content -Path $ladrunoSpPath -Value ((New-LadrunoWrapperLines -Mode "sp") -join "`r`n") -Encoding Ascii
    Write-Log "Created launcher: $ladrunoSpPath"

    $mpLauncherPath = Join-Path $artifactDir "OpenSeesMP-Launch.cmd"
    $mpLauncherLines = @($baseLauncherLines)
    $mpLauncherLines += @(
        "if `"%OPS_MPI_N%`"==`"`" set `"OPS_MPI_N=2`"",
        "set `"OPS_LOG=%~dp0OpenSeesMP-launch.log`"",
        "where mpiexec >nul 2>&1",
        "if errorlevel 1 (",
        "  echo ERROR: mpiexec not found. Install/initialize Intel MPI or run from oneAPI shell.",
        "  echo ERROR: mpiexec not found.>> `"%OPS_LOG%`"",
        "  pause",
        "  endlocal & exit /b 1",
        ")",
        "echo [%date% %time%] OpenSeesMP launcher start > `"%OPS_LOG%`"",
        "echo OPS_MPI_N=%OPS_MPI_N%>> `"%OPS_LOG%`""
    )
    if ($parallelExampleDir) {
        $mpLauncherLines += @(
            "if not `"%~1`"==`"`" goto run_args",
            "if not exist `"$parallelExampleDir`" goto no_default",
            "pushd `"$parallelExampleDir`"",
            "mpiexec -localonly -n %OPS_MPI_N% `"$($openSeesMpExe.FullName)`" exampleMP.tcl >> `"%OPS_LOG%`" 2>&1",
            "set `"OPS_EXIT=%ERRORLEVEL%`"",
            "if not `"%OPS_EXIT%`"==`"0`" (",
            "  mpiexec -launcher service -localonly -n %OPS_MPI_N% `"$($openSeesMpExe.FullName)`" exampleMP.tcl >> `"%OPS_LOG%`" 2>&1",
            "  set `"OPS_EXIT=%ERRORLEVEL%`"",
            ")",
            "popd",
            "goto finalize",
            ":run_args",
            "mpiexec -localonly -n %OPS_MPI_N% `"$($openSeesMpExe.FullName)`" %* >> `"%OPS_LOG%`" 2>&1",
            "set `"OPS_EXIT=%ERRORLEVEL%`"",
            "if not `"%OPS_EXIT%`"==`"0`" (",
            "  mpiexec -launcher service -localonly -n %OPS_MPI_N% `"$($openSeesMpExe.FullName)`" %* >> `"%OPS_LOG%`" 2>&1",
            "  set `"OPS_EXIT=%ERRORLEVEL%`"",
            ")",
            "goto finalize",
            ":no_default",
            "echo ERROR: default MP example folder not found: $parallelExampleDir",
            "echo Usage: OpenSeesMP-Launch.cmd your_parallel_script.tcl",
            "set `"OPS_EXIT=1`"",
            "goto finalize",
            ":finalize",
            "if `"%OPS_EXIT%`"==`"0`" (",
            "  endlocal & exit /b 0",
            ")",
            "echo ERROR: OpenSeesMP MPI launch failed (exit %OPS_EXIT%).",
            "echo Full log: %OPS_LOG%",
            "echo Hint: run from Intel oneAPI command prompt and test: mpiexec -n 2 hostname",
            "echo Hint: PowerShell profile errors (for example missing 'mise') can break Intel MPI bootstrap.",
            "echo Hint: run as Administrator: powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\fix_intel_mpi_windows11.ps1 -Ranks 2 -SetWinRmAutomatic",
            "if exist `"%OPS_LOG%`" type `"%OPS_LOG%`"",
            "pause",
            "endlocal & exit /b %OPS_EXIT%"
        )
    } else {
        $mpLauncherLines += @(
            "if not `"%~1`"==`"`" goto run_args",
            "echo Usage: OpenSeesMP-Launch.cmd your_parallel_script.tcl",
            "echo No default MP example path was detected.",
            "set `"OPS_EXIT=1`"",
            "goto finalize",
            ":run_args",
            "mpiexec -localonly -n %OPS_MPI_N% `"$($openSeesMpExe.FullName)`" %* >> `"%OPS_LOG%`" 2>&1",
            "set `"OPS_EXIT=%ERRORLEVEL%`"",
            "if not `"%OPS_EXIT%`"==`"0`" (",
            "  mpiexec -launcher service -localonly -n %OPS_MPI_N% `"$($openSeesMpExe.FullName)`" %* >> `"%OPS_LOG%`" 2>&1",
            "  set `"OPS_EXIT=%ERRORLEVEL%`"",
            ")",
            "goto finalize",
            ":finalize",
            "if `"%OPS_EXIT%`"==`"0`" (",
            "  endlocal & exit /b 0",
            ")",
            "echo ERROR: OpenSeesMP MPI launch failed (exit %OPS_EXIT%).",
            "echo Full log: %OPS_LOG%",
            "echo Hint: run from Intel oneAPI command prompt and test: mpiexec -n 2 hostname",
            "echo Hint: PowerShell profile errors (for example missing 'mise') can break Intel MPI bootstrap.",
            "echo Hint: run as Administrator: powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\fix_intel_mpi_windows11.ps1 -Ranks 2 -SetWinRmAutomatic",
            "if exist `"%OPS_LOG%`" type `"%OPS_LOG%`"",
            "pause",
            "endlocal & exit /b %OPS_EXIT%"
        )
    }
    Set-Content -Path $mpLauncherPath -Value ($mpLauncherLines -join "`r`n") -Encoding Ascii
    Write-Log "Created launcher: $mpLauncherPath"

    $ladrunoMpPath = Join-Path $artifactDir "opensees_mp_ladruno.bat"
    Set-Content -Path $ladrunoMpPath -Value ((New-LadrunoWrapperLines -Mode "mp") -join "`r`n") -Encoding Ascii
    Write-Log "Created launcher: $ladrunoMpPath"

    Write-Step "Copying oneAPI runtime DLLs"
    $dependencyDlls = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
    $knownRuntimeDlls = @(
        "libiomp5md.dll",
        "libifcoremd.dll",
        "libifportmd.dll",
        "libmmd.dll",
        "svml_dispmd.dll",
        "impi.dll"
    )
    foreach ($name in $knownRuntimeDlls) {
        $null = $dependencyDlls.Add($name)
    }

    foreach ($binary in @($openSeesExe.FullName, $openSeesSpExe.FullName, $openSeesMpExe.FullName, $openSeesPyPyd)) {
        foreach ($dep in (Get-DumpbinDependencies -BinaryPath $binary)) {
            $null = $dependencyDlls.Add($dep)
        }
    }

    if (-not $env:ONEAPI_ROOT) {
        throw "ONEAPI_ROOT was not set by setvars.bat"
    }

    $oneApiRoots = @(
        (Join-Path $env:ONEAPI_ROOT "compiler\latest\bin"),
        (Join-Path $env:ONEAPI_ROOT "compiler\latest\redist\intel64_win\compiler"),
        (Join-Path $env:ONEAPI_ROOT "mpi\latest\bin"),
        (Join-Path $env:ONEAPI_ROOT "mpi\latest\bin\release"),
        (Join-Path $env:ONEAPI_ROOT "mkl\latest\redist\intel64"),
        (Join-Path $env:MKLROOT "redist\intel64")
    ) | Where-Object { $_ -and (Test-Path -Path $_) }

    foreach ($dllName in $dependencyDlls) {
        $destPath = Join-Path $artifactDir $dllName
        if (Test-Path -Path $destPath) {
            continue
        }

        $copied = $false
        foreach ($root in $oneApiRoots) {
            $match = Get-ChildItem -Path $root -Recurse -File -Filter $dllName -ErrorAction SilentlyContinue |
                Select-Object -First 1
            if ($match) {
                Copy-Item -Path $match.FullName -Destination $destPath -Force
                $copied = $true
                break
            }
        }

        if (-not $copied -and ($knownRuntimeDlls -contains $dllName)) {
            Write-Log "Could not locate runtime DLL '$dllName' in oneAPI folders." "WARN"
        }
    }

    if (-not $SkipTests) {
        Write-Step "Running smoke tests (mode: $SmokeMode)"
        Stop-StaleOpenSeesProcesses

        $examplePython = Join-Path $RepoRoot "EXAMPLES\ExamplePython\example_variable_analysis.py"
        $serialSmokeScript = Join-Path $BuildDir "smoke_serial.tcl"

        if ($tclTargetDir) {
            $env:TCL_LIBRARY = $tclTargetDir
            Write-Log "TCL_LIBRARY set to: $tclTargetDir"
        } else {
            Write-Log "TCL_LIBRARY was not set because Tcl runtime was not found." "WARN"
        }

        Set-Content -Path $serialSmokeScript -Encoding ascii -Value @'
model basic -ndm 1 -ndf 1
node 1 0.0
node 2 1.0
fix 1 1
uniaxialMaterial Elastic 1 1000.0
element truss 1 1 2 1.0 1
timeSeries Linear 1
pattern Plain 1 1 {
    load 2 1.0
}
constraints Plain
numberer RCM
system BandGeneral
test NormUnbalance 1.0e-8 10 0
algorithm Newton
integrator LoadControl 1.0
analysis Static
analyze 1
puts "serial-smoke-ok"
wipe
exit
'@
        Invoke-Checked -FilePath $openSeesExe.FullName -Arguments @($serialSmokeScript)

        if ($env:PYTHONPATH) {
            $env:PYTHONPATH = "$artifactDir;$env:PYTHONPATH"
        } else {
            $env:PYTHONPATH = $artifactDir
        }

        Invoke-Checked -FilePath $pythonCmd -Arguments @("-c", "import opensees; print('opensees import ok')")
        if ($SmokeMode -eq "full") {
            try {
                Invoke-Checked -FilePath $pythonCmd -Arguments @($examplePython) -TimeoutSec $SmokeTimeoutSec
            } catch {
                Write-Log "Optional Python variable analysis smoke test failed: $($_.Exception.Message)" "WARN"
            }
        } else {
            Write-Log "Skipping optional Python variable analysis smoke test in quick mode."
        }

        $mpiexec = Get-Command mpiexec -ErrorAction SilentlyContinue
        if ($mpiexec) {
            if ($SmokeMode -eq "quick") {
                Write-Step "Running quick MPI smoke tests (ParallelModelMP)"
                $mpiWorkDir = Join-Path $RepoRoot "EXAMPLES\ParallelModelMP"

                Stop-StaleOpenSeesProcesses
                Invoke-MpiSmoke -FilePath $mpiexec.Source -Arguments @("-n", "2", $openSeesSpExe.FullName, "exampleSP.tcl") -WorkingDirectory $mpiWorkDir -TimeoutSec $SmokeTimeoutSec

                Stop-StaleOpenSeesProcesses
                Invoke-MpiSmoke -FilePath $mpiexec.Source -Arguments @("-n", "2", $openSeesMpExe.FullName, "exampleMP.tcl") -WorkingDirectory $mpiWorkDir -TimeoutSec $SmokeTimeoutSec
            } else {
                Write-Step "Running full MPI smoke tests (SmallMP)"
                $spWorkDir = Join-Path $RepoRoot "EXAMPLES\ParallelModelMP"
                $mpWorkDir = Join-Path $RepoRoot "EXAMPLES\SmallMP"

                Stop-StaleOpenSeesProcesses
                Invoke-MpiSmoke -FilePath $mpiexec.Source -Arguments @("-n", "2", $openSeesSpExe.FullName, "exampleSP.tcl") -WorkingDirectory $spWorkDir -TimeoutSec $SmokeTimeoutSec

                Stop-StaleOpenSeesProcesses
                Invoke-MpiSmoke -FilePath $mpiexec.Source -Arguments @("-n", "2", $openSeesMpExe.FullName, "Example.tcl") -WorkingDirectory $mpWorkDir -TimeoutSec $SmokeTimeoutSec
            }
        } else {
            Write-Log "mpiexec was not found. Skipping OpenSeesSP/OpenSeesMP runtime tests." "WARN"
        }
    }

    Write-Step "Build completed"
    Write-Log "OpenSees     : $($openSeesExe.FullName)"
    Write-Log "OpenSeesSP   : $($openSeesSpExe.FullName)"
    Write-Log "OpenSeesMP   : $($openSeesMpExe.FullName)"
    Write-Log "OpenSeesPy   : $openSeesPyPyd"
    if ($tclTargetDir) {
        Write-Log "Tcl runtime  : $tclTargetDir"
    }
} finally {
    Stop-StaleOpenSeesProcesses
    Write-Log "Cleanup complete."
}
