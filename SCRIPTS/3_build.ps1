#==============================================================================
#  OpenSees Windows 11 -- Step 3: Build
#
#  Compiles OpenSees (serial + parallel + Python module) using the build
#  harness that was copied into the source tree by Step 2.
#
#  This script automatically initializes the Visual Studio and Intel oneAPI
#  environment if cl/ifx are not already available.  It does this by
#  re-launching itself inside a cmd shell that first calls
#  init_oneapi_windows11.cmd.
#
#  Usage (normal PowerShell from the OpenSees source root):
#    powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\3_build.ps1
#
#  Common overrides:
#    -SkipMumps          MUMPS already built from a previous run
#    -SkipTests          Skip smoke tests after build
#    -SmokeMode full     Run the full smoke-test suite instead of quick
#    -Parallel 4         Limit parallel compilation jobs
#    -DryRun             Print commands without executing
#==============================================================================
[CmdletBinding()]
param(
    [string]$BuildDir     = "build-win11",
    [string]$Triplet      = "x64-windows-static",
    [string]$VcpkgRoot    = "third_party\vcpkg",
    [string]$MumpsRoot    = "third_party\mumps",
    [ValidateSet("quick", "full")][string]$SmokeMode = "quick",
    [int]$SmokeTimeoutSec = 600,
    [int]$Parallel        = 0,
    [switch]$SkipMumps,
    [switch]$SkipBuild,
    [switch]$SkipTests,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# Resolve the repo root (one level above SCRIPTS/).
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
$RepoRoot  = Split-Path -Parent $ScriptDir
$buildScript = Join-Path $ScriptDir "build_windows11_full.ps1"

function Convert-ToPowerShellArgumentString {
    param([System.Collections.IDictionary]$Parameters)

    $parts = @()
    foreach ($entry in $Parameters.GetEnumerator()) {
        $key = [string]$entry.Key
        $value = $entry.Value

        if ($value -is [switch]) {
            if ($value.IsPresent) {
                $parts += "-$key"
            }
            continue
        }

        if ($value -is [bool]) {
            if ($value) {
                $parts += "-$key"
            }
            continue
        }

        if ($null -eq $value) {
            continue
        }

        $escaped = [string]$value -replace '"', '\"'
        $parts += "-$key `"$escaped`""
    }

    return ($parts -join " ")
}

# --------------------------------------------------------------------------
# Check if the build environment (cl, ifx) is already available.
# If not, re-launch this entire script inside a cmd shell that first
# calls init_oneapi_windows11.cmd to set up VS + oneAPI.
# --------------------------------------------------------------------------
$clFound  = $null -ne (Get-Command cl  -ErrorAction SilentlyContinue)
$ifxFound = $null -ne (Get-Command ifx -ErrorAction SilentlyContinue)

if ((-not $clFound) -or (-not $ifxFound)) {
    $initCmd = Join-Path $ScriptDir "init_oneapi_windows11.cmd"
    if (-not (Test-Path $initCmd)) {
        Write-Host "[FAIL] cl/ifx not in PATH and init_oneapi_windows11.cmd not found." -ForegroundColor Red
        Write-Host "  Initialize your VS + oneAPI environment manually, then re-run." -ForegroundColor Yellow
        exit 1
    }

    Write-Host "cl/ifx not found in PATH -- initializing VS + oneAPI environment..." -ForegroundColor Yellow
    Write-Host ""

    $thisScript = $MyInvocation.MyCommand.Definition
    $fwdArgsStr = Convert-ToPowerShellArgumentString -Parameters $PSBoundParameters

    # Launch: cmd -> init_oneapi -> powershell -> this script (with env loaded)
    $cmdLine = if ([string]::IsNullOrWhiteSpace($fwdArgsStr)) {
        "call `"$initCmd`" && powershell -NoProfile -ExecutionPolicy Bypass -File `"$thisScript`""
    } else {
        "call `"$initCmd`" && powershell -NoProfile -ExecutionPolicy Bypass -File `"$thisScript`" $fwdArgsStr"
    }

    if ($DryRun) {
        Write-Host "[DRY RUN] Would execute:" -ForegroundColor Magenta
        Write-Host "  cmd /d /s /c $cmdLine" -ForegroundColor Magenta
        exit 0
    }

    cmd.exe /d /s /c $cmdLine
    exit $LASTEXITCODE
}

# --------------------------------------------------------------------------
# Environment is available -- proceed with the build.
# --------------------------------------------------------------------------
$fwdArgs = @{}
if (-not $PSBoundParameters.ContainsKey("VcpkgRoot")) {
    $fwdArgs["VcpkgRoot"] = $VcpkgRoot
}
if (-not $PSBoundParameters.ContainsKey("MumpsRoot")) {
    $fwdArgs["MumpsRoot"] = $MumpsRoot
}
foreach ($entry in $PSBoundParameters.GetEnumerator()) {
    $fwdArgs[$entry.Key] = $entry.Value
}

if (-not (Test-Path $buildScript)) {
    Write-Host "[FAIL] Build script not found: $buildScript" -ForegroundColor Red
    Write-Host "  Run SCRIPTS\2_fetch_source.ps1 first to copy the build harness." -ForegroundColor Yellow
    exit 1
}

Write-Host ""
Write-Host "Launching build_windows11_full.ps1 with:" -ForegroundColor Cyan
Write-Host "  BuildDir   = $BuildDir"
Write-Host "  Triplet    = $Triplet"
Write-Host "  VcpkgRoot  = $($fwdArgs['VcpkgRoot'])"
Write-Host "  MumpsRoot  = $($fwdArgs['MumpsRoot'])"
Write-Host "  SmokeMode  = $SmokeMode"
Write-Host "  Parallel   = $Parallel"
Write-Host ""

& $buildScript @fwdArgs
$exitCode = $LASTEXITCODE

if ($exitCode -ne 0) {
    Write-Host ""
    Write-Host "[FAIL] Build finished with exit code $exitCode" -ForegroundColor Red
    exit $exitCode
}

Write-Host ""
Write-Host "Build completed successfully." -ForegroundColor Green
Write-Host ""
Write-Host "Next step:" -ForegroundColor Green
Write-Host "  powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\4_package.ps1" -ForegroundColor White
Write-Host ""
