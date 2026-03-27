#==============================================================================
#  OpenSees Windows 11 — Step 3: Build
#
#  Compiles OpenSees (serial + parallel + Python module) using the build
#  harness that was copied into the source tree by Step 2.
#
#  Usage (normal PowerShell from the OpenSees source root):
#    powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\3_build.ps1
#
#  All parameters are forwarded to build_windows11_full.ps1.
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
    [string]$VcpkgRoot    = "",
    [string]$MumpsRoot    = "",
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

# Default vcpkg and MUMPS locations match Step 2 layout.
if ([string]::IsNullOrWhiteSpace($VcpkgRoot)) {
    $VcpkgRoot = Join-Path $RepoRoot "third_party\vcpkg"
}
if ([string]::IsNullOrWhiteSpace($MumpsRoot)) {
    $MumpsRoot = Join-Path $RepoRoot "third_party\mumps"
}

# Build the forwarded argument list.
$fwdArgs = @(
    "-BuildDir",       $BuildDir,
    "-Triplet",        $Triplet,
    "-VcpkgRoot",      $VcpkgRoot,
    "-MumpsRoot",      $MumpsRoot,
    "-SmokeMode",      $SmokeMode,
    "-SmokeTimeoutSec", $SmokeTimeoutSec
)

if ($Parallel -gt 0) { $fwdArgs += @("-Parallel", $Parallel) }
if ($SkipMumps)      { $fwdArgs += "-SkipMumps" }
if ($SkipBuild)      { $fwdArgs += "-SkipBuild" }
if ($SkipTests)      { $fwdArgs += "-SkipTests" }
if ($DryRun)         { $fwdArgs += "-DryRun" }

$buildScript = Join-Path $ScriptDir "build_windows11_full.ps1"

if (-not (Test-Path $buildScript)) {
    Write-Host "[FAIL] Build script not found: $buildScript" -ForegroundColor Red
    Write-Host "  Run SCRIPTS\2_fetch_source.ps1 first to copy the build harness." -ForegroundColor Yellow
    exit 1
}

Write-Host ""
Write-Host "Launching build_windows11_full.ps1 with:" -ForegroundColor Cyan
Write-Host "  VcpkgRoot  = $VcpkgRoot"
Write-Host "  MumpsRoot  = $MumpsRoot"
Write-Host "  BuildDir   = $BuildDir"
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
