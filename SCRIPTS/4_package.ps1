#==============================================================================
#  OpenSees Windows 11 — Step 4: Package / Installer
#
#  Creates a portable zip and (optionally) a Windows installer from the
#  compiled build output.
#
#  Usage (normal PowerShell from the OpenSees source root):
#    powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\4_package.ps1
#
#  Common overrides:
#    -BuildDir .\build-win11          Where the build output lives
#    -OutputDir .\dist                Where to write the installer / zip
#    -AppName "My Custom OpenSees"    Change the product name
#    -SkipInnoCompile                 Only zip, no .exe installer
#    -DryRun                          Print actions without executing
#==============================================================================
[CmdletBinding()]
param(
    [string]$BuildDir  = "build-win11",
    [string]$OutputDir = "dist",
    [string]$AppName   = "El Ladruno OpenSees",
    [string]$Publisher = "El Ladruno",
    [string]$AppVersion = "",
    [switch]$IncludeExamples,
    [switch]$IncludeAllExamples,
    [switch]$SkipZip,
    [switch]$SkipInnoCompile,
    [string]$InnoCompilerPath = "",
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition

$installerScript = Join-Path $ScriptDir "create_el_ladruno_installer.ps1"

if (-not (Test-Path $installerScript)) {
    Write-Host "[FAIL] Installer script not found: $installerScript" -ForegroundColor Red
    Write-Host "  Run SCRIPTS\2_fetch_source.ps1 first to copy the build harness." -ForegroundColor Yellow
    exit 1
}

# Verify that the build directory has the expected artifacts.
if (-not $DryRun) {
    $requiredExe = Join-Path $BuildDir "OpenSees.exe"
    if (-not (Test-Path $requiredExe)) {
        Write-Host "[FAIL] $requiredExe not found." -ForegroundColor Red
        Write-Host "  Run SCRIPTS\3_build.ps1 first to compile OpenSees." -ForegroundColor Yellow
        exit 1
    }
}

# Build forwarded argument list.
$fwdArgs = @(
    "-BuildDir",  $BuildDir,
    "-OutputDir", $OutputDir,
    "-AppName",   $AppName,
    "-Publisher",  $Publisher
)

if ($AppVersion)         { $fwdArgs += @("-AppVersion", $AppVersion) }
if ($IncludeExamples)    { $fwdArgs += "-IncludeExamples" }
if ($IncludeAllExamples) { $fwdArgs += "-IncludeAllExamples" }
if ($SkipZip)            { $fwdArgs += "-SkipZip" }
if ($SkipInnoCompile)    { $fwdArgs += "-SkipInnoCompile" }
if ($InnoCompilerPath)   { $fwdArgs += @("-InnoCompilerPath", $InnoCompilerPath) }

Write-Host ""
Write-Host "Launching create_el_ladruno_installer.ps1 with:" -ForegroundColor Cyan
Write-Host "  BuildDir  = $BuildDir"
Write-Host "  OutputDir = $OutputDir"
Write-Host "  AppName   = $AppName"
Write-Host ""

& $installerScript @fwdArgs
$exitCode = $LASTEXITCODE

if ($exitCode -ne 0) {
    Write-Host ""
    Write-Host "[FAIL] Packaging finished with exit code $exitCode" -ForegroundColor Red
    exit $exitCode
}

Write-Host ""
Write-Host "Packaging completed.  Output in: $OutputDir" -ForegroundColor Green
Write-Host ""
