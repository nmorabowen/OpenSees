<#
.SYNOPSIS
    Build the Ladruno OpenSees Inno Setup installer (setup.exe).

.DESCRIPTION
    Wraps `iscc.exe` (Inno Setup compiler). Run this AFTER
    `Ladruno_scripts\build.bat` has produced dist/.

    Output: Ladruno_files\Ladruno_OpenSees_<version>_setup.exe -- a single
    self-contained installer with a venv-picker wizard.

    Prerequisite (one-time):
        winget install JRSoftware.InnoSetup

.PARAMETER Version
    Version string baked into the installer filename / Add-Remove-Programs
    entry. Default: today's date in YYYYMMDD form.

.PARAMETER Iscc
    Path to iscc.exe. Default: C:\Program Files (x86)\Inno Setup 6\iscc.exe.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File Ladruno_scripts\build_inno_installer.ps1
    powershell -ExecutionPolicy Bypass -File Ladruno_scripts\build_inno_installer.ps1 -Version 2026.05
#>
param(
    [string]$Version = (Get-Date -Format "yyyyMMdd"),
    [string]$Iscc    = "C:\Program Files (x86)\Inno Setup 6\iscc.exe"
)

$ErrorActionPreference = "Stop"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

$scriptDir = $PSScriptRoot
$root      = Split-Path -Parent $scriptDir
$dist      = Join-Path $root "dist"
$out       = Join-Path $root "Ladruno_files"
$iss       = Join-Path $scriptDir "installer.iss"
$banner    = Join-Path $scriptDir "banner_ASCII.txt"
$bannerBom = Join-Path $scriptDir "banner_ASCII_utf8bom.txt"

Write-Host "Ladruno OpenSees - Inno Setup builder" -ForegroundColor Cyan
Write-Host "  version : $Version"
Write-Host "  source  : $dist"
Write-Host "  output  : $out"
Write-Host "  iscc    : $Iscc"
Write-Host ""

# ---- preflight -----------------------------------------------------------
if (-not (Test-Path $Iscc)) {
    Write-Error @"
Inno Setup compiler not found at:
    $Iscc

Install it via:
    winget install JRSoftware.InnoSetup
or pass -Iscc <path-to-iscc.exe>.
"@
    exit 1
}
if (-not (Test-Path $dist)) {
    Write-Error "dist/ not found at $dist. Run Ladruno_scripts\build.bat first."
    exit 1
}
if (-not (Test-Path $iss)) {
    Write-Error "installer.iss not found at $iss"
    exit 1
}
if (-not (Test-Path $banner)) {
    Write-Error "banner_ASCII.txt not found at $banner"
    exit 1
}
if (-not (Test-Path $out)) {
    New-Item -ItemType Directory -Path $out | Out-Null
}

# ---- ensure banner has UTF-8 BOM so Inno's TStringList decodes it --------
# banner_ASCII.txt is the source-of-truth (no BOM). Inno's TStringList only
# auto-detects encoding when a BOM is present; without one it falls back to
# the system ANSI codepage and the box-drawing glyphs render as garbage.
# We write a parallel BOM-prefixed copy that the .iss bundles via dontcopy.
$bannerText = [System.IO.File]::ReadAllText($banner, [System.Text.Encoding]::UTF8)
$utf8WithBom = New-Object System.Text.UTF8Encoding($true)
[System.IO.File]::WriteAllText($bannerBom, $bannerText, $utf8WithBom)

# ---- run iscc ------------------------------------------------------------
$env:LADRUNO_VERSION = $Version
$env:LADRUNO_DIST    = $dist
$env:LADRUNO_OUT     = $out

Write-Host "Compiling installer ..." -ForegroundColor Yellow
& $Iscc /Q $iss
$exitCode = $LASTEXITCODE
if ($exitCode -ne 0) {
    Write-Error "iscc.exe exited with code $exitCode"
    exit $exitCode
}

# ---- summary -------------------------------------------------------------
$produced = Join-Path $out "Ladruno_OpenSees_${Version}_setup.exe"
if (-not (Test-Path $produced)) {
    Write-Warning "iscc reported success but $produced was not produced. Inspect output above."
    exit 1
}

$mb = [math]::Round((Get-Item $produced).Length / 1MB, 1)
Write-Host ""
Write-Host "Built: $produced ($mb MB)" -ForegroundColor Green
Write-Host ""
Write-Host "Distribute this single .exe to recipients." -ForegroundColor Yellow
Write-Host "For automation / power users, install.ps1 + .zip in $out remain available."
