<#
.SYNOPSIS
    Wire the local OpenSeesPy build (dist/bin/opensees.pyd + bundled MKL
    DLLs) into a Python virtualenv so `import opensees` works without
    fiddling with sys.path.

.DESCRIPTION
    Creates the venv at -VenvPath if it doesn't already exist, then drops
    a .pth file in its site-packages pointing at dist/bin/. That puts
    opensees.pyd on sys.path automatically every time the venv is
    activated. The MKL/iomp DLLs co-located with opensees.pyd are found
    by Windows because they live in the same directory as the .pyd.

    No pip install / wheel build needed; this is a developer convenience
    so edits to OpenSeesPy show up the next time you re-run the
    interpreter — no re-install step.

.PARAMETER VenvPath
    Where to create / find the venv. Default: <project-root>\opensees_venv

.PARAMETER Python
    Python interpreter to use as the venv base. Default: the Python 3.11
    we built OpenSeesPy against. Mismatching this will produce
    `ImportError: DLL load failed` on import.

.PARAMETER Force
    Recreate the venv from scratch (deletes -VenvPath first).

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File Ladruno_scripts\wire_pyenv.ps1
    powershell -ExecutionPolicy Bypass -File Ladruno_scripts\wire_pyenv.ps1 -VenvPath D:\envs\ops -Force
#>
param(
    [string]$VenvPath = "",
    [string]$Python   = "C:\Users\nmora\AppData\Local\Programs\Python\Python311\python.exe",
    [switch]$Force
)

$ErrorActionPreference = "Stop"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

$scriptDir = $PSScriptRoot
$root      = Split-Path -Parent $scriptDir
$distBin   = Join-Path $root "dist\bin"
$logDir    = Join-Path $root "Ladruno_files"

if ([string]::IsNullOrWhiteSpace($VenvPath)) {
    $VenvPath = Join-Path $root "opensees_venv"
}

Write-Host "Ladruno OpenSeesPy venv wiring" -ForegroundColor Cyan
Write-Host "  python    : $Python"
Write-Host "  venv      : $VenvPath"
Write-Host "  pyd dir   : $distBin"
Write-Host ""

# ---------- preflight checks ----------------------------------------------
if (-not (Test-Path $Python)) {
    Write-Error "Python interpreter not found at: $Python`nOverride with -Python <path-to-python.exe>."
    exit 1
}
if (-not (Test-Path "$distBin\opensees.pyd")) {
    Write-Error "$distBin\opensees.pyd not found. Run Ladruno_scripts\build.bat OpenSeesPy first."
    exit 1
}

# Make sure the chosen Python matches the .pyd's ABI (cpython 3.11 == cp311)
$pyVersion = & $Python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
if ($pyVersion -ne "3.11") {
    Write-Warning "Selected Python is $pyVersion; opensees.pyd was built for 3.11."
    Write-Warning "It will likely fail to import. Pass -Python <path-to-python311.exe> to override."
}

# ---------- create / refresh venv -----------------------------------------
if ($Force -and (Test-Path $VenvPath)) {
    Write-Host "Removing existing venv (Force) ..."
    Remove-Item -Recurse -Force $VenvPath
}

if (-not (Test-Path $VenvPath)) {
    Write-Host "Creating venv ..."
    & $Python -m venv $VenvPath
    if ($LASTEXITCODE -ne 0) { Write-Error "venv creation failed."; exit 1 }
} else {
    Write-Host "Reusing existing venv."
}

# ---------- locate site-packages ------------------------------------------
$venvPython = Join-Path $VenvPath "Scripts\python.exe"
if (-not (Test-Path $venvPython)) {
    Write-Error "Expected $venvPython to exist after venv creation."
    exit 1
}
$sitePackages = & $venvPython -c "import sysconfig; print(sysconfig.get_paths()['purelib'])"
Write-Host "  site-pkg  : $sitePackages"

# ---------- write the .pth file -------------------------------------------
# A .pth file is read by site.py at interpreter startup; each non-blank
# non-comment line is added to sys.path. We use a single-line pth so
# `import opensees` finds opensees.pyd in dist/bin/.
$pthPath = Join-Path $sitePackages "ladruno_opensees.pth"
Set-Content -Path $pthPath -Value $distBin -Encoding ASCII
Write-Host "Wrote .pth file: $pthPath"
Write-Host "                 -> $distBin"

# ---------- smoke test ----------------------------------------------------
Write-Host ""
Write-Host "Smoke test: import opensees in the venv ..." -NoNewline
$smokeOut = & $venvPython -c @"
import sys, os
import opensees as ops
ops.wipe()
ops.model('BasicBuilder', '-ndm', 2, '-ndf', 2)
ops.node(1, 0.0, 0.0); ops.node(2, 1.0, 0.0)
ops.fix(1, 1, 1); ops.fix(2, 0, 1)
ops.uniaxialMaterial('Elastic', 1, 1000.0)
ops.element('truss', 1, 1, 2, 1.0, 1)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1); ops.load(2, 10.0, 0.0)
ops.constraints('Plain'); ops.numberer('Plain')
ops.system('BandSPD'); ops.test('NormDispIncr', 1.0e-9, 10)
ops.algorithm('Linear'); ops.integrator('LoadControl', 1.0)
ops.analysis('Static'); ops.analyze(1)
print(f'ux={ops.nodeDisp(2, 1):.6f}')
"@ 2>&1
if ($LASTEXITCODE -ne 0) {
    Write-Host " FAILED" -ForegroundColor Red
    Write-Host $smokeOut
    exit 1
}
Write-Host " ok ($($smokeOut.Trim()))" -ForegroundColor Green

# ---------- log to Ladruno_files/ -----------------------------------------
if (-not (Test-Path $logDir)) { New-Item -ItemType Directory -Path $logDir | Out-Null }
$log = Join-Path $logDir "wire_pyenv.log"
@"
Wired at:    $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")
Venv:        $VenvPath
Python:      $Python ($pyVersion)
.pth file:   $pthPath
.pth target: $distBin
Smoke test:  $($smokeOut.Trim())
"@ | Set-Content -Path $log -Encoding UTF8

# ---------- usage hint ----------------------------------------------------
Write-Host ""
Write-Host "Activate the venv:" -ForegroundColor Yellow
Write-Host "  cmd:        $VenvPath\Scripts\activate.bat"
Write-Host "  powershell: & '$VenvPath\Scripts\Activate.ps1'"
Write-Host ""
Write-Host "Then in Python:"
Write-Host "  import opensees as ops"
Write-Host ""
Write-Host "Log: $log" -ForegroundColor DarkGray
