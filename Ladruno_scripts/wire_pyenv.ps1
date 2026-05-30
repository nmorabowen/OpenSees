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
    Python interpreter to use as the venv base. Default: the Python 3.12
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
    [string]$Python   = "C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe",
    [switch]$Force
)

$ErrorActionPreference = "Stop"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

$scriptDir = $PSScriptRoot
$root      = Split-Path -Parent $scriptDir
$distBin   = Join-Path $root "dist\bin"
$distMp    = Join-Path $root "dist\openseesmp"
$wirePy    = Join-Path $scriptDir "wire_venv_pth.py"
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

# Make sure the chosen Python matches the .pyd's ABI (cpython 3.12 == cp312)
$pyVersion = & $Python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
if ($pyVersion -ne "3.12") {
    Write-Warning "Selected Python is $pyVersion; opensees.pyd was built for 3.12."
    Write-Warning "It will likely fail to import. Pass -Python <path-to-python312.exe> to override."
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

# ---------- write the .pth via the shared wirer ---------------------------
# Delegate to wire_venv_pth.py (the SAME helper the Inno installer runs) so
# the standalone and installer paths can never diverge. It drops a
# _ladruno_opensees_boot.py + one-line .pth that put BOTH dist\bin and
# dist\openseesmp on sys.path, register their bundled-DLL dirs via
# os.add_dll_directory (process-local; no global PATH change), and alias
# openseespy -> our sequential opensees (skipped under MPI).
if (-not (Test-Path $wirePy)) {
    Write-Error "wire_venv_pth.py not found at $wirePy"
    exit 1
}
$mpArg = if (Test-Path "$distMp\openseesmp.pyd") { $distMp } else { "" }
& $venvPython $wirePy $distBin $mpArg
$wireExit = $LASTEXITCODE
$pthPath = Join-Path $sitePackages "ladruno_opensees.pth"
if ($wireExit -eq 0) {
    Write-Host "Wrote .pth + _ladruno_opensees_boot.py in $sitePackages"
    Write-Host "                 -> import opensees   from $distBin"
    if ($mpArg) { Write-Host "                 -> import openseesmp from $distMp" }
    else        { Write-Host "                 (openseesmp not built; run build.bat OpenSeesPyMP to add it)" }
    Write-Host "                 + openseespy/.opensees aliased (skipped under MPI)"
} elseif ($wireExit -eq 3) {
    Write-Warning "venv Python is not 3.12; the .pyd modules will fail to import."
} else {
    Write-Error "wire_venv_pth.py exited $wireExit"
    exit 1
}

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
Write-Host "  import opensees as ops          # sequential"
Write-Host ""
Write-Host "For MPI-parallel Python, launch with the bundled mpiexec:" -ForegroundColor Yellow
Write-Host "  $distMp\mpiexec.exe -n 4 $venvPython your_script.py"
Write-Host "  (your_script.py does:  import openseesmp as ops)"
Write-Host ""
Write-Host "Log: $log" -ForegroundColor DarkGray
