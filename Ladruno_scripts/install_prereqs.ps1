<#
.SYNOPSIS
    install_prereqs.ps1 -- install the OpenSees (Ladruno fork) build toolchain.

.DESCRIPTION
    One-shot, idempotent installer for everything build.bat needs on a fresh
    Windows machine. Mirrors the prereq list documented in build.bat:

        1. Visual Studio 2022 Community  (C++ "Desktop development" workload)
        2. Intel oneAPI Base Toolkit     (MKL: BLAS / LAPACK / ScaLAPACK)
        3. Intel oneAPI HPC Toolkit      (ifx Fortran compiler + Intel MPI)
        4. CMake, Ninja, Python 3.12     (build drivers)
        5. Conan 2.x                     (pip --user; HDF5 / TCL / ZLIB / Eigen3)

    MUMPS is NOT installed here -- build.bat downloads and builds it on first run.

    Most installs need administrator rights, so winget raises a UAC prompt per
    package. Either accept each prompt, or run this whole script from an
    elevated PowerShell to be asked once.

.NOTES
    After this completes, open a fresh shell and run:
        Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP

.PARAMETER WhatIf
    Print the actions without installing anything.
#>
[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [switch]$IncludeVSWorkloadRecommended = $true
)

$ErrorActionPreference = 'Stop'

function Test-WingetPackage {
    param([string]$Id)
    # `winget list --id <id>` exits 0 and echoes a row when present.
    $out = winget list --id $Id --exact --source winget 2>$null
    return ($LASTEXITCODE -eq 0 -and ($out -join "`n") -match [regex]::Escape($Id))
}

function Install-WingetPackage {
    param(
        [string]$Id,
        [string]$Name,
        [string[]]$ExtraArgs = @()
    )
    Write-Host "`n=== $Name ($Id) ===" -ForegroundColor Cyan
    if (Test-WingetPackage -Id $Id) {
        Write-Host "  already installed -- skipping" -ForegroundColor Green
        return
    }
    if (-not $PSCmdlet.ShouldProcess($Name, "winget install")) { return }

    $args = @(
        'install', '--id', $Id, '--exact', '--source', 'winget',
        '--accept-package-agreements', '--accept-source-agreements',
        '--disable-interactivity'
    ) + $ExtraArgs

    Write-Host "  installing (a UAC prompt may appear -- approve it) ..." -ForegroundColor Yellow
    winget @args
    # winget returns 0 on success; -1978335189 == "no applicable upgrade / already installed"
    if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne -1978335189) {
        throw "winget install $Id failed (exit $LASTEXITCODE)"
    }
    Write-Host "  done" -ForegroundColor Green
}

# ---- 0. Preconditions ----------------------------------------------------
if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
    throw "winget not found. Install 'App Installer' from the Microsoft Store, then re-run."
}
Write-Host "winget: $((Get-Command winget).Source)  ($(winget --version))"

# ---- 1. Visual Studio 2022 Community + C++ Desktop workload --------------
# The bare winget package installs the shell only; --override hands the VS
# bootstrapper the workload + recommended components it actually needs to
# provide cl.exe, link.exe, and the Windows SDK.
$vsOverride = '--quiet --wait --norestart --add Microsoft.VisualStudio.Workload.NativeDesktop'
if ($IncludeVSWorkloadRecommended) { $vsOverride += ' --includeRecommended' }
Install-WingetPackage -Id 'Microsoft.VisualStudio.2022.Community' `
    -Name 'Visual Studio 2022 Community (C++ Desktop workload)' `
    -ExtraArgs @('--override', $vsOverride)

# ---- 2. Intel oneAPI Base Toolkit (MKL) ---------------------------------
Install-WingetPackage -Id 'Intel.OneAPI.BaseToolkit' `
    -Name 'Intel oneAPI Base Toolkit (MKL)'

# ---- 3. Intel oneAPI HPC Toolkit (ifx Fortran + Intel MPI) ---------------
# Installed AFTER Base: the HPC toolkit layers onto the Base oneAPI root.
Install-WingetPackage -Id 'Intel.OneAPI.HPCToolkit' `
    -Name 'Intel oneAPI HPC Toolkit (ifx + MPI)'

# ---- 4. Build drivers ----------------------------------------------------
Install-WingetPackage -Id 'Kitware.CMake'        -Name 'CMake'
Install-WingetPackage -Id 'Ninja-build.Ninja'    -Name 'Ninja'
Install-WingetPackage -Id 'Python.Python.3.12'   -Name 'Python 3.12' `
    -ExtraArgs @('--scope', 'user')

# ---- 5. Conan 2.x (pip --user) ------------------------------------------
# setup_env.bat puts %APPDATA%\Python\Python312\Scripts on PATH, which is
# exactly where `pip install --user conan` drops conan.exe.
Write-Host "`n=== Conan 2.x (pip --user) ===" -ForegroundColor Cyan
$py = @(
    "$env:LOCALAPPDATA\Programs\Python\Python312\python.exe",
    "$env:ProgramFiles\Python312\python.exe"
) | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $py) {
    $cmd = Get-Command python -ErrorAction SilentlyContinue
    if ($cmd) { $py = $cmd.Source }
}
if (-not $py) {
    Write-Warning "Python 3.12 not found on PATH yet (new shell may be needed). Run 'pip install --user conan' manually."
} elseif ($PSCmdlet.ShouldProcess('conan', 'pip install --user')) {
    & $py -m pip install --user --upgrade pip conan
    if ($LASTEXITCODE -ne 0) { throw "pip install conan failed (exit $LASTEXITCODE)" }
    Write-Host "  done" -ForegroundColor Green
}

# ---- 6. Summary ----------------------------------------------------------
Write-Host "`n=== Prerequisites installed ===" -ForegroundColor Cyan
Write-Host @"
Next steps (open a FRESH shell so PATH/installer changes take effect):

  cd $($PSScriptRoot | Split-Path)
  Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP

build.bat loads the toolchain via setup_env.bat, builds MUMPS once
(~15-20 min), runs Conan, configures CMake, and builds the targets.

If a Visual Studio or oneAPI UAC prompt was dismissed, re-run this script --
it skips what is already installed and retries the rest.
"@ -ForegroundColor Gray
