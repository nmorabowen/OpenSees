#==============================================================================
#  OpenSees Windows 11 -- Step 1: Install Dependencies
#
#  Run this script ONCE on a clean Windows 11 machine (elevated PowerShell).
#  It installs every tool needed to compile OpenSees from source.
#
#  Usage (elevated / "Run as Administrator"):
#    powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\1_install_dependencies.ps1
#
#  Optional flags:
#    -SkipInnoSetup         Skip Inno Setup (only needed for installer packaging)
#    -SkipVerification      Skip the post-install verification step
#    -DryRun                Print what would be installed without installing
#==============================================================================
[CmdletBinding()]
param(
    [switch]$SkipInnoSetup,
    [switch]$SkipVerification,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$script:DryRun = $DryRun.IsPresent

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
function Write-Banner {
    param([string]$Message)
    Write-Host ""
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "  $Message" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
}

function Write-Step {
    param([string]$Message)
    Write-Host ""
    Write-Host "[STEP] $Message" -ForegroundColor Green
}

function Write-Warn {
    param([string]$Message)
    Write-Host "[WARN] $Message" -ForegroundColor Yellow
}

function Write-Fail {
    param([string]$Message)
    Write-Host "[FAIL] $Message" -ForegroundColor Red
}

function Test-Elevated {
    $identity  = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = New-Object Security.Principal.WindowsPrincipal($identity)
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Invoke-WingetInstall {
    param(
        [string]$Id,
        [string]$DisplayName,
        [string]$Override = ""
    )

    Write-Step "Installing $DisplayName ($Id)"

    if ($script:DryRun) {
        Write-Host "  [DRY RUN] winget install --id $Id"
        if ($Override) { Write-Host "  [DRY RUN]   --override: $Override" }
        return
    }

    $wingetArgs = @(
        "install", "--id", $Id, "-e",
        "--accept-source-agreements", "--accept-package-agreements"
    )
    if ($Override) {
        $wingetArgs += "--override"
        $wingetArgs += $Override
    }

    & winget @wingetArgs
    if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne -1978335189) {
        # -1978335189 = "already installed" -- not an error
        Write-Warn "$DisplayName install returned exit code $LASTEXITCODE (may already be installed)"
    }
}

# --------------------------------------------------------------------------
# Pre-flight
# --------------------------------------------------------------------------
Write-Banner "OpenSees Windows 11 -- Dependency Installer"

if (-not (Test-Elevated)) {
    Write-Fail "This script must run in an elevated (Administrator) PowerShell."
    Write-Host "  Right-click PowerShell -> Run as Administrator, then re-run." -ForegroundColor Yellow
    exit 1
}

# Check that winget is available
try {
    $null = Get-Command winget -ErrorAction Stop
} catch {
    Write-Fail "winget is not available on this system."
    Write-Host "  Install App Installer from the Microsoft Store, or upgrade Windows." -ForegroundColor Yellow
    exit 1
}

if ($DryRun) {
    Write-Host ""
    Write-Host "  *** DRY RUN MODE -- nothing will be installed ***" -ForegroundColor Magenta
}

# --------------------------------------------------------------------------
# 1. Git
# --------------------------------------------------------------------------
Invoke-WingetInstall -Id "Git.Git" -DisplayName "Git"

# --------------------------------------------------------------------------
# 2. Visual Studio 2022 with C++ workload
# --------------------------------------------------------------------------
Invoke-WingetInstall `
    -Id "Microsoft.VisualStudio.2022.Community" `
    -DisplayName "Visual Studio 2022 Community (C++ workload)" `
    -Override "--wait --quiet --add Microsoft.VisualStudio.Workload.NativeDesktop --includeRecommended"

# --------------------------------------------------------------------------
# 3. Intel oneAPI Base Toolkit (MKL)
# --------------------------------------------------------------------------
Invoke-WingetInstall -Id "Intel.OneAPI.BaseToolkit" -DisplayName "Intel oneAPI Base Toolkit (MKL)"

# --------------------------------------------------------------------------
# 4. Intel oneAPI HPC Toolkit (ifx, Intel MPI)
# --------------------------------------------------------------------------
Invoke-WingetInstall -Id "Intel.OneAPI.HPCToolkit" -DisplayName "Intel oneAPI HPC Toolkit (ifx, MPI)"

# --------------------------------------------------------------------------
# 5. CMake
# --------------------------------------------------------------------------
Invoke-WingetInstall -Id "Kitware.CMake" -DisplayName "CMake"

# --------------------------------------------------------------------------
# 6. Ninja
# --------------------------------------------------------------------------
Invoke-WingetInstall -Id "Ninja-build.Ninja" -DisplayName "Ninja"

# --------------------------------------------------------------------------
# 7. Python 3.11 (x64) -- pinned for OpenSeesPy ABI compatibility
# --------------------------------------------------------------------------
Invoke-WingetInstall -Id "Python.Python.3.11" -DisplayName "Python 3.11 (x64)"

# --------------------------------------------------------------------------
# 8. Inno Setup (optional -- for installer packaging only)
# --------------------------------------------------------------------------
if (-not $SkipInnoSetup) {
    Invoke-WingetInstall -Id "JRSoftware.InnoSetup" -DisplayName "Inno Setup 6 (optional, for installer)"
} else {
    Write-Step "Skipping Inno Setup (-SkipInnoSetup)"
}

# --------------------------------------------------------------------------
# Refresh PATH so newly installed tools are visible
# --------------------------------------------------------------------------
if (-not $DryRun) {
    Write-Step "Refreshing PATH from registry"
    $machinePath = [Environment]::GetEnvironmentVariable("Path", "Machine")
    $userPath    = [Environment]::GetEnvironmentVariable("Path", "User")
    $env:Path    = "$machinePath;$userPath"
}

# --------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------
if ($SkipVerification -or $DryRun) {
    Write-Step "Skipping verification"
} else {
    Write-Banner "Verifying installations"

    $allOk = $true

    # --- Git ---
    Write-Step "Checking Git"
    try   { & git --version } catch { Write-Fail "Git not found"; $allOk = $false }

    # --- CMake ---
    Write-Step "Checking CMake"
    try   { & cmake --version } catch { Write-Fail "CMake not found"; $allOk = $false }

    # --- Ninja ---
    Write-Step "Checking Ninja"
    try   { & ninja --version } catch { Write-Fail "Ninja not found"; $allOk = $false }

    # --- Python 3.11 ---
    Write-Step "Checking Python 3.11"
    $pyFound = $false
    foreach ($pyCandidate in @("py -3.11 --version", "python3.11 --version", "python --version")) {
        try {
            $pyOut = Invoke-Expression $pyCandidate 2>&1
            if ($pyOut -match "3\.11") { $pyFound = $true; Write-Host "  $pyOut"; break }
        } catch {}
    }
    if (-not $pyFound) { Write-Fail "Python 3.11 not found"; $allOk = $false }

    # --- Visual Studio 2022 ---
    Write-Step "Checking Visual Studio 2022"
    $vswhereExe = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
    if (Test-Path $vswhereExe) {
        $vsPath = & $vswhereExe -latest -products * -property installationPath 2>&1
        if ($vsPath) { Write-Host "  Found: $vsPath" }
        else         { Write-Fail "vswhere found but no VS installation detected"; $allOk = $false }
    } else {
        Write-Fail "vswhere.exe not found -- Visual Studio may not be installed"
        $allOk = $false
    }

    # --- Intel oneAPI (ifx, MKL, mpiexec) ---
    Write-Step "Checking Intel oneAPI (ifx and mpiexec)"
    $oneapiSetvars = "${env:ProgramFiles(x86)}\Intel\oneAPI\setvars.bat"
    if (Test-Path $oneapiSetvars) {
        Write-Host "  setvars.bat found: $oneapiSetvars"
        $checkCmdParts = @()
        if ($vsPath) {
            $checkCmdParts += "set `"VS2022INSTALLDIR=$($vsPath.Trim())`""
        }
        $checkCmdParts += @(
            "call `"$oneapiSetvars`" intel64 >nul 2>&1",
            "where ifx",
            "where mpiexec"
        )
        $checkCmd = $checkCmdParts -join " && "
        $result = cmd /c $checkCmd 2>&1
        if ($LASTEXITCODE -eq 0) {
            Write-Host "  ifx and mpiexec found after oneAPI init"
        } else {
            Write-Warn "oneAPI setvars.bat exists but ifx/mpiexec not found after init"
            Write-Warn "  You may need to reboot or re-run the oneAPI installer"
            $allOk = $false
        }
    } else {
        Write-Fail "Intel oneAPI setvars.bat not found"
        $allOk = $false
    }

    # --- Inno Setup ---
    if (-not $SkipInnoSetup) {
        Write-Step "Checking Inno Setup"
        $isccPaths = @(
            "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe",
            "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
            "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
        )
        $isccFound = $false
        foreach ($p in $isccPaths) {
            if (Test-Path $p) { Write-Host "  Found: $p"; $isccFound = $true; break }
        }
        if (-not $isccFound) {
            Write-Warn "Inno Setup not found (optional, only needed for installer packaging)"
        }
    }

    # --- Summary ---
    Write-Host ""
    if ($allOk) {
        Write-Banner "All mandatory dependencies verified successfully"
    } else {
        Write-Banner "Some dependencies are missing -- review the warnings above"
        Write-Host ""
        Write-Host "  If tools were just installed, try closing and reopening PowerShell," -ForegroundColor Yellow
        Write-Host "  or reboot, then re-run this script to verify." -ForegroundColor Yellow
    }
}

# --------------------------------------------------------------------------
# Next steps
# --------------------------------------------------------------------------
Write-Host ""
Write-Host "Next step:" -ForegroundColor Green
Write-Host "  Run SCRIPTS\2_fetch_source.ps1 to clone OpenSees and the build harness." -ForegroundColor White
Write-Host ""
