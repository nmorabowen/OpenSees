#==============================================================================
#  OpenSees Windows 11 — Step 2: Fetch Source & Build Harness
#
#  Clones an OpenSees source tree, the Windows build-harness files, vcpkg,
#  and MUMPS, then copies the harness into the source tree so it is ready
#  to compile.
#
#  Usage (normal PowerShell, no elevation needed):
#    powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\2_fetch_source.ps1
#
#  Common overrides:
#    -OpenSeesRepo "https://github.com/YourFork/OpenSees.git"
#    -OpenSeesBranch "my-feature"
#    -WorkDir "D:\builds\opensees"
#    -DryRun
#==============================================================================
[CmdletBinding()]
param(
    # The OpenSees source repository to compile.
    [string]$OpenSeesRepo = "https://github.com/OpenSees/OpenSees.git",

    # Branch / tag / commit to check out (empty = default branch).
    [string]$OpenSeesBranch = "",

    # Root working directory.  The source tree will be created inside it.
    [string]$WorkDir = "C:\work",

    # Name of the folder that holds the OpenSees source.
    [string]$SourceDirName = "OpenSees-src",

    # Windows build-harness repository (contains CMakeLists.txt, scripts, etc.).
    [string]$HarnessRepo = "https://github.com/nmorabowen/opensees_windows_compilation_files.git",

    # MUMPS repository and pinned commit.
    [string]$MumpsRepo = "https://github.com/OpenSees/mumps.git",
    [string]$MumpsCommit = "ec5f340",

    # Print commands without executing.
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

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

function Invoke-Git {
    param([string[]]$Arguments)

    $argText = $Arguments -join " "
    Write-Host "  >> git $argText"

    if ($DryRun) {
        Write-Host "  [DRY RUN] skipped"
        return
    }

    & git @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "git $argText failed with exit code $LASTEXITCODE"
    }
}

# --------------------------------------------------------------------------
# Resolve paths
# --------------------------------------------------------------------------
$SourceRoot  = Join-Path $WorkDir $SourceDirName
$HarnessDir  = Join-Path $WorkDir "opensees_windows_compilation_files"
$VcpkgDir    = Join-Path $SourceRoot "third_party\vcpkg"
$MumpsDir    = Join-Path $SourceRoot "third_party\mumps"

Write-Banner "OpenSees Windows 11 — Fetch Source & Build Harness"
Write-Host ""
Write-Host "  OpenSees repo  : $OpenSeesRepo"
Write-Host "  Branch/tag     : $(if ($OpenSeesBranch) { $OpenSeesBranch } else { '(default)' })"
Write-Host "  Working dir    : $WorkDir"
Write-Host "  Source dir     : $SourceRoot"
Write-Host "  Harness repo   : $HarnessRepo"
Write-Host ""

if ($DryRun) {
    Write-Host "  *** DRY RUN — nothing will be written to disk ***" -ForegroundColor Magenta
}

# --------------------------------------------------------------------------
# 1. Create working directory
# --------------------------------------------------------------------------
Write-Step "Creating working directory: $WorkDir"
if (-not $DryRun) {
    New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
}

# --------------------------------------------------------------------------
# 2. Clone OpenSees source
# --------------------------------------------------------------------------
Write-Step "Cloning OpenSees source into $SourceRoot"
if (Test-Path $SourceRoot) {
    Write-Host "  Directory already exists — skipping clone."
} else {
    $cloneArgs = @("clone", $OpenSeesRepo, $SourceRoot)
    if ($OpenSeesBranch) {
        $cloneArgs = @("clone", "--branch", $OpenSeesBranch, $OpenSeesRepo, $SourceRoot)
    }
    Invoke-Git -Arguments $cloneArgs
}

# If a specific branch was requested and the dir already existed, check it out.
if ($OpenSeesBranch -and (Test-Path $SourceRoot) -and (-not $DryRun)) {
    Push-Location $SourceRoot
    try {
        $currentBranch = git rev-parse --abbrev-ref HEAD 2>&1
        if ($currentBranch -ne $OpenSeesBranch) {
            Write-Host "  Switching to branch: $OpenSeesBranch"
            Invoke-Git -Arguments @("-C", $SourceRoot, "checkout", $OpenSeesBranch)
        }
    } finally { Pop-Location }
}

# --------------------------------------------------------------------------
# 3. Clone vcpkg
# --------------------------------------------------------------------------
Write-Step "Cloning vcpkg into $VcpkgDir"
if (Test-Path $VcpkgDir) {
    Write-Host "  Directory already exists — skipping clone."
} else {
    if (-not $DryRun) {
        New-Item -ItemType Directory -Force -Path (Join-Path $SourceRoot "third_party") | Out-Null
    }
    Invoke-Git -Arguments @("clone", "https://github.com/microsoft/vcpkg.git", $VcpkgDir)
}

# --------------------------------------------------------------------------
# 4. Clone MUMPS (pinned commit)
# --------------------------------------------------------------------------
Write-Step "Cloning MUMPS into $MumpsDir (commit $MumpsCommit)"
if (Test-Path $MumpsDir) {
    Write-Host "  Directory already exists — skipping clone."
} else {
    Invoke-Git -Arguments @("clone", $MumpsRepo, $MumpsDir)
    if ($MumpsCommit -and (-not $DryRun)) {
        Invoke-Git -Arguments @("-C", $MumpsDir, "checkout", $MumpsCommit)
    }
}

# --------------------------------------------------------------------------
# 5. Clone the Windows build-harness repository
# --------------------------------------------------------------------------
Write-Step "Cloning Windows build harness into $HarnessDir"
if (Test-Path $HarnessDir) {
    Write-Host "  Directory already exists — pulling latest."
    if (-not $DryRun) {
        Invoke-Git -Arguments @("-C", $HarnessDir, "pull", "--ff-only")
    }
} else {
    Invoke-Git -Arguments @("clone", $HarnessRepo, $HarnessDir)
}

# --------------------------------------------------------------------------
# 6. Copy build-harness files into the OpenSees source tree
# --------------------------------------------------------------------------
Write-Step "Copying build-harness files into $SourceRoot"

if (-not $DryRun) {
    # Ensure destination directories exist.
    New-Item -ItemType Directory -Force -Path (Join-Path $SourceRoot "SCRIPTS")       | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $SourceRoot "cmake\cmake")   | Out-Null
}

# Map: source (in harness repo) -> destination (in OpenSees source tree).
$filesToCopy = @(
    @{ Src = "CMakeLists.txt";                               Dst = "CMakeLists.txt" },
    @{ Src = "vcpkg.json";                                   Dst = "vcpkg.json" },
    @{ Src = "cmake\cmake\OpenSeesDependencies.cmake";       Dst = "cmake\cmake\OpenSeesDependencies.cmake" },
    @{ Src = "cmake\cmake\OpenSeesDependenciesWin.cmake";    Dst = "cmake\cmake\OpenSeesDependenciesWin.cmake" },
    @{ Src = "SCRIPTS\build_windows11_full.ps1";             Dst = "SCRIPTS\build_windows11_full.ps1" },
    @{ Src = "SCRIPTS\1_install_dependencies.ps1";           Dst = "SCRIPTS\1_install_dependencies.ps1" },
    @{ Src = "SCRIPTS\2_fetch_source.ps1";                   Dst = "SCRIPTS\2_fetch_source.ps1" },
    @{ Src = "SCRIPTS\3_build.ps1";                          Dst = "SCRIPTS\3_build.ps1" },
    @{ Src = "SCRIPTS\4_package.ps1";                        Dst = "SCRIPTS\4_package.ps1" },
    @{ Src = "SCRIPTS\init_oneapi_windows11.cmd";            Dst = "SCRIPTS\init_oneapi_windows11.cmd" },
    @{ Src = "SCRIPTS\fix_intel_mpi_windows11.ps1";          Dst = "SCRIPTS\fix_intel_mpi_windows11.ps1" },
    @{ Src = "SCRIPTS\create_el_ladruno_installer.ps1";      Dst = "SCRIPTS\create_el_ladruno_installer.ps1" }
)

foreach ($entry in $filesToCopy) {
    $srcPath = Join-Path $HarnessDir $entry.Src
    $dstPath = Join-Path $SourceRoot $entry.Dst

    if (Test-Path $srcPath) {
        Write-Host "  $($entry.Src)  ->  $($entry.Dst)"
        if (-not $DryRun) {
            Copy-Item -Path $srcPath -Destination $dstPath -Force
        }
    } else {
        Write-Host "  [skip] $($entry.Src) not found in harness repo" -ForegroundColor Yellow
    }
}

# --------------------------------------------------------------------------
# 7. Also copy source-code patches if present
# --------------------------------------------------------------------------
$sourcePatches = @(
    @{ Src = "SRC\element\UWelements\Tcl_generateInterfacePoints.cpp";
       Dst = "SRC\element\UWelements\Tcl_generateInterfacePoints.cpp" },
    @{ Src = "SRC\modelbuilder\tcl\myCommands.cpp";
       Dst = "SRC\modelbuilder\tcl\myCommands.cpp" }
)

$hasSrcPatches = $false
foreach ($entry in $sourcePatches) {
    $srcPath = Join-Path $HarnessDir $entry.Src
    if (Test-Path $srcPath) { $hasSrcPatches = $true; break }
}

if ($hasSrcPatches) {
    Write-Step "Copying source-code patches (MSVC linker fixes)"
    foreach ($entry in $sourcePatches) {
        $srcPath = Join-Path $HarnessDir $entry.Src
        $dstPath = Join-Path $SourceRoot $entry.Dst
        if (Test-Path $srcPath) {
            Write-Host "  $($entry.Src)"
            if (-not $DryRun) {
                Copy-Item -Path $srcPath -Destination $dstPath -Force
            }
        }
    }
}

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
Write-Banner "Source tree is ready"
Write-Host ""
Write-Host "  Source root     : $SourceRoot" -ForegroundColor White
Write-Host "  vcpkg           : $VcpkgDir" -ForegroundColor White
Write-Host "  MUMPS           : $MumpsDir" -ForegroundColor White
Write-Host "  Build harness   : copied from $HarnessDir" -ForegroundColor White
Write-Host ""
Write-Host "Next step:" -ForegroundColor Green
Write-Host "  cd $SourceRoot" -ForegroundColor White
Write-Host "  powershell -NoProfile -ExecutionPolicy Bypass -File SCRIPTS\3_build.ps1" -ForegroundColor White
Write-Host ""
