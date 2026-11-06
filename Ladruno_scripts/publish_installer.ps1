<#
.SYNOPSIS
    Attach the locally built Inno installer to an existing GitHub release (WP-166).

.DESCRIPTION
    Refuses unless dist\ is a FULL 5-target build (same rule as build.bat
    `installer`: a partial dist\ ships a mixed install). Uploads the newest
    Ladruno_files\*_setup.exe with `gh release upload`. The release (tag
    ladruno-v*) must already exist; this script never creates tags or releases.

.EXAMPLE
    Ladruno_scripts\build.bat clean installer
    powershell -ExecutionPolicy Bypass -File Ladruno_scripts\publish_installer.ps1 -Tag ladruno-v2026.10
#>
param(
    [Parameter(Mandatory = $true)][string]$Tag,
    [string]$Repo = "nmorabowen/OpenSees",
    [switch]$Clobber
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$dist = Join-Path $root "dist"

# The 5 build targets' artifacts (see AGENTS.md "Artifacts").
$required = @(
    "bin\OpenSees.exe", "bin\OpenSeesSP.exe", "bin\OpenSeesMP.exe",
    "bin\opensees.pyd", "openseesmp\openseesmp.pyd"
)
$missing = @($required | Where-Object { -not (Test-Path (Join-Path $dist $_)) })
if ($missing.Count -gt 0) {
    Write-Error ("dist\ is not a full 5-target build; missing: " + ($missing -join ", ") +
                 ". Run Ladruno_scripts\build.bat installer (no explicit targets).")
    exit 1
}

$setup = Get-ChildItem (Join-Path $root "Ladruno_files") -Filter "*_setup.exe" -ErrorAction SilentlyContinue |
         Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $setup) { Write-Error "No Ladruno_files\*_setup.exe; run build.bat installer first."; exit 1 }

# The installer must not be older than the binaries it wraps.
$newest = $required | ForEach-Object { (Get-Item (Join-Path $dist $_)).LastWriteTime } |
          Sort-Object -Descending | Select-Object -First 1
if ($setup.LastWriteTime -lt $newest) {
    Write-Error "$($setup.Name) is older than the dist\ binaries; rebuild the installer."; exit 1
}

gh release view $Tag -R $Repo *> $null
if ($LASTEXITCODE -ne 0) {
    Write-Error "Release $Tag does not exist on $Repo (the maintainer creates it by tagging)."; exit 1
}

Write-Host "Uploading $($setup.Name) to $Repo $Tag"
$ghArgs = @("release", "upload", $Tag, $setup.FullName, "-R", $Repo)
if ($Clobber) { $ghArgs += "--clobber" }
& gh @ghArgs
exit $LASTEXITCODE
