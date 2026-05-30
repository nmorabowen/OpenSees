<#
.SYNOPSIS
  Run the Zone-A test bed against a built OpenSees Python module without the
  PYTHONPATH dance. Fork-local convenience (not upstreamed).

.DESCRIPTION
  Points PYTHONPATH/PATH at a build directory containing opensees.pyd (+ its MKL
  DLLs), then runs `pytest -m zone_a` (or any marker) in tests/ with py 3.12.

.EXAMPLE
  ./Ladruno_scripts/run_zone_a.ps1
  ./Ladruno_scripts/run_zone_a.ps1 -Markers t0 -DistBin C:\...\OpenSees\dist\bin
  ./Ladruno_scripts/run_zone_a.ps1 -- -k bezier        # extra args after -- go to pytest

.NOTES
  In a git worktree the build usually lives in the MAIN checkout, so pass
  -DistBin <main-checkout>\dist\bin explicitly.
#>
param(
    [string]$DistBin = (Join-Path $PSScriptRoot '..\dist\bin'),
    [string]$Markers = 'zone_a',
    [string]$Python  = '3.12',
    [Parameter(ValueFromRemainingArguments = $true)] $Extra
)
$ErrorActionPreference = 'Stop'

if (-not (Test-Path (Join-Path $DistBin 'opensees.pyd'))) {
    throw "No opensees.pyd under '$DistBin'. Build first, or pass -DistBin <path-with-opensees.pyd> (e.g. the main checkout's dist\bin when running from a worktree)."
}
$DistBin = (Resolve-Path $DistBin).Path
$tests   = (Resolve-Path (Join-Path $PSScriptRoot '..\tests')).Path

$env:PYTHONPATH = $DistBin
$env:PATH       = "$DistBin;$env:PATH"

Write-Host "opensees.pyd : $DistBin"
Write-Host "tests        : $tests"
Write-Host "markers      : $Markers`n"

Push-Location $tests
try {
    & py "-$Python" -m pytest -v -m $Markers @Extra
    exit $LASTEXITCODE
} finally {
    Pop-Location
}
