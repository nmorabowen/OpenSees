# Run the parallel lumped-SMS shared-node validation: np=1 (reference) vs np=2 (split).
# Uses the openseesmp.pyd Python module (CentralDifferenceSMS is registered in the
# interpreter layer, NOT the legacy Tcl integrator parser).
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = (Resolve-Path "$here\..\..").Path
$dist = Join-Path $root 'dist\openseesmp'   # openseesmp.pyd lives here, not dist\bin
$py   = 'C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe'
$mpiexec = 'C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin\mpiexec.exe'

Set-Location $here
Remove-Item -ErrorAction SilentlyContinue tip_np1.out, tip_np2.out

Write-Host "=== np=1 (reference: whole bar on one rank) ==="
& $mpiexec -n 1 $py sms_bar_mp.py $dist
Write-Host "=== np=2 (split at shared node 11) ==="
& $mpiexec -n 2 $py sms_bar_mp.py $dist

Write-Host "=== compare ==="
& $py compare.py
exit $LASTEXITCODE
