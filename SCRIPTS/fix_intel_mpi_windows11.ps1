[CmdletBinding()]
param(
    [int]$Ranks = 2,
    [switch]$SetWinRmAutomatic
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Write-Log {
    param(
        [Parameter(Mandatory = $true)][string]$Message,
        [ValidateSet("INFO", "WARN", "ERROR")][string]$Level = "INFO"
    )
    $ts = (Get-Date).ToString("yyyy-MM-dd HH:mm:ss")
    $line = "[$ts] [$Level] $Message"
    switch ($Level) {
        "WARN" { Write-Host $line -ForegroundColor Yellow }
        "ERROR" { Write-Host $line -ForegroundColor Red }
        default { Write-Host $line }
    }
}

function Test-IsAdmin {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = New-Object Security.Principal.WindowsPrincipal($identity)
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Invoke-CmdChecked {
    param(
        [Parameter(Mandatory = $true)][string]$Command
    )
    Write-Log ">> $Command"
    cmd /c $Command
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $Command"
    }
}

function Try-Cmd {
    param(
        [Parameter(Mandatory = $true)][string]$Command,
        [string]$FailureMessage = ""
    )
    Write-Log ">> $Command"
    cmd /c $Command
    if ($LASTEXITCODE -ne 0) {
        if ($FailureMessage) {
            Write-Log $FailureMessage "WARN"
        } else {
            Write-Log "Command failed with exit code ${LASTEXITCODE}: $Command" "WARN"
        }
    }
}

if (-not (Test-IsAdmin)) {
    throw "This script must be run as Administrator."
}

$oneApiRoot = ${env:ONEAPI_ROOT}
if (-not $oneApiRoot) {
    $defaultOneApiRoot = "C:\Program Files (x86)\Intel\oneAPI"
    if (Test-Path -Path $defaultOneApiRoot -PathType Container) {
        $oneApiRoot = $defaultOneApiRoot
    }
}

if (-not $oneApiRoot) {
    throw "ONEAPI_ROOT was not set and default oneAPI root was not found."
}

$setvarsBat = Join-Path $oneApiRoot "setvars.bat"
if (-not (Test-Path -Path $setvarsBat -PathType Leaf)) {
    throw "Could not find setvars.bat at '$setvarsBat'"
}

$hydraServiceExe = Join-Path $oneApiRoot "mpi\latest\bin\hydra_service.exe"
if (-not (Test-Path -Path $hydraServiceExe -PathType Leaf)) {
    throw "Could not find hydra_service.exe at '$hydraServiceExe'"
}

Write-Log "Using ONEAPI_ROOT: $oneApiRoot"

Write-Log "Configuring WinRM (required by Intel MPI powershell launcher on Windows client OS)."
Try-Cmd -Command "winrm quickconfig -quiet" -FailureMessage "winrm quickconfig reported an error. Check WinRM logs and rerun."

if ($SetWinRmAutomatic) {
    Try-Cmd -Command "sc config WinRM start= auto" -FailureMessage "Could not set WinRM start type to auto."
}

Try-Cmd -Command "net start WinRM" -FailureMessage "WinRM may already be running or failed to start."

Write-Log "Installing/starting Intel MPI Hydra service."
Try-Cmd -Command "`"$hydraServiceExe`" -remove" -FailureMessage "Hydra service remove failed (continuing)."
Invoke-CmdChecked -Command "`"$hydraServiceExe`" -install"
Invoke-CmdChecked -Command "`"$hydraServiceExe`" -start"

Write-Log "Verifying MPI runtime with mpiexec."
$mpiCheck = "call `"$setvarsBat`" intel64 >nul 2>&1 && mpiexec -n $Ranks hostname"
Invoke-CmdChecked -Command $mpiCheck

Write-Log "Intel MPI bootstrap check passed. OpenSeesSP/OpenSeesMP should now be runnable." "INFO"
