# ADR-86b T8 -- re-run ADR-90 GATE U's COARSE legs on a capped engine.
#
# Four concurrent legs, one 45-minute chunk (the GATE U note's own infrastructure
# rule: budget any re-run in <= 50-minute chunks, or run in the foreground -- both
# of that campaign's launches were killed from outside at ~1 h).
#
#   h1.0_e0.6944  x  -maxSubsteps {1000, 10000}
#   h1.0_e0.60    x  -maxSubsteps {1000, 10000}
#
# LADRUNO_DIST_BIN pins the engine to THIS worktree's dist/bin and
# LADRUNO_A2_EXPECT_BUILD pins its git hash, so a stale .pyd fails loudly rather
# than producing a number attributed to the wrong build.
param(
  [string]$Root = "C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\agent-aeef311decf9449d1",
  [string]$OutRoot,
  [double]$Wall = 2700,
  [int[]]$Caps = @(1000, 10000)
)

$ErrorActionPreference = "Stop"
if (-not $OutRoot) { $OutRoot = Join-Path $Root "Ladruno_files\testbed\hypo_bearing\adr86b_t8" }
New-Item -ItemType Directory -Force -Path $OutRoot | Out-Null

# The engine's OWN compiled-in stamp, not `git rev-parse HEAD`: a commit that
# touches only tests/ moves HEAD without rebuilding, and the driver refuses to
# run when the two disagree. Read it from the binary under test.
$env:PYTHONPATH = Join-Path $Root "dist\bin"
$env:LADRUNO_OPENSEES_QUIET = "1"
$hash = (& $(if (Test-Path "C:\Users\nmb\AppData\Local\Programs\Python\Python312\python.exe")
              { "C:\Users\nmb\AppData\Local\Programs\Python\Python312\python.exe" }
              else { (Get-Command python3.12).Source }) `
         -c "import opensees as o;print(o.ladrunoBuild())" 2>$null | Select-Object -Last 1).Trim()
$py   = "C:\Users\nmb\AppData\Local\Programs\Python\Python312\python.exe"
if (-not (Test-Path $py)) { $py = (Get-Command python3.12).Source }
$drv  = Join-Path $Root "Ladruno_files\testbed\hypo_bearing\sanisand_tau0_band.py"

"engine hash pinned to $hash"
"driver          $drv"
"out root        $OutRoot"

$jobs = @()
foreach ($leg in @("h1.0_e0.6944", "h1.0_e0.60")) {
  foreach ($cap in $Caps) {
    $dir = Join-Path $OutRoot ("{0}_cap{1}" -f $leg, $cap)
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $log = Join-Path $dir "run.log"
    $env:LADRUNO_DIST_BIN = Join-Path $Root "dist\bin"
    $env:LADRUNO_A2_EXPECT_BUILD = $hash
    $env:OMP_NUM_THREADS = "4"
    $args = @($drv, "--out", $dir, "--legs", $leg, "--wall", $Wall,
              "--maxsubsteps", $cap)
    $p = Start-Process -FilePath $py -ArgumentList $args -PassThru `
           -RedirectStandardOutput $log -RedirectStandardError ($log + ".err") `
           -WorkingDirectory $Root -WindowStyle Hidden
    "launched $leg cap=$cap  pid=$($p.Id)  -> $dir"
    $jobs += $p
  }
}
"all launched at $(Get-Date -Format o); pids: $($jobs.Id -join ',')"
