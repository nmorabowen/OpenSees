@echo off
REM ADR-95: relaunch the long legs OUTSIDE the Claude session (they die when it restarts). Run from a plain cmd window.
cd /d "%~dp0"
powershell -NoProfile -Command "Get-CimInstance Win32_Process | Where-Object { $_.Name -like 'python*' -and $_.CommandLine -like '*path_diag*' } | ForEach-Object { Stop-Process -Id $_.ProcessId -Force }"
set SC=C:\Users\nmb\AppData\Local\Temp\claude\C--Users-nmb-Documents-Github-OpenSees--claude-worktrees-prandtl-reissner-bezier-954a90\5408bcef-bf62-4967-964e-a179fa61676a\scratchpad
set P4=%SC%\dist_p4\bin
set FX=%SC%\dist_fixed\bin
start "adr95 h20 long"  /min cmd /c "set ADR95_DIST=%P4%&& python3.12 -u quad_path_diag.py --elem h20uri --h0 1.0 --cond --cond-at 2e-3 --branch --tmax 14400 --suffix _p4long > p4_h20uri_long.log 2>&1 & echo P1_DONE exit=%%ERRORLEVEL%% >> p4_h20uri_long.log"
start "adr95 tet10 long" /min cmd /c "set ADR95_DIST=%P4%&& python3.12 -u tet_path_diag.py --elem tet10 --branch --cond-at 5e-3 --cond-every 50 --sfrac 0.15 --budget 200 --tmax 14400 --suffix _p4long > p4_tet10_long.log 2>&1 & echo P1_DONE exit=%%ERRORLEVEL%% >> p4_tet10_long.log"
start "adr95 bezstd long" /min cmd /c "set ADR95_DIST=%P4%&& python3.12 -u tet_path_diag.py --elem beziertet10 --branch --cond-at 5e-3 --cond-every 50 --sfrac 0.15 --budget 200 --tmax 14400 --suffix _p4long > p4_bezstd_long.log 2>&1 & echo P1_DONE exit=%%ERRORLEVEL%% >> p4_bezstd_long.log"
start "adr95 sy2"  /min cmd /c "set ADR95_DIST=%FX%&& python3.12 -u quad_path_diag.py --elem h20uri --h0 1.0 --sy 2.0 --cond --cond-at 1e-3 --branch --tmax 7200 --suffix _p2sy2 > p2_h20uri_sy2.log 2>&1 & echo P1_DONE exit=%%ERRORLEVEL%% >> p2_h20uri_sy2.log"
start "adr95 sy20" /min cmd /c "set ADR95_DIST=%FX%&& python3.12 -u quad_path_diag.py --elem h20uri --h0 1.0 --sy 20.0 --cond --cond-at 1e-3 --branch --tmax 7200 --suffix _p2sy20 > p2_h20uri_sy20.log 2>&1 & echo P1_DONE exit=%%ERRORLEVEL%% >> p2_h20uri_sy20.log"
echo launched 5 legs; logs: p4_h20uri_long p4_tet10_long p4_bezstd_long p2_h20uri_sy2 p2_h20uri_sy20 (look for "MODE =" / P1_DONE)
