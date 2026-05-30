@echo off
rem ---------------------------------------------------------------------------
rem Energy-balance (D8) gate: copy the fresh dev build, run the model with the
rem BUILD python (no boot .pth), then validate with the venv python (has h5py).
rem Usage: run_energy.bat
rem ---------------------------------------------------------------------------
setlocal
set WT=C:\Users\nmora\Github\OpenSees_Compile\mpco-ladruno-wt
set BUILDPY=C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe
set VENVPY=C:\Users\nmora\venv\opensees_venv\Scripts\python.exe
set TESTS=%WT%\Ladruno_scripts\ladruno_recorder_tests
set TMP=%TESTS%\_devrun
if not exist "%TMP%" mkdir "%TMP%"
copy /Y "%WT%\build\build\Release\OpenSeesPy.dll" "%TMP%\opensees.pyd" >nul
if errorlevel 1 (echo COPY FAILED & exit /b 1)

echo === [1/2] energy_model.py (BUILD python) ===
"%BUILDPY%" "%TESTS%\energy_model.py" "%TMP%" "%TMP%"
if errorlevel 1 (echo MODEL FAILED & exit /b 1)

echo === [2/2] energy_check.py (venv python) ===
cd /d "%TESTS%"
"%VENVPY%" "%TESTS%\energy_check.py" "%TMP%\energy.ladruno" "%TMP%\energy_sidecar.txt"
exit /b %errorlevel%
