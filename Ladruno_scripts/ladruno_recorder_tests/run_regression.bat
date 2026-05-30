@echo off
rem Regression gates: nodal parity, element parity, multi-stage — all vs frozen mpco.
setlocal
set WT=C:\Users\nmora\Github\OpenSees_Compile\mpco-ladruno-wt
set BUILDPY=C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe
set VENVPY=C:\Users\nmora\venv\opensees_venv\Scripts\python.exe
set TESTS=%WT%\Ladruno_scripts\ladruno_recorder_tests
set TMP=%TESTS%\_devrun
if not exist "%TMP%" mkdir "%TMP%"
copy /Y "%WT%\build\build\Release\OpenSeesPy.dll" "%TMP%\opensees.pyd" >nul
cd /d "%TESTS%"

echo === NODAL PARITY ===
"%BUILDPY%" "%TESTS%\parity_model.py" "%TMP%" "%TMP%"
if errorlevel 1 (echo NODAL MODEL FAILED & exit /b 1)
"%VENVPY%" "%TESTS%\parity_check.py" "%TMP%\ref.mpco" "%TMP%\test.ladruno"
if errorlevel 1 (echo NODAL PARITY FAILED & exit /b 1)

echo === ELEMENT PARITY ===
"%BUILDPY%" "%TESTS%\element_parity_model.py" "%TMP%" "%TMP%"
if errorlevel 1 (echo ELEM MODEL FAILED & exit /b 1)
"%VENVPY%" "%TESTS%\element_parity_check.py" "%TMP%\eref.mpco" "%TMP%\etest.ladruno"
if errorlevel 1 (echo ELEM PARITY FAILED & exit /b 1)

echo === MULTI-STAGE ===
"%BUILDPY%" "%TESTS%\multistage_model.py" "%TMP%" "%TMP%"
if errorlevel 1 (echo MS MODEL FAILED & exit /b 1)
"%VENVPY%" "%TESTS%\multistage_check.py" "%TMP%\ms_ref.mpco" "%TMP%\ms_test.ladruno"
if errorlevel 1 (echo MS CHECK FAILED & exit /b 1)

echo === ALL REGRESSION GATES PASSED ===
exit /b 0
