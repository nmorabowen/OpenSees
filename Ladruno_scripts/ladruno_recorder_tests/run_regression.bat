@echo off
REM ==========================================================================
REM  Ladruno recorder regression — full gate battery vs frozen MPCO.
REM
REM  Runs every model/check pair against the CURRENT worktree build. The model
REM  scripts run under the BUILD python (matches the opensees.pyd ABI); the
REM  check scripts run under the venv python (h5py + numpy). TMP points at the
REM  freshly-built dist\bin so the pyd and its MKL DLLs are co-located -- that
REM  is why ALL gates run here, not just the 3 that hardcode a separate MKL dir.
REM
REM  Usage (from any cwd):
REM    Ladruno_scripts\ladruno_recorder_tests\run_regression.bat
REM
REM  Prereq: build OpenSeesPy first (Ladruno_scripts\build.bat OpenSeesPy) so
REM  %ROOT%\dist\bin\opensees.pyd is current.
REM ==========================================================================
setlocal enabledelayedexpansion

REM ----- resolve paths (this script lives in Ladruno_scripts\ladruno_recorder_tests) -----
set "TESTS=%~dp0"
if "%TESTS:~-1%"=="\" set "TESTS=%TESTS:~0,-1%"
pushd "%TESTS%\..\.."
set "ROOT=%CD%"
popd

set "DISTBIN=%ROOT%\dist\bin"
set "BUILDPY=C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
set "VENVPY=C:\Users\nmora\venv\opensees_venv\Scripts\python.exe"
set "OUT=%TESTS%\_devrun"

if not exist "%DISTBIN%\opensees.pyd" (
    echo ERROR: %DISTBIN%\opensees.pyd not found. Build first:
    echo     Ladruno_scripts\build.bat OpenSeesPy
    exit /b 1
)
if not exist "%OUT%" mkdir "%OUT%"
cd /d "%TESTS%"

set /a FAILS=0

echo ==========================================================
echo  Ladruno recorder regression
echo  build : %DISTBIN%\opensees.pyd
echo  out   : %OUT%
echo ==========================================================

REM ---- existing gates -------------------------------------------------------
call :gate2 "NODAL PARITY"        parity_model.py          parity_check.py          ref.mpco        test.ladruno
call :gate2 "ELEMENT PARITY quad" element_parity_model.py  element_parity_check.py  eref.mpco       etest.ladruno
call :checkonly "ELEMENT PARITY beam" element_parity_check.py eref_beam.mpco etest_beam.ladruno
call :gate2 "MULTI-STAGE"         multistage_model.py      multistage_check.py      ms_ref.mpco     ms_test.ladruno
call :gate1 "STANDARD QUAD"       standard_quad_model.py   standard_quad_check.py
call :gate1 "NODAL ENVELOPE"      envelope_model.py        envelope_check.py
call :gate1 "ELEMENT ENVELOPE"    element_envelope_model.py element_envelope_check.py
call :gateA "LOCAL AXES 3D"       localaxes_model.py       localaxes_check.py       "%OUT%" localaxes.ladruno
call :gateA "LOCAL AXES 2D"       localaxes_model_2d.py    localaxes_check.py       "%OUT%" localaxes_2d.ladruno
call :gateA "ENERGY BALANCE"      energy_model.py          energy_check.py          "%OUT%\energy.ladruno" "%OUT%\energy_sidecar.txt"
call :gate1 "BEZIER TRI6"         bezier_model.py          bezier_check.py

REM ---- new element-type coverage -------------------------------------------
call :gate2 "TRUSS/ZEROLENGTH FORCE" truss_force_model.py element_parity_check.py  tref.mpco       ttest.ladruno
call :checkonly "TRUSS/ZEROLENGTH NODAL" parity_check.py  tref.mpco       ttest.ladruno
call :gate2 "FRAME3D localForce"  frame3d_model.py         element_parity_check.py  fref.mpco       ftest.ladruno
call :checkonly "FRAME3D NODAL"   parity_check.py          fref.mpco       ftest.ladruno
call :gate2 "SHELL forces"        shell_model.py           element_parity_check.py  shref.mpco      shtest.ladruno
call :checkonly "SHELL NODAL"     parity_check.py          shref.mpco      shtest.ladruno
call :gate2 "EIGEN modes"         eigen_model.py           eigen_check.py           eig_ref.mpco    eig_test.ladruno
call :gateA "PRECISION f32"       precision_model.py       precision_check.py       "%OUT%\prec_f64.ladruno" "%OUT%\prec_f32.ladruno"

echo ==========================================================
if %FAILS% GTR 0 (
    echo  REGRESSION: %FAILS% gate failures
    exit /b 1
)
echo  REGRESSION: ALL GATES PASSED
echo  note: mp_parallel gate is openseesmp-only; run it separately
exit /b 0

REM ==========================================================================
REM  :gate2  label model check refName newName   -> model writes ref+new, 2-arg check
REM ==========================================================================
:gate2
echo.
echo --- %~1 ---
"%BUILDPY%" "%TESTS%\%~2" "%DISTBIN%" "%OUT%"
if errorlevel 1 (echo [FAIL] %~1 : model runner failed & set /a FAILS+=1 & goto :eof)
"%VENVPY%" "%TESTS%\%~3" "%OUT%\%~4" "%OUT%\%~5"
if errorlevel 1 (echo [FAIL] %~1 & set /a FAILS+=1) else (echo [PASS] %~1)
goto :eof

REM  :gate1  label model check   -> model writes, check takes OUT dir only
:gate1
echo.
echo --- %~1 ---
"%BUILDPY%" "%TESTS%\%~2" "%DISTBIN%" "%OUT%"
if errorlevel 1 (echo [FAIL] %~1 : model runner failed & set /a FAILS+=1 & goto :eof)
"%VENVPY%" "%TESTS%\%~3" "%OUT%"
if errorlevel 1 (echo [FAIL] %~1 & set /a FAILS+=1) else (echo [PASS] %~1)
goto :eof

REM  :gateA  label model check arg2 arg3  -> model writes, check takes 2 explicit args
:gateA
echo.
echo --- %~1 ---
"%BUILDPY%" "%TESTS%\%~2" "%DISTBIN%" "%OUT%"
if errorlevel 1 (echo [FAIL] %~1 : model runner failed & set /a FAILS+=1 & goto :eof)
"%VENVPY%" "%TESTS%\%~3" %~4 %~5
if errorlevel 1 (echo [FAIL] %~1 & set /a FAILS+=1) else (echo [PASS] %~1)
goto :eof

REM  :checkonly label check refName newName  -> reuse files from a prior model run
:checkonly
echo.
echo --- %~1 ---
"%VENVPY%" "%TESTS%\%~2" "%OUT%\%~3" "%OUT%\%~4"
if errorlevel 1 (echo [FAIL] %~1 & set /a FAILS+=1) else (echo [PASS] %~1)
goto :eof
