@echo off
setlocal enabledelayedexpansion
REM ==========================================================================
REM  mutation_build.bat -- build an ADR-87 D2 MUTANT OpenSeesPy on Windows.
REM
REM  WHY THIS EXISTS
REM    .github/workflows/ladruno_mutation.yml builds mutants on Ubuntu with a
REM    bare `cmake -S . -B build/zero -DCMAKE_TOOLCHAIN_FILE=...`. That recipe
REM    does NOT work on this machine for two reasons:
REM      1. build.bat configures via `cmake --preset conan-release`, and a
REM         preset PINS its own binaryDir -- you cannot redirect it with -B, so
REM         a mutant would overwrite the normal build tree.
REM      2. add_compile_definitions(LADRUNO_MUTATE_<FAMILY>=...) is GLOBAL, so
REM         toggling the define in the shared tree recompiles all ~1969 objects
REM         -- twice, once each way. A separate build dir costs one build and
REM         leaves the normal tree untouched.
REM    So this mirrors build.bat's Windows configure flags (ifx, MUMPS, forced
REM    response files, CMS off) into build\mut_<FAMILY>_<MODE>\.
REM
REM  Usage (from a fresh cmd.exe, at the fork root):
REM    Ladruno_scripts\mutation_build.bat SHELLMOD ZERO
REM    Ladruno_scripts\mutation_build.bat CONTINUUM SCALE
REM
REM  Output: a self-contained module dir  dist\mut_<FAMILY>_<MODE>\  holding the
REM  mutant opensees.pyd beside a copy of dist\bin's runtime DLLs, so it can be
REM  handed straight to mutation_gate.py --module-dir without disturbing the
REM  normal dist\bin.
REM
REM  Then:
REM    python Ladruno_scripts\mutation_gate.py run --family <FAMILY> ^
REM        --expect <FAMILY>=<MODE> --module-dir dist\mut_<FAMILY>_<MODE> ^
REM        --out artifacts\<mode>.json
REM ==========================================================================

if "%~1"=="" (echo ERROR: usage: mutation_build.bat ^<FAMILY^> ^<MODE^>  e.g. SHELLMOD ZERO & exit /b 1)
if "%~2"=="" (echo ERROR: usage: mutation_build.bat ^<FAMILY^> ^<MODE^>  e.g. SHELLMOD ZERO & exit /b 1)
set "FAMILY=%~1"
set "MODE=%~2"

set "SCRIPT_DIR=%~dp0"
if "%SCRIPT_DIR:~-1%"=="\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
pushd "%SCRIPT_DIR%\.."
set "ROOT=%CD%"
popd

set "SRC=%ROOT%"
set "MUMPS_INSTALL=%ROOT%\mumps-install"
set "MUT_BUILD=%ROOT%\build\mut_%FAMILY%_%MODE%"
set "MUT_DIST=%ROOT%\dist\mut_%FAMILY%_%MODE%"
set "TOOLCHAIN=%ROOT%\build\build\Release\generators\conan_toolchain.cmake"

REM The mutant reuses the NORMAL build's conan toolchain, so a normal build must
REM have been configured at least once. Fail loud rather than half-configure.
if not exist "%TOOLCHAIN%" (
    echo ERROR: no conan toolchain at "%TOOLCHAIN%".
    echo Run `Ladruno_scripts\build.bat` once first -- the mutant build reuses it.
    exit /b 1
)

REM ----- Locate Python 3.12 (same probe order as build.bat) ----------------
if defined PYEXE if exist "%PYEXE%" goto :pyexe_ok
set "PYEXE="
for %%P in (
    "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
    "C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
    "%ProgramFiles%\Python312\python.exe"
) do if not defined PYEXE if exist "%%~P" set "PYEXE=%%~P"
if not defined PYEXE for /f "delims=" %%P in ('where python 2^>nul') do if not defined PYEXE set "PYEXE=%%P"
if not defined PYEXE (echo ERROR: could not locate a Python 3.12 interpreter; set PYEXE explicitly & exit /b 1)
:pyexe_ok
echo === Using Python: %PYEXE% ===

REM ----- Toolchain (same guard idiom as build.bat) -------------------------
if defined SKIP_SETUP_ENV (
    echo === Toolchain reload skipped ^(SKIP_SETUP_ENV set^) ===
) else if defined LADRUNO_TOOLCHAIN_LOADED (
    echo === Toolchain already loaded in this shell; skipping setup_env.bat ===
) else (
    echo === Loading toolchain ===
    call "%SCRIPT_DIR%\setup_env.bat" || (echo setup_env.bat failed & exit /b 1)
)

echo.
echo === MUTANT configure: %FAMILY%=%MODE% -^> %MUT_BUILD% ===
pushd "%SRC%"
cmake -S "%SRC%" -B "%MUT_BUILD%" -G Ninja ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_TOOLCHAIN_FILE="%TOOLCHAIN%" ^
    -DCMAKE_Fortran_COMPILER=ifx ^
    -DPython_EXECUTABLE="%PYEXE%" ^
    -DMUMPS_DIR="%MUMPS_INSTALL%/lib" ^
    -DMUMPS_INCLUDE_DIR="%MUMPS_INSTALL%/include" ^
    -DLADRUNO_CMS=OFF -DLADRUNO_CMS_BUILD_TESTS=OFF ^
    -DCMAKE_NINJA_FORCE_RESPONSE_FILE=ON ^
    -DLADRUNO_MUTATE_FAMILY=%FAMILY% ^
    -DLADRUNO_MUTATE_MODE=%MODE%
set "RC=%errorlevel%"
popd
if not "%RC%"=="0" (echo MUTANT CMake configure failed & exit /b 1)

echo.
echo === MUTANT build: OpenSeesPy ===
cmake --build "%MUT_BUILD%" --target OpenSeesPy
if not "%errorlevel%"=="0" (echo MUTANT build failed & exit /b 1)

REM ----- Assemble a self-contained module dir ------------------------------
REM Copy the NORMAL dist\bin (for the co-located MKL/Tcl/etc runtime DLLs),
REM then overwrite opensees.pyd with the mutant. Copying rather than adding the
REM mutant to dist\bin keeps the normal module pristine -- a mutant left in
REM dist\bin is exactly the silent-wrong-binary failure mutation_gate.py's
REM --expect check exists to catch.
echo.
echo === Assembling %MUT_DIST% ===
if not exist "%ROOT%\dist\bin\opensees.pyd" (
    echo ERROR: no normal dist\bin\opensees.pyd to source runtime DLLs from.
    echo Run `Ladruno_scripts\build.bat OpenSees OpenSeesPy` first.
    exit /b 1
)
if exist "%MUT_DIST%" rmdir /s /q "%MUT_DIST%"
mkdir "%MUT_DIST%"
xcopy /q /y "%ROOT%\dist\bin\*.dll" "%MUT_DIST%\" >nul 2>nul
copy /y "%ROOT%\dist\bin\opensees.pyd" "%MUT_DIST%\opensees.pyd" >nul

if not exist "%MUT_BUILD%\OpenSeesPy.dll" (
    echo ERROR: mutant OpenSeesPy.dll not found in "%MUT_BUILD%"
    exit /b 1
)
copy /y "%MUT_BUILD%\OpenSeesPy.dll" "%MUT_DIST%\opensees.pyd" >nul

echo.
echo === Done. Mutant module: %MUT_DIST%\opensees.pyd ===
echo Verify it reports itself before trusting any score:
echo   set PYTHONPATH=%MUT_DIST%
echo   "%PYEXE%" -c "import opensees as o; print(o.ladrunoMutation())"
echo Expect: %FAMILY%=%MODE%
echo.
endlocal
