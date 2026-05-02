@echo off
setlocal enabledelayedexpansion
REM ==========================================================================
REM  build.bat — drive the OpenSees build end to end
REM
REM  Prereqs (one-time, on this machine):
REM    1. Visual Studio 2022 Community (with the C++ workload + Windows SDK)
REM    2. Intel oneAPI Base + HPC toolkit (compiler, MKL, MPI)
REM    3. CMake 3.23+, Ninja, Python 3.12, Git
REM    4. Conan 2.x  (pip install --user conan)
REM    5. MUMPS 5.5.1 source archive at mumps-archive\mumps_src.tar.gz
REM       (downloaded once from https://mumps-solver.org/MUMPS_5.5.1.tar.gz)
REM
REM  Usage (from a fresh cmd.exe):
REM    cd <fork-clone-root>                   (e.g. C:\Users\nmora\Github\OpenSees)
REM    Ladruno_scripts\build.bat              -> build everything
REM    Ladruno_scripts\build.bat clean        -> wipe build/install/dist, rebuild
REM    Ladruno_scripts\build.bat rebuild      -> wipe build/, rebuild
REM    Ladruno_scripts\build.bat <target>     -> rebuild only that target
REM                                              (OpenSees, OpenSeesSP, OpenSeesMP, OpenSeesPy)
REM
REM  Incremental rebuilds after source edits in OpenSees\: just re-run this
REM  script. CMake + Ninja handle dependency tracking.
REM ==========================================================================

REM This script lives in Ladruno_scripts/. Resolve the project root one level up.
set "SCRIPT_DIR=%~dp0"
if "%SCRIPT_DIR:~-1%"=="\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
pushd "%SCRIPT_DIR%\.."
set "ROOT=%CD%"
popd

REM After the fork-restructure, Ladruno_scripts/ lives inside OpenSees/, so
REM ROOT itself is the OpenSees source root.
set "SRC=%ROOT%"
set "BUILD=%ROOT%\build"
set "BUILD_DIR=%BUILD%\build\Release"
set "INSTALL=%ROOT%\install"
set "DIST=%ROOT%\dist"
set "PROFILE=%SCRIPT_DIR%\opensees-msvc-static.profile"
set "PYEXE=C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"

set "MUMPS_SRC=%ROOT%\mumps-src"
set "MUMPS_BUILD=%ROOT%\mumps-build"
set "MUMPS_INSTALL=%ROOT%\mumps-install"
set "MUMPS_ARCHIVE=%ROOT%\mumps-archive"

set "MKL_BIN=C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin"
set "ICOMP_BIN=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"
set "IMPI_BIN=C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin"

REM ----- Argument parsing --------------------------------------------------
set "MODE=%1"
if /i "%MODE%"=="clean" (
    echo Wiping build/, install/, dist/ ...
    rmdir /s /q "%BUILD%"   2>nul
    rmdir /s /q "%INSTALL%" 2>nul
    rmdir /s /q "%DIST%"    2>nul
    set "MODE="
) else if /i "%MODE%"=="rebuild" (
    echo Wiping build/ ...
    rmdir /s /q "%BUILD%" 2>nul
    set "MODE="
)

set "TARGETS=OpenSees OpenSeesSP OpenSeesMP OpenSeesPy"
if not "%MODE%"=="" set "TARGETS=%MODE%"

REM ----- Load toolchain ----------------------------------------------------
echo === Loading toolchain ===
call "%SCRIPT_DIR%\setup_env.bat" || (echo setup_env.bat failed & exit /b 1)

REM ----- 1. MUMPS: build once if not already installed ---------------------
if not exist "%MUMPS_INSTALL%\lib\dmumps.lib" (
    echo.
    echo === Step 1: Building MUMPS 5.5.1 (one-time, ~15-20 min) ===
    if not exist "%MUMPS_ARCHIVE%\mumps_src.tar.gz" (
        echo.
        echo MUMPS source archive not found. Downloading from upstream:
        if not exist "%MUMPS_ARCHIVE%" mkdir "%MUMPS_ARCHIVE%"
        curl -L -o "%MUMPS_ARCHIVE%\mumps_src.tar.gz" ^
             "https://mumps-solver.org/MUMPS_5.5.1.tar.gz"
        if errorlevel 1 (
            echo ERROR: failed to download MUMPS source. Manually place the
            echo archive at %MUMPS_ARCHIVE%\mumps_src.tar.gz
            exit /b 1
        )
    )
    if not exist "%MUMPS_SRC%" (
        echo Cloning scivision/mumps build system...
        git clone --depth 1 --branch v5.5.1.5 ^
            https://github.com/scivision/mumps.git "%MUMPS_SRC%"
        if errorlevel 1 (echo MUMPS clone failed & exit /b 1)
    )
    if exist "%MUMPS_BUILD%" rmdir /s /q "%MUMPS_BUILD%"
    REM intsize64=OFF: 32-bit Fortran INTEGER, matches Intel MPI's default
    REM LP64 layer and OpenSees' LP64 ScaLAPACK linkage. Setting this ON
    REM produces "Instance Error 1 in DMUMPS_F77" at runtime under Intel MPI.
    cmake -S "%MUMPS_SRC%" -B "%MUMPS_BUILD%" -G Ninja ^
        -DCMAKE_BUILD_TYPE=Release ^
        -DCMAKE_C_COMPILER=cl -DCMAKE_Fortran_COMPILER=ifx ^
        -Dintsize64=OFF -Darith=d ^
        -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF ^
        -Dlocal="%MUMPS_ARCHIVE%" ^
        -DCMAKE_INSTALL_PREFIX="%MUMPS_INSTALL%" ^
        -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded
    if errorlevel 1 (echo MUMPS configure failed & exit /b 1)
    cmake --build "%MUMPS_BUILD%" -j 8
    if errorlevel 1 (echo MUMPS build failed & exit /b 1)
    cmake --install "%MUMPS_BUILD%"
    if errorlevel 1 (echo MUMPS install failed & exit /b 1)
) else (
    echo === Step 1: MUMPS already built at %MUMPS_INSTALL% ===
)

REM ----- 2. Conan: fetch / verify deps -------------------------------------
echo.
echo === Step 2: Conan install (HDF5, TCL, ZLIB, Eigen3) ===
conan install "%SRC%\conanfile.py" ^
    --output-folder="%BUILD%" ^
    --profile:host="%PROFILE%" ^
    --profile:build="%PROFILE%" ^
    --build=missing
if errorlevel 1 (echo Conan install failed & exit /b 1)

REM ----- 3. CMake configure ------------------------------------------------
echo.
echo === Step 3: CMake configure ===
pushd "%SRC%"
cmake --preset conan-release ^
    -DCMAKE_Fortran_COMPILER=ifx ^
    -DCMAKE_INSTALL_PREFIX="%INSTALL%" ^
    -DPython_EXECUTABLE="%PYEXE%" ^
    -DMUMPS_DIR="%MUMPS_INSTALL%/lib" ^
    -DMUMPS_INCLUDE_DIR="%MUMPS_INSTALL%/include"
set "RC=%errorlevel%"
popd
if not "%RC%"=="0" (echo CMake configure failed & exit /b 1)

REM ----- 4. Build each target ----------------------------------------------
echo.
echo === Step 4: Building targets: %TARGETS% ===
for %%T in (%TARGETS%) do (
    echo.
    echo --- Building %%T ---
    cmake --build "%BUILD_DIR%" --target %%T -j 8
    if errorlevel 1 (echo Build of %%T failed & exit /b 1)
)

REM ----- 5. Curate dist ----------------------------------------------------
echo.
echo === Step 5: Copying binaries to dist\ ===
if not exist "%DIST%\bin"         mkdir "%DIST%\bin"
if not exist "%DIST%\lib\tcl8.6"  mkdir "%DIST%\lib\tcl8.6"

REM Executables and the OpenSeesPy module
for %%F in (OpenSees.exe OpenSeesSP.exe OpenSeesMP.exe) do (
    if exist "%BUILD_DIR%\%%F" (
        echo   copying %%F
        copy /y "%BUILD_DIR%\%%F" "%DIST%\bin\" >nul
    )
)
REM OpenSeesPy is built as OpenSeesPy.dll; Python expects opensees.pyd
if exist "%BUILD_DIR%\OpenSeesPy.dll" (
    echo   copying OpenSeesPy.dll -^> opensees.pyd
    copy /y "%BUILD_DIR%\OpenSeesPy.dll" "%DIST%\bin\opensees.pyd" >nul
)

REM Intel MKL runtime DLLs. mkl_intel_thread + mkl_core are link-time
REM dependencies; mkl_def / mkl_avx* / mkl_mc3 are CPU kernel DLLs that MKL
REM dlopens at runtime (FATAL "mkl_def.2.dll not found" otherwise).
REM mkl_scalapack_lp64 + mkl_blacs_intelmpi_lp64 are needed by OpenSeesSP/MP.
echo   copying Intel MKL runtime DLLs
for %%D in (
    mkl_intel_thread.2.dll mkl_core.2.dll mkl_def.2.dll
    mkl_avx2.2.dll mkl_avx512.2.dll mkl_mc3.2.dll
    mkl_scalapack_lp64.2.dll mkl_blacs_intelmpi_lp64.2.dll
) do (
    if exist "%MKL_BIN%\%%D" copy /y "%MKL_BIN%\%%D" "%DIST%\bin\" >nul
)
if exist "%ICOMP_BIN%\libiomp5md.dll" copy /y "%ICOMP_BIN%\libiomp5md.dll" "%DIST%\bin\" >nul

REM Intel MPI runtime is NOT copied here. To launch OpenSeesSP / OpenSeesMP
REM the user must run them via mpiexec, which only works after
REM `setup_env.bat` has activated Intel oneAPI in the shell — that puts the
REM oneAPI version of impi.dll on PATH. A dist copy would load first and
REM subtly mismatch what mpiexec/MUMPS were built against (manifests as
REM MUMPS "Instance Error 1" under -n>1).

REM Tcl init.tcl (looked up relative to OpenSees.exe at runtime).
REM `conan cache path tcl/8.6.11` returns the export folder, not the package
REM folder where libs live. The package lives at
REM   %USERPROFILE%\.conan2\p\b\tcl<hash>\p\lib\tcl8.6\init.tcl
REM `dir /s /b` doesn't accept wildcards in directory parts, so we recurse
REM from the conan root and grep the result for the pattern we want.
for /f "delims=" %%I in ('dir /s /b /a-d "%USERPROFILE%\.conan2\init.tcl" 2^>nul ^| findstr "\\p\\lib\\tcl8.6\\init.tcl"') do (
    set "TCL_INIT=%%I"
    goto :got_tcl
)
:got_tcl
if defined TCL_INIT (
    for %%D in ("!TCL_INIT!") do set "TCL_LIB_DIR=%%~dpD"
    echo   copying Tcl library from !TCL_LIB_DIR!
    xcopy /e /q /y "!TCL_LIB_DIR!" "%DIST%\lib\tcl8.6\" >nul
) else (
    echo   WARNING: init.tcl not found in Conan cache; OpenSees.exe will warn at startup
)

echo.
echo === Done. Built artifacts ===
dir /b "%DIST%\bin"
echo.
echo To use:
echo   set PATH=%DIST%\bin;%%PATH%%        (for OpenSees/SP/MP exe variants)
echo   In Python:  sys.path.insert(0, r"%DIST%\bin"); import opensees
echo.
endlocal
