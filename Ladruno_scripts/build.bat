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
REM                                              (OpenSees, OpenSeesSP, OpenSeesMP,
REM                                               OpenSeesPy, OpenSeesPyMP)
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
set "IMPI_LIBFABRIC=C:\Program Files (x86)\Intel\oneAPI\mpi\latest\opt\mpi\libfabric\bin"

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

set "TARGETS=OpenSees OpenSeesSP OpenSeesMP OpenSeesPy OpenSeesPyMP"
if not "%MODE%"=="" set "TARGETS=%MODE%"

REM ----- Load toolchain ----------------------------------------------------
REM setup_env.bat (vcvars64 + 3x oneAPI vars.bat) is the slow part of a build.
REM Skip it when the toolchain is already active in this shell:
REM   - `set SKIP_SETUP_ENV=1` before running, to force-skip (you manage the
REM     toolchain yourself), OR
REM   - `call Ladruno_scripts\setup_env.bat` ONCE in your cmd window, then run
REM     build.bat repeatedly -- the LADRUNO_TOOLCHAIN_LOADED sentinel is
REM     inherited and we auto-skip the reload.
REM Within a single build.bat run the sentinel can't carry across runs (this
REM script's setlocal discards it at endlocal), so the pre-load-once workflow
REM is what actually saves time.
if defined SKIP_SETUP_ENV (
    echo === Toolchain reload skipped ^(SKIP_SETUP_ENV set^) ===
) else if defined LADRUNO_TOOLCHAIN_LOADED (
    echo === Toolchain already loaded in this shell; skipping setup_env.bat ===
) else (
    echo === Loading toolchain ===
    call "%SCRIPT_DIR%\setup_env.bat" || (echo setup_env.bat failed & exit /b 1)
)

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

REM ----- 2b. Pin CMake to the conan-provided copy --------------------------
REM A system-wide CMake 4.x on PATH shadows the conan-pinned 3.31.x and makes
REM the `conan-release` preset configure fail (ZLIB-not-found / preset
REM mismatch) WITHOUT regenerating build.ninja, so the build silently dies at
REM Step 3 with stale binaries left in dist\. Locate the conan cache cmake.exe
REM (now guaranteed present: conan install just ran) and put its dir FIRST on
REM PATH so every `cmake` call below resolves to the matching version.
REM Discovery mirrors the init.tcl cache-grep idiom further down; no hash is
REM hardcoded. See project_build_env_cmake43_conan_zlib build note.
echo.
echo === Step 2b: Pinning CMake to the conan-provided copy ===
set "CONAN_CMAKE="
for /f "delims=" %%I in ('dir /s /b /a-d "%USERPROFILE%\.conan2\cmake.exe" 2^>nul ^| findstr /i "\\p\\bin\\cmake.exe"') do set "CONAN_CMAKE=%%I"
if defined CONAN_CMAKE (
    for %%D in ("!CONAN_CMAKE!") do set "CONAN_CMAKE_DIR=%%~dpD"
    set "PATH=!CONAN_CMAKE_DIR!;!PATH!"
    echo   pinned: !CONAN_CMAKE!
) else (
    echo   WARNING: conan CMake not found in cache; falling back to PATH cmake.
    echo            If configure fails with ZLIB-not-found, a system CMake 4.x
    echo            is likely shadowing the conan 3.31.x — see the build notes.
)
for /f "delims=" %%V in ('cmake --version 2^>^&1 ^| findstr /i "version"') do echo   using %%V

REM ----- 2c. Heal a cache written by a foreign cmake ------------------------
REM If a shadowing cmake ever got as far as configure, it overwrote
REM CMakeCache.txt; reconfiguring that cache with the conan cmake is asking
REM for stale/poisoned results. The cache records which binary wrote it
REM (CMAKE_COMMAND). If that isn't the conan-cache copy, delete JUST the
REM cache file: the next configure is clean, while compiled objects survive
REM (ninja re-runs only what its command hashes say changed).
if defined CONAN_CMAKE if exist "%BUILD_DIR%\CMakeCache.txt" (
    findstr /b /i /c:"CMAKE_COMMAND:INTERNAL=" "%BUILD_DIR%\CMakeCache.txt" | findstr /i ".conan2" >nul
    if errorlevel 1 (
        echo   CMakeCache.txt was written by a non-conan cmake -- deleting it
        echo   ^(one-time clean reconfigure; compiled objects are preserved^)
        del /q "%BUILD_DIR%\CMakeCache.txt"
    )
)

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
    REM `|| (...)` reliably catches a non-zero exit of the immediately
    REM preceding command. The old `if errorlevel 1 (...)` inside this
    REM parenthesized for-body failed to trip on a cmake/ninja link failure
    REM (openseesmp.dll LNK1120) and let Step 5 run on a non-built target.
    cmake --build "%BUILD_DIR%" --target %%T -j 8 || (echo Build of %%T failed & exit /b 1)
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

REM ----- OpenSeesPyMP: self-contained MPI Python package ------------------
REM Ladruno Patch 9. Built as openseesmp.dll (CMake OUTPUT_NAME); Python
REM expects openseesmp.pyd. Shipped to a DEDICATED dir, NOT %DIST%\bin: a
REM dist\bin-local impi.dll shadows the oneAPI one and breaks the co-located
REM Tcl OpenSeesSP/MP.exe (MUMPS "Instance Error 1" under -n>1; see the note
REM below). openseesmp.pyd links impi.dll at load time, so it gets a
REM co-located, same-oneAPI-version Intel MPI runtime + Hydra launcher here.
REM It imports as `openseesmp` (rename shim) so it never collides with the
REM sequential `import opensees`. Launch with this dir's own mpiexec, e.g.
REM   %DIST%\openseesmp\mpiexec -n 4 python driver.py
if not exist "%BUILD_DIR%\openseesmp.dll" goto :no_openseesmp
if not exist "%DIST%\openseesmp" mkdir "%DIST%\openseesmp"
echo   copying openseesmp.dll -^> openseesmp\openseesmp.pyd
copy /y "%BUILD_DIR%\openseesmp.dll" "%DIST%\openseesmp\openseesmp.pyd" >nul
echo   bundling Intel MPI runtime + Hydra launcher into openseesmp\
for %%M in (impi.dll mpiexec.exe hydra_bstrap_proxy.exe hydra_pmi_proxy.exe hydra_service.exe) do (
    if exist "%IMPI_BIN%\%%M" copy /y "%IMPI_BIN%\%%M" "%DIST%\openseesmp\" >nul
)
if exist "%IMPI_LIBFABRIC%\libfabric.dll" copy /y "%IMPI_LIBFABRIC%\libfabric.dll" "%DIST%\openseesmp\" >nul
echo   mirroring MKL runtime into openseesmp\ (self-contained off oneAPI shell)
for %%D in (
    mkl_intel_thread.2.dll mkl_core.2.dll mkl_def.2.dll
    mkl_avx2.2.dll mkl_avx512.2.dll mkl_mc3.2.dll
    mkl_scalapack_lp64.2.dll mkl_blacs_intelmpi_lp64.2.dll
) do (
    if exist "%MKL_BIN%\%%D" copy /y "%MKL_BIN%\%%D" "%DIST%\openseesmp\" >nul
)
if exist "%ICOMP_BIN%\libiomp5md.dll" copy /y "%ICOMP_BIN%\libiomp5md.dll" "%DIST%\openseesmp\" >nul
:no_openseesmp

REM (OpenSeesPySP packaging removed — Python has no SP subsystem; SP is
REM  Tcl-only via OpenSeesSP.exe. See Ladruno Patch 9 docs.)

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

REM Intel MPI runtime is deliberately NOT copied into %DIST%\bin. To launch
REM the Tcl OpenSeesSP / OpenSeesMP exes the user runs them via mpiexec,
REM which only works after `setup_env.bat` has activated Intel oneAPI in the
REM shell — that puts the oneAPI impi.dll on PATH. A dist\bin copy would
REM load first and subtly mismatch what mpiexec/MUMPS were built against
REM (manifests as MUMPS "Instance Error 1" under -n>1). The Python MP
REM module does NOT reintroduce this hazard: its Intel MPI runtime lives in
REM the separate, self-contained %DIST%\openseesmp\ (same oneAPI version as
REM its bundled mpiexec, never on the Tcl exes' PATH) — see the OpenSeesPyMP
REM block above.

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
if exist "%DIST%\openseesmp\openseesmp.pyd" (
    echo.
    echo === openseesmp package ===
    dir /b "%DIST%\openseesmp"
)
echo.
echo To use:
echo   set PATH=%DIST%\bin;%%PATH%%        (for OpenSees/SP/MP exe variants)
echo   In Python:  sys.path.insert(0, r"%DIST%\bin"); import opensees
echo   MPI Python:  set PATH=%DIST%\openseesmp;%%PATH%%
echo                %DIST%\openseesmp\mpiexec -n 4 python driver.py   (import openseesmp)
echo.
endlocal
