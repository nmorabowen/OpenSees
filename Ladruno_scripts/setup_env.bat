@echo off
REM ==========================================================================
REM  setup_env.bat — load the OpenSees build toolchain into this shell
REM
REM  We don't use Intel's setvars.bat because it iterates over every installed
REM  oneAPI component and the iteration breaks on this machine. We instead
REM  call the three component env\vars.bat scripts we actually need:
REM    - compiler  (ifx Fortran, plus its runtime DLLs)
REM    - mkl       (BLAS / LAPACK / ScaLAPACK static libs)
REM    - mpi       (Intel MPI: mpiexec, mpi.lib, mpi_f.lib, includes)
REM
REM  Order:
REM    1. vcvars64.bat → MSVC cl.exe, link.exe, Windows SDK
REM    2. compiler\env\vars.bat
REM    3. mkl\env\vars.bat
REM    4. mpi\env\vars.bat
REM    5. PATH for Conan
REM
REM  Usage (from a fresh cmd.exe):  call setup_env.bat
REM ==========================================================================

REM ---- 1. MSVC toolchain (x64) -------------------------------------------
REM vswhere.exe is shipped with the VS installer; add it to PATH so that any
REM downstream tool (e.g. Conan's per-package conanvcvars.bat, TCL's nmake
REM build) can locate Visual Studio.
set "VSWHERE_DIR=C:\Program Files (x86)\Microsoft Visual Studio\Installer"
if exist "%VSWHERE_DIR%\vswhere.exe" set "PATH=%VSWHERE_DIR%;%PATH%"

set "VS_VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if not exist "%VS_VCVARS%" goto :err_vs
call "%VS_VCVARS%"
if errorlevel 1 goto :err_vcvars

REM ---- 2-4. Intel oneAPI components --------------------------------------
set "ONEAPI_ROOT=C:\Program Files (x86)\Intel\oneAPI"

call "%ONEAPI_ROOT%\compiler\latest\env\vars.bat" intel64
if errorlevel 1 goto :err_compiler

call "%ONEAPI_ROOT%\mkl\latest\env\vars.bat" intel64
if errorlevel 1 goto :err_mkl

call "%ONEAPI_ROOT%\mpi\latest\env\vars.bat"
if errorlevel 1 goto :err_mpi

REM Intel MPI sometimes needs I_MPI_ROOT explicitly (some installers don't set it)
if "%I_MPI_ROOT%"=="" set "I_MPI_ROOT=%ONEAPI_ROOT%\mpi\latest"

REM ---- 5. Conan (installed via `py -3.12 -m pip install --user conan`) --
set "CONAN_DIR=%APPDATA%\Python\Python312\Scripts"
if exist "%CONAN_DIR%\conan.exe" set "PATH=%CONAN_DIR%;%PATH%"

REM ---- 6. Workaround: TCL 8.6.11's makefile.vc assumes that cmd will
REM        search the current directory for executables (it builds
REM        nmakehlp.exe in cwd and invokes it as bare `nmakehlp`). Modern
REM        Windows sets NoDefaultCurrentDirectoryInExePath=1 by default for
REM        security, which breaks the TCL build. We clear it for THIS shell
REM        session only — the user/machine setting is untouched.
set "NoDefaultCurrentDirectoryInExePath="

REM ---- 6. Sanity check ---------------------------------------------------
echo.
echo === Toolchain loaded ===
for %%T in (cl.exe ifx.exe mpiexec.exe cmake.exe ninja.exe conan.exe) do call :checktool %%T
echo.
echo MKL_ROOT=%MKLROOT%
echo I_MPI_ROOT=%I_MPI_ROOT%
echo.

REM ---- 7. Sentinel: mark the toolchain as loaded in THIS shell ------------
REM build.bat checks this so a second build in the same window can skip the
REM slow vcvars64 + 3x oneAPI vars.bat reload. Persists only when this script
REM is `call`ed directly from an interactive shell (build.bat's own setlocal
REM discards it at endlocal, which is why pre-loading once is what saves time).
set "LADRUNO_TOOLCHAIN_LOADED=1"
goto :eof

:checktool
where %1 >nul 2>&1
if errorlevel 1 (echo   [MISSING] %1) else (echo   [OK]      %1)
exit /b 0

:err_vs
echo ERROR: Cannot find Visual Studio 2022 vcvars64.bat at:
echo   %VS_VCVARS%
exit /b 1
:err_vcvars
echo ERROR: vcvars64.bat failed.
exit /b 1
:err_compiler
echo ERROR: compiler\env\vars.bat failed.
exit /b 1
:err_mkl
echo ERROR: mkl\env\vars.bat failed.
exit /b 1
:err_mpi
echo ERROR: mpi\env\vars.bat failed.
exit /b 1
