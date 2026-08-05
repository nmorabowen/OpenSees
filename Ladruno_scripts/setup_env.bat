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

REM Ladruno: locate vcvars64.bat machine-agnostically via vswhere instead of a
REM hardcoded VS 2022 path. Works across VS editions/years (2022, 2026/v18, ...)
REM by asking vswhere for the newest install carrying the C++ x64 toolset. Falls
REM back to the historical hardcoded 2022 Community path if vswhere is absent.
set "VS_VCVARS="
if exist "%VSWHERE_DIR%\vswhere.exe" (
    for /f "usebackq delims=" %%I in (`"%VSWHERE_DIR%\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do (
        if exist "%%I\VC\Auxiliary\Build\vcvars64.bat" set "VS_VCVARS=%%I\VC\Auxiliary\Build\vcvars64.bat"
    )
)
if not defined VS_VCVARS set "VS_VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if not exist "%VS_VCVARS%" goto :err_vs
echo === Using MSVC: %VS_VCVARS% ===
call "%VS_VCVARS%"
if errorlevel 1 goto :err_vcvars

REM ---- 2-4. Intel oneAPI components --------------------------------------
set "ONEAPI_ROOT=C:\Program Files (x86)\Intel\oneAPI"

REM ---- 2a. Pick a FORTRAN-CAPABLE compiler (Ladruno ADR-75 / oneAPI 2026) --
REM  oneAPI 2026.1 ships NO Fortran: `compiler\2026.1` has neither ifx.exe nor
REM  ifconsol.lib (Intel split Fortran into its own package). The moment an
REM  update re-points the `compiler\latest` junction to 2026.x, this fork's
REM  builds break in two different places, and the second one only bites MPI:
REM     * no ifx      -> "CMAKE_Fortran_COMPILER: ifx ... was not found in the
REM                       PATH" at the MUMPS configure (Step 1 of build.bat)
REM     * no ifconsol -> "LINK : fatal error LNK1104: cannot open file
REM                       'ifconsol.lib'" when linking OpenSeesMP/SP/PyMP,
REM                       because MKL's ScaLAPACK/BLACS import libraries still
REM                       carry a /DEFAULTLIB:ifconsol directive. Serial targets
REM                       never link ScaLAPACK, so they keep building fine --
REM                       which is exactly why this hides until an MP build.
REM  This was observed LIVE: an update flipped `latest` 2025.1 -> 2026.1 between
REM  two builds in one session (the earlier log shows
REM  "Fortran compiler identification is IntelLLVM 2025.1.1" from latest\bin).
REM  Fix: prefer `latest`, but if it has no ifx.exe fall back to the
REM  newest installed compiler that does, and source THAT vars.bat.
REM  NB: written FLAT (labels, no parenthesized read-after-write) because this
REM  script must not `setlocal` -- it exports the environment to its caller --
REM  so delayed expansion (!VAR!) is unavailable here.
set "OPS_ICC_DIR=%ONEAPI_ROOT%\compiler\latest"
if exist "%OPS_ICC_DIR%\bin\ifx.exe" goto :icc_chosen
set "OPS_ICC_DIR="
for /f "delims=" %%D in ('dir /b /ad /o-n "%ONEAPI_ROOT%\compiler" 2^>nul') do if not defined OPS_ICC_DIR if exist "%ONEAPI_ROOT%\compiler\%%D\bin\ifx.exe" set "OPS_ICC_DIR=%ONEAPI_ROOT%\compiler\%%D"
if defined OPS_ICC_DIR goto :icc_fallback
echo ERROR: no installed oneAPI compiler provides ifx.exe ^(Fortran^).
echo        oneAPI 2026+ ships Fortran separately - install the Fortran
echo        package, or keep a 2025.x compiler alongside it.
goto :err_compiler
:icc_fallback
echo   [fortran] compiler\latest has no ifx.exe - falling back to %OPS_ICC_DIR%
:icc_chosen
call "%OPS_ICC_DIR%\env\vars.bat" intel64
if errorlevel 1 goto :err_compiler

REM ---- 2b. ifconsol.lib fallback -----------------------------------------
REM  Same root cause; kept separate because a compiler CAN have ifx but not
REM  ifconsol. Purely additive to LIB.
if exist "%OPS_ICC_DIR%\lib\ifconsol.lib" goto :ifconsol_ok
set "_IFCONSOL_DIR="
for /f "delims=" %%D in ('dir /b /ad /o-n "%ONEAPI_ROOT%\compiler" 2^>nul') do if not defined _IFCONSOL_DIR if exist "%ONEAPI_ROOT%\compiler\%%D\lib\ifconsol.lib" set "_IFCONSOL_DIR=%ONEAPI_ROOT%\compiler\%%D\lib"
if defined _IFCONSOL_DIR goto :ifconsol_add
echo   [ifconsol] WARNING: no installed compiler provides ifconsol.lib;
echo              MPI targets ^(OpenSeesMP/SP/PyMP^) will fail to link ^(LNK1104^).
goto :ifconsol_ok
:ifconsol_add
echo   [ifconsol] active compiler lacks ifconsol.lib - adding %_IFCONSOL_DIR%
set "LIB=%_IFCONSOL_DIR%;%LIB%"
:ifconsol_ok

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
