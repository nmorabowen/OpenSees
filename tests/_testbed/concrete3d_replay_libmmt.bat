@echo off
REM concrete3d_replay_libmmt.bat <replay|replay-ucrt|selfcheck|selfcheck-ucrt> [outdir]   (run from anywhere; apeGmsh #2621 F1)
REM
REM Builds the LadrunoConcrete3D kernel testbeds the way the shipped Windows opensees.pyd is built: MSVC with the fork's
REM Release flags (/O2 /Ob2 /DNDEBUG /EHsc /std:c++17 /MT, no /arch, no /fp) and the math library the pyd links. The
REM pyd's atan2/cos/sin/exp/log/pow/acos resolve to Intel's libmmt, pulled in by the Fortran objects'
REM /DEFAULTLIB:libmmt directives, NOT to the UCRT a plain MSVC link uses; eig3sym differs by 1-13 ulp between the two.
REM   replay          concrete3d_replay.cpp vs concrete3d_replay_records.txt (recorded by the pyd under libmmt), built
REM                   against the kernel the records were written by (the fixture's kernel-blob); exit 0 iff bit-exact
REM   replay-ucrt     the same replay with a plain (UCRT) link and --cross: shows what a wrong libm does to the records
REM   selfcheck       concrete3d_kernel_check.cpp under libmmt, against this tree's kernel. KNOWN ISSUE: the pristine
REM                   d14fb6958 kernel fails C7 (1 of 120000 CDPM2-potential fuzz trials) under libmmt and passes under
REM                   the UCRT; it is NOT skipped here. Unset LADRUNO_FUZZ_TRIALS for the full trial count.
REM   selfcheck-ucrt  the same self-check with a plain (UCRT) link, for the comparison
REM Needs Visual Studio (vswhere) and the oneAPI compiler's static libraries (ONEAPI_ROOT or the default install).
setlocal EnableExtensions EnableDelayedExpansion
set "MODE=%~1"
set "OUT=%~2"
if "%OUT%"=="" set "OUT=%TEMP%\c3d_replay_libmmt"
set "REPO=%~dp0..\.."
pushd "%REPO%" || exit /b 2
set "REPO=%CD%"
if not exist "%OUT%" mkdir "%OUT%"

set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" (echo ERROR: vswhere not found & popd & exit /b 2)
for /f "usebackq delims=" %%i in (`"%VSWHERE%" -latest -products * -property installationPath`) do set "VSDIR=%%i"
if not exist "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" (echo ERROR: vcvars64.bat not found under "%VSDIR%" & popd & exit /b 2)
call "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
if "%ONEAPI_ROOT%"=="" set "ONEAPI_ROOT=%ProgramFiles(x86)%\Intel\oneAPI"
set "INTEL_LIB=%ONEAPI_ROOT%\compiler\latest\lib"
set "LIBMMT=/link /LIBPATH:"%INTEL_LIB%" libmmt.lib libircmt.lib svml_dispmt.lib libirc.lib"
set "CLFLAGS=/nologo /O2 /Ob2 /DNDEBUG /EHsc /std:c++17 /MT"

if /i "%MODE%"=="replay" goto replay
if /i "%MODE%"=="replay-ucrt" goto replay_ucrt
if /i "%MODE%"=="selfcheck" goto selfcheck
if /i "%MODE%"=="selfcheck-ucrt" goto selfcheck_ucrt
echo usage: %~nx0 ^<replay^|replay-ucrt^|selfcheck^|selfcheck-ucrt^> [outdir]
popd & exit /b 2

:needmmt
if not exist "%INTEL_LIB%\libmmt.lib" (echo ERROR: libmmt.lib not found in "%INTEL_LIB%": this build would not link the pyd's libm & popd & exit /b 2)
exit /b 0

:replay_ucrt
set "TAG=win-ucrt"
set "LINKLIBS="
set "CROSS=--cross"
goto replay_build

:replay
call :needmmt || exit /b 2
set "TAG=win-libmmt"
set "LINKLIBS=%LIBMMT%"
set "CROSS="

:replay_build
set "FIX=%REPO%\tests\_testbed\concrete3d_replay_records.txt"
set "WANT="
for /f "usebackq tokens=*" %%l in (`findstr /b /c:"#@ " "%FIX%"`) do set "META=%%l"
set "META=!META:*kernel-blob=!"
for /f "tokens=1" %%a in ("!META:~1!") do set "WANT=%%a"
if "%WANT%"=="" (echo ERROR: the fixture names no kernel-blob & popd & exit /b 2)
for /f "usebackq delims=" %%h in (`git hash-object SRC\material\nD\LadrunoConcrete3DKernel.h`) do set "HAVE=%%h"
set "INC=%REPO%"
if not "%HAVE%"=="%WANT%" (
  echo kernel header %HAVE% differs from the recording kernel %WANT%: building the replay against the recorded blob
  if not exist "%OUT%\ref\SRC\material\nD" mkdir "%OUT%\ref\SRC\material\nD"
  git cat-file blob %WANT% > "%OUT%\ref\SRC\material\nD\LadrunoConcrete3DKernel.h" || (echo ERROR: blob %WANT% not in this clone & popd & exit /b 2)
  set "INC=%OUT%\ref"
)
cl %CLFLAGS% /DC3D_REPLAY_LIBM=\"%TAG%\" /I"%INC%" /Fo:"%OUT%\replay_%TAG%.obj" /Fe:"%OUT%\concrete3d_replay_%TAG%.exe" tests\_testbed\concrete3d_replay.cpp %LINKLIBS% >"%OUT%\build_replay_%TAG%.log" 2>&1
if errorlevel 1 (type "%OUT%\build_replay_%TAG%.log" & popd & exit /b 2)
"%OUT%\concrete3d_replay_%TAG%.exe" check "%FIX%" %CROSS%
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%

:selfcheck
call :needmmt || exit /b 2
cl %CLFLAGS% /I"%REPO%" /Fo:"%OUT%\c3dchk_libmmt.obj" /Fe:"%OUT%\c3dchk_libmmt.exe" tests\_testbed\concrete3d_kernel_check.cpp %LIBMMT% >"%OUT%\build_selfcheck.log" 2>&1
if errorlevel 1 (type "%OUT%\build_selfcheck.log" & popd & exit /b 2)
echo kernel self-check under Intel libmmt (the pyd's libm); LADRUNO_FUZZ_TRIALS=%LADRUNO_FUZZ_TRIALS%
"%OUT%\c3dchk_libmmt.exe" tests\_testbed\concrete3d_oracle_fixture.txt
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%

:selfcheck_ucrt
cl %CLFLAGS% /I"%REPO%" /Fo:"%OUT%\c3dchk_ucrt.obj" /Fe:"%OUT%\c3dchk_ucrt.exe" tests\_testbed\concrete3d_kernel_check.cpp >"%OUT%\build_selfcheck_ucrt.log" 2>&1
if errorlevel 1 (type "%OUT%\build_selfcheck_ucrt.log" & popd & exit /b 2)
echo kernel self-check under the UCRT (NOT the pyd's libm; comparison only); LADRUNO_FUZZ_TRIALS=%LADRUNO_FUZZ_TRIALS%
"%OUT%\c3dchk_ucrt.exe" tests\_testbed\concrete3d_oracle_fixture.txt
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%
