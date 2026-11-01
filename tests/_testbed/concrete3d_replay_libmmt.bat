@echo off
REM concrete3d_replay_libmmt.bat <replay|replay-ucrt|contract|record|selfcheck|selfcheck-ucrt> [outdir]   (apeGmsh #2621 F1)
REM concrete3d_replay_libmmt.bat replay-kernel <outdir> <root>   (the kernel-mutant proof, see below)
REM
REM Builds the LadrunoConcrete3D kernel testbeds the way the shipped Windows opensees.pyd is built: MSVC with the fork's
REM Release flags (/O2 /Ob2 /DNDEBUG /EHsc /std:c++17 /MT, no /arch, no /fp) and the math library the pyd links. The
REM pyd's atan2/cos/sin/exp/log/pow/acos resolve to Intel's libmmt, pulled in by the Fortran objects'
REM /DEFAULTLIB:libmmt directives, NOT to the UCRT a plain MSVC link uses; eig3sym differs by 1-13 ulp between the two.
REM Every build uses THIS TREE's kernel header; there is no fallback to the header the records were written by.
REM   replay          concrete3d_replay.cpp vs concrete3d_replay_records.txt (written under libmmt); exit 0 iff every
REM                   record replays bit-exactly. A kernel change that moves any value FAILS here by design: that PR runs
REM                   `record`, commits the re-recorded fixture, and its diff (exactly the values that moved) is what the
REM                   reviewers judge. No field or record class is exempt.
REM   replay-ucrt     the same replay with a plain (UCRT) link and --cross: shows what a wrong libm does to the records
REM   contract        the field-classification contract under libmmt (ci.sh --contract's mirror, on the libmmt
REM                   fixture): the header's FIELD TABLE block must be the compiled table (`fields --check`), and every
REM                   row `contract` generates from it must give its class's outcome; exit 0 iff all do. Each CONTRACT
REM                   FAIL line names the field and its class.
REM   replay-kernel   replay under libmmt built against <root>\SRC\material\nD\LadrunoConcrete3DKernel.h instead of
REM                   this tree's (a mutated copy: ci.sh's kernel mutants), with that header's own blob. Must FAIL (exit 1).
REM   record          re-record the fixture under libmmt on this tree's kernel into <outdir>\concrete3d_replay_records.txt
REM   selfcheck       concrete3d_kernel_check.cpp under libmmt, full trial count. The d14fb6958 kernel failed C7 (1 of
REM                   120000 CDPM2-potential fuzz trials) under libmmt and passed under the UCRT; ladruno's kernel (#877)
REM                   passes under both (ALL PASS, measured on 6b7ebf3a8). The script still refuses to run when
REM                   LADRUNO_FUZZ_TRIALS is set: a smaller count (e.g. the Zone-A 200) never reaches a rare trial.
REM   selfcheck-ucrt  the same self-check with a plain (UCRT) link, for the comparison (same refusal)
REM Needs Visual Studio (vswhere), git, and the oneAPI compiler's static libraries (ONEAPI_ROOT or the default install).
setlocal EnableExtensions
REM only replay-kernel sets INCROOT on purpose: one inherited from the caller's environment must not redirect a build
set "INCROOT="
set "MODE=%~1"
set "OUT=%~2"
if "%OUT%"=="" set "OUT=%TEMP%\c3d_replay_libmmt"
set "REPO=%~dp0..\.."
pushd "%REPO%" || exit /b 2
set "REPO=%CD%"
if not exist "%OUT%" mkdir "%OUT%"
set "FIX=%REPO%\tests\_testbed\concrete3d_replay_records.txt"

set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" (echo ERROR: vswhere not found & popd & exit /b 2)
for /f "usebackq delims=" %%i in (`"%VSWHERE%" -latest -products * -property installationPath`) do set "VSDIR=%%i"
if not exist "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" (echo ERROR: vcvars64.bat not found under "%VSDIR%" & popd & exit /b 2)
call "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
if "%ONEAPI_ROOT%"=="" set "ONEAPI_ROOT=%ProgramFiles(x86)%\Intel\oneAPI"
set "INTEL_LIB=%ONEAPI_ROOT%\compiler\latest\lib"
set "LIBMMT=/link /LIBPATH:"%INTEL_LIB%" libmmt.lib libircmt.lib svml_dispmt.lib libirc.lib"
set "CLFLAGS=/nologo /O2 /Ob2 /DNDEBUG /EHsc /std:c++17 /MT"
set "HAVE="
for /f "usebackq delims=" %%h in (`git hash-object SRC\material\nD\LadrunoConcrete3DKernel.h`) do set "HAVE=%%h"
if "%HAVE%"=="" (echo ERROR: git hash-object of the kernel header failed & popd & exit /b 2)

if /i "%MODE%"=="replay" goto replay
if /i "%MODE%"=="replay-ucrt" goto replay_ucrt
if /i "%MODE%"=="contract" goto contract
if /i "%MODE%"=="replay-kernel" goto replay_kernel
if /i "%MODE%"=="record" goto record
if /i "%MODE%"=="selfcheck" goto selfcheck
if /i "%MODE%"=="selfcheck-ucrt" goto selfcheck_ucrt
echo usage: %~nx0 ^<replay^|replay-ucrt^|contract^|record^|selfcheck^|selfcheck-ucrt^> [outdir]   ^|   replay-kernel ^<outdir^> ^<root^>
popd & exit /b 2

:needmmt
if not exist "%INTEL_LIB%\libmmt.lib" (echo ERROR: libmmt.lib not found in "%INTEL_LIB%": this build would not link the pyd's libm & exit /b 2)
exit /b 0

:fulltrials
if defined LADRUNO_FUZZ_TRIALS (echo ERROR: LADRUNO_FUZZ_TRIALS=%LADRUNO_FUZZ_TRIALS% is set: the self-check must run the full trial count: a reduced count never reaches a rare trial. Unset it. & exit /b 2)
exit /b 0

:build_replay
REM %1 = libm tag; LINKLIBS = the link libraries (empty for the UCRT); INCROOT = the kernel header's root (default the repo)
if "%INCROOT%"=="" set "INCROOT=%REPO%"
cl %CLFLAGS% /DC3D_REPLAY_LIBM=\"%~1\" /I"%INCROOT%" /Fo:"%OUT%\replay_%~1.obj" /Fe:"%OUT%\concrete3d_replay_%~1.exe" tests\_testbed\concrete3d_replay.cpp %LINKLIBS% >"%OUT%\build_replay_%~1.log" 2>&1
if errorlevel 1 (type "%OUT%\build_replay_%~1.log" & exit /b 2)
echo built against the kernel header under %INCROOT%, blob %HAVE%
exit /b 0

:replay
call :needmmt || (popd & exit /b 2)
set "LINKLIBS=%LIBMMT%"
call :build_replay win-libmmt || (popd & exit /b 2)
"%OUT%\concrete3d_replay_win-libmmt.exe" check "%FIX%" --kernel-blob %HAVE%
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%

:replay_ucrt
set "LINKLIBS="
call :build_replay win-ucrt || (popd & exit /b 2)
"%OUT%\concrete3d_replay_win-ucrt.exe" check "%FIX%" --cross --kernel-blob %HAVE%
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%

:contract
call :needmmt || (popd & exit /b 2)
set "LINKLIBS=%LIBMMT%"
call :build_replay win-libmmt || (popd & exit /b 2)
"%OUT%\concrete3d_replay_win-libmmt.exe" fields --check tests\_testbed\concrete3d_replay.cpp
if errorlevel 1 (echo ERROR: the header's FIELD TABLE block is not the compiled field table & popd & exit /b 1)
"%OUT%\concrete3d_replay_win-libmmt.exe" contract "%FIX%" --kernel-blob %HAVE%
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" echo ERROR: the field-classification contract failed (exit %RC%): each CONTRACT FAIL line names the field and its class
popd & exit /b %RC%

:replay_kernel
call :needmmt || (popd & exit /b 2)
set "INCROOT=%~3"
if not exist "%INCROOT%\SRC\material\nD\LadrunoConcrete3DKernel.h" (echo ERROR: no kernel header under "%INCROOT%" & popd & exit /b 2)
set "HAVE="
for /f "usebackq delims=" %%h in (`git hash-object "%INCROOT%\SRC\material\nD\LadrunoConcrete3DKernel.h"`) do set "HAVE=%%h"
if "%HAVE%"=="" (echo ERROR: git hash-object of the mutated kernel header failed & popd & exit /b 2)
set "LINKLIBS=%LIBMMT%"
call :build_replay win-libmmt || (popd & exit /b 2)
"%OUT%\concrete3d_replay_win-libmmt.exe" check "%FIX%" --kernel-blob %HAVE%
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%

:record
call :needmmt || (popd & exit /b 2)
set "LINKLIBS=%LIBMMT%"
call :build_replay win-libmmt || (popd & exit /b 2)
"%OUT%\concrete3d_replay_win-libmmt.exe" record "%FIX%" "%OUT%\concrete3d_replay_records.txt" --kernel-blob %HAVE%
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%

:selfcheck
call :fulltrials || (popd & exit /b 2)
call :needmmt || (popd & exit /b 2)
cl %CLFLAGS% /I"%REPO%" /Fo:"%OUT%\c3dchk_libmmt.obj" /Fe:"%OUT%\c3dchk_libmmt.exe" tests\_testbed\concrete3d_kernel_check.cpp %LIBMMT% >"%OUT%\build_selfcheck.log" 2>&1
if errorlevel 1 (type "%OUT%\build_selfcheck.log" & popd & exit /b 2)
echo kernel self-check under Intel libmmt (the pyd's libm), full trial count, kernel blob %HAVE%
"%OUT%\c3dchk_libmmt.exe" tests\_testbed\concrete3d_oracle_fixture.txt
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%

:selfcheck_ucrt
call :fulltrials || (popd & exit /b 2)
cl %CLFLAGS% /I"%REPO%" /Fo:"%OUT%\c3dchk_ucrt.obj" /Fe:"%OUT%\c3dchk_ucrt.exe" tests\_testbed\concrete3d_kernel_check.cpp >"%OUT%\build_selfcheck_ucrt.log" 2>&1
if errorlevel 1 (type "%OUT%\build_selfcheck_ucrt.log" & popd & exit /b 2)
echo kernel self-check under the UCRT (NOT the pyd's libm; comparison only), full trial count, kernel blob %HAVE%
"%OUT%\c3dchk_ucrt.exe" tests\_testbed\concrete3d_oracle_fixture.txt
set "RC=%ERRORLEVEL%"
popd & exit /b %RC%
