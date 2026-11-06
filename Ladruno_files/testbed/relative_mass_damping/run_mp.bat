@echo off
call "C:\Program Files (x86)\Intel\oneAPI\mpi\latest\env\vars.bat"
set "BIN=%~1"
set "PATH=%BIN%;%PATH%"
py -3.12 "%~dp0rr_mp.py" "%BIN%\OpenSeesMP.exe" mpiexec "%~2" %3
