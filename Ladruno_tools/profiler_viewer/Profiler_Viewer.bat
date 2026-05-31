@echo off
REM One-click launcher for the Ladruno profiler viewer (Windows).
REM Double-click to open the viewer, or drag a profile.h5 file onto this .bat.
REM It finds Python, then hands off to launch.py (which provisions a venv,
REM installs deps, builds the UI if needed, and opens the browser).

setlocal
cd /d "%~dp0"

REM Locate a Python launcher / interpreter.
set "PYEXE="
where py >nul 2>nul && set "PYEXE=py -3"
if not defined PYEXE (
    where python >nul 2>nul && set "PYEXE=python"
)
if not defined PYEXE (
    echo.
    echo   Python 3 was not found on PATH.
    echo   Install Python 3.10+ from https://www.python.org/downloads/
    echo   (tick "Add python.exe to PATH" in the installer^), then run this again.
    echo.
    pause
    exit /b 1
)

REM %~1 is an optional profile.h5 (e.g. drag-and-drop onto the .bat).
%PYEXE% "%~dp0launch.py" %1
if errorlevel 1 (
    echo.
    echo   The viewer exited with an error. See the messages above.
    pause
)
endlocal
