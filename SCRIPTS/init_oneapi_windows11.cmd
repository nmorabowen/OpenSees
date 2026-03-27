@echo off
set "VS2022INSTALLDIR="

if exist "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools" set "VS2022INSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
if not defined VS2022INSTALLDIR if exist "C:\Program Files\Microsoft Visual Studio\2022\Community" set "VS2022INSTALLDIR=C:\Program Files\Microsoft Visual Studio\2022\Community"
if not defined VS2022INSTALLDIR if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional" set "VS2022INSTALLDIR=C:\Program Files\Microsoft Visual Studio\2022\Professional"
if not defined VS2022INSTALLDIR if exist "C:\Program Files\Microsoft Visual Studio\2022\Enterprise" set "VS2022INSTALLDIR=C:\Program Files\Microsoft Visual Studio\2022\Enterprise"

if defined VS2022INSTALLDIR echo Using VS2022INSTALLDIR=%VS2022INSTALLDIR%
if not defined VS2022INSTALLDIR echo WARNING: Visual Studio 2022 was not found. oneAPI will continue without MSVC setup.

call "%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat" intel64
