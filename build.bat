@echo off
cd common
call build.bat
cd %~dp0
mkdir msvc_build
cd msvc_build
cmake -DCMAKE_GENERATOR_PLATFORM=x64 ..
cmake --build . --config Release
pause
exit