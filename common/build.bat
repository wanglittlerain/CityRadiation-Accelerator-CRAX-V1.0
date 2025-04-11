@echo off
mkdir msvc_build
cd msvc_build
cmake .. -G "Visual Studio 17" -Wno-dev
cmake --build . --config Release