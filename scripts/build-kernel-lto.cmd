@echo off
setlocal
set "BOTTLENECK_ENABLE_LTO=1"
set "BOTTLENECK_BUILD_DIR=build\lto"
call "%~dp0build-kernel.cmd"
exit /b %errorlevel%
