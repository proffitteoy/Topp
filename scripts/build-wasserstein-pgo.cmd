@echo off
setlocal

if not defined BOTTLENECK_BUILD_DIR set "BOTTLENECK_BUILD_DIR=build\pgo-msvc"
set "BOTTLENECK_PGO_MODE=%~1"
set "BOTTLENECK_PGO_TARGET=%~2"
if not defined BOTTLENECK_PGO_TARGET set "BOTTLENECK_PGO_TARGET=bench"
set "BOTTLENECK_PGO_DATABASE=%~3"

if /I "%BOTTLENECK_PGO_MODE%"=="instrument" (
  set "BOTTLENECK_PGO_LINK=/GENPROFILE"
) else if /I "%BOTTLENECK_PGO_MODE%"=="optimize" (
  set "BOTTLENECK_PGO_LINK=/USEPROFILE"
) else (
  echo Usage: build-wasserstein-pgo.cmd instrument^|optimize [bench^|tests] [profile-database-name]
  exit /b 2
)

if /I "%BOTTLENECK_PGO_TARGET%"=="bench" (
  set "BOTTLENECK_PGO_SOURCE=benchmarks\wasserstein_core_bench.cpp"
  set "BOTTLENECK_PGO_NAME=wasserstein_core_bench"
) else if /I "%BOTTLENECK_PGO_TARGET%"=="tests" (
  set "BOTTLENECK_PGO_SOURCE=tests\wasserstein_core_tests.cpp"
  set "BOTTLENECK_PGO_NAME=wasserstein_core_tests"
) else (
  echo Unknown target: %BOTTLENECK_PGO_TARGET%
  exit /b 2
)
if not defined BOTTLENECK_PGO_DATABASE set "BOTTLENECK_PGO_DATABASE=%BOTTLENECK_PGO_NAME%"

call "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64
if errorlevel 1 exit /b %errorlevel%

if not exist "%BOTTLENECK_BUILD_DIR%" mkdir "%BOTTLENECK_BUILD_DIR%"

cl.exe /nologo /std:c++20 /O2 /GL /EHsc /W4 /permissive- /Zc:__cplusplus /arch:AVX2 /c ^
  src\distance_avx2.cpp /Fo:%BOTTLENECK_BUILD_DIR%\distance_avx2.obj
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 /O2 /GL /EHsc /W4 /permissive- /Zc:__cplusplus /arch:AVX2 /c ^
  src\wasserstein_avx2.cpp /Fo:%BOTTLENECK_BUILD_DIR%\wasserstein_avx2.obj
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 /O2 /GL /EHsc /W4 /permissive- /Zc:__cplusplus ^
  /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /DBOTTLENECK_HAVE_WASSERSTEIN_AVX2=1 /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp src\wasserstein.cpp ^
  %BOTTLENECK_PGO_SOURCE% %BOTTLENECK_BUILD_DIR%\distance_avx2.obj ^
  %BOTTLENECK_BUILD_DIR%\wasserstein_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\%BOTTLENECK_PGO_NAME%.exe ^
  /link /LTCG %BOTTLENECK_PGO_LINK% ^
  /PGD:%BOTTLENECK_BUILD_DIR%\%BOTTLENECK_PGO_DATABASE%.pgd
if errorlevel 1 exit /b %errorlevel%

if /I "%BOTTLENECK_PGO_MODE%"=="instrument" (
  copy /Y "%VCToolsInstallDir%bin\Hostx64\x64\pgort140.dll" ^
    "%BOTTLENECK_BUILD_DIR%\pgort140.dll" >nul
  if errorlevel 1 exit /b %errorlevel%
)

exit /b 0
