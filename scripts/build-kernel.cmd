@echo off
setlocal

if not defined BOTTLENECK_BUILD_DIR set "BOTTLENECK_BUILD_DIR=build\manual"
set "BOTTLENECK_COMPILE_OPT=/O2"
set "BOTTLENECK_LINK_OPT="
if /I "%BOTTLENECK_ENABLE_LTO%"=="1" (
  set "BOTTLENECK_COMPILE_OPT=/O2 /GL"
  set "BOTTLENECK_LINK_OPT=/link /LTCG"
)

call "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64
if errorlevel 1 exit /b %errorlevel%

if not exist "%BOTTLENECK_BUILD_DIR%" mkdir "%BOTTLENECK_BUILD_DIR%"

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /arch:AVX2 /c ^
  src\distance_avx2.cpp /Fo:%BOTTLENECK_BUILD_DIR%\distance_avx2.obj
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /arch:AVX2 /c ^
  src\wasserstein_avx2.cpp /Fo:%BOTTLENECK_BUILD_DIR%\wasserstein_avx2.obj
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp tests\bottleneck_core_tests.cpp %BOTTLENECK_BUILD_DIR%\distance_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\bottleneck_core_tests.exe %BOTTLENECK_LINK_OPT%
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp benchmarks\bottleneck_core_bench.cpp %BOTTLENECK_BUILD_DIR%\distance_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\bottleneck_core_bench.exe %BOTTLENECK_LINK_OPT%
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp benchmarks\bottleneck_grid_bench.cpp %BOTTLENECK_BUILD_DIR%\distance_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\bottleneck_grid_bench.exe %BOTTLENECK_LINK_OPT%
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp benchmarks\bottleneck_batch_bench.cpp %BOTTLENECK_BUILD_DIR%\distance_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\bottleneck_batch_bench.exe %BOTTLENECK_LINK_OPT%
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp benchmarks\bottleneck_large_bench.cpp %BOTTLENECK_BUILD_DIR%\distance_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\bottleneck_large_bench.exe %BOTTLENECK_LINK_OPT%
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /DBOTTLENECK_HAVE_WASSERSTEIN_AVX2=1 /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp src\wasserstein.cpp tests\wasserstein_core_tests.cpp %BOTTLENECK_BUILD_DIR%\distance_avx2.obj %BOTTLENECK_BUILD_DIR%\wasserstein_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\wasserstein_core_tests.exe %BOTTLENECK_LINK_OPT%
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /DBOTTLENECK_HAVE_WASSERSTEIN_AVX2=1 /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp src\wasserstein.cpp benchmarks\wasserstein_core_bench.cpp %BOTTLENECK_BUILD_DIR%\distance_avx2.obj %BOTTLENECK_BUILD_DIR%\wasserstein_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\wasserstein_core_bench.exe %BOTTLENECK_LINK_OPT%
if errorlevel 1 exit /b %errorlevel%

cl.exe /nologo /std:c++20 %BOTTLENECK_COMPILE_OPT% /EHsc /W4 /permissive- /Zc:__cplusplus /DBOTTLENECK_HAVE_AVX2_KERNEL=1 /LD /Iinclude ^
  src\bottleneck_core.cpp src\geometric_backend.cpp tests\bottleneck_core_c_api.cpp %BOTTLENECK_BUILD_DIR%\distance_avx2.obj ^
  /Fe:%BOTTLENECK_BUILD_DIR%\bottleneck_core_c.dll %BOTTLENECK_LINK_OPT%
if errorlevel 1 exit /b %errorlevel%

exit /b 0
