cmake_minimum_required(VERSION 3.24)

project(bottleneck_core VERSION 0.1.0 LANGUAGES CXX)

option(BOTTLENECK_BUILD_TESTS "Build kernel correctness tests" ON)
option(BOTTLENECK_BUILD_BENCHMARKS "Build kernel microbenchmarks" ON)
option(BOTTLENECK_ENABLE_AVX2 "Build runtime-dispatched AVX2 kernels on MSVC x64" ON)
option(TOPP_BUILD_PYTHON "Build the topp Python extension" OFF)

add_library(bottleneck_core STATIC
    src/bottleneck_core.cpp
    src/geometric_backend.cpp
    src/wasserstein.cpp
)
add_library(bottleneck::core ALIAS bottleneck_core)
set_target_properties(bottleneck_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(bottleneck_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(bottleneck_core PUBLIC cxx_std_20)

if(BOTTLENECK_ENABLE_AVX2 AND MSVC AND CMAKE_SIZEOF_VOID_P EQUAL 8)
    target_sources(bottleneck_core PRIVATE src/distance_avx2.cpp src/wasserstein_avx2.cpp)
    set_source_files_properties(src/distance_avx2.cpp src/wasserstein_avx2.cpp
        PROPERTIES COMPILE_OPTIONS /arch:AVX2)
    target_compile_definitions(bottleneck_core PRIVATE
        BOTTLENECK_HAVE_AVX2_KERNEL=1
        BOTTLENECK_HAVE_WASSERSTEIN_AVX2=1)
endif()

if(TOPP_BUILD_PYTHON)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)

    pybind11_add_module(_core MODULE python/bindings.cpp)
    target_link_libraries(_core PRIVATE bottleneck::core)
    target_compile_features(_core PRIVATE cxx_std_20)
    install(TARGETS _core LIBRARY DESTINATION topp RUNTIME DESTINATION topp)
endif()

if(MSVC)
    target_compile_options(bottleneck_core PRIVATE /W4 /permissive- /EHsc)
else()
    target_compile_options(bottleneck_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(BOTTLENECK_BUILD_TESTS)
    enable_testing()
    add_executable(bottleneck_core_tests tests/bottleneck_core_tests.cpp)
    target_link_libraries(bottleneck_core_tests PRIVATE bottleneck::core)
    target_compile_features(bottleneck_core_tests PRIVATE cxx_std_20)
    add_test(NAME bottleneck_core_tests COMMAND bottleneck_core_tests)

    add_executable(wasserstein_core_tests tests/wasserstein_core_tests.cpp)
    target_link_libraries(wasserstein_core_tests PRIVATE bottleneck::core)
    target_compile_features(wasserstein_core_tests PRIVATE cxx_std_20)
    add_test(NAME wasserstein_core_tests COMMAND wasserstein_core_tests)

    add_library(bottleneck_core_c SHARED tests/bottleneck_core_c_api.cpp)
    target_link_libraries(bottleneck_core_c PRIVATE bottleneck::core)
    target_compile_features(bottleneck_core_c PRIVATE cxx_std_20)
endif()

if(BOTTLENECK_BUILD_BENCHMARKS)
    add_executable(bottleneck_core_bench benchmarks/bottleneck_core_bench.cpp)
    target_link_libraries(bottleneck_core_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_core_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_grid_bench benchmarks/bottleneck_grid_bench.cpp)
    target_link_libraries(bottleneck_grid_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_grid_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_batch_bench benchmarks/bottleneck_batch_bench.cpp)
    target_link_libraries(bottleneck_batch_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_batch_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_large_bench benchmarks/bottleneck_large_bench.cpp)
    target_link_libraries(bottleneck_large_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_large_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_multiplicity_bench benchmarks/bottleneck_multiplicity_bench.cpp)
    target_link_libraries(bottleneck_multiplicity_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_multiplicity_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_mandatory_bench benchmarks/bottleneck_mandatory_bench.cpp)
    target_link_libraries(bottleneck_mandatory_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_mandatory_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_component_bench benchmarks/bottleneck_component_bench.cpp)
    target_link_libraries(bottleneck_component_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_component_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_geometry_bench benchmarks/bottleneck_geometry_bench.cpp)
    target_link_libraries(bottleneck_geometry_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_geometry_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_search_bench benchmarks/bottleneck_search_bench.cpp)
    target_link_libraries(bottleneck_search_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_search_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_router_bench benchmarks/bottleneck_router_bench.cpp)
    target_link_libraries(bottleneck_router_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_router_bench PRIVATE cxx_std_20)

    add_executable(bottleneck_wasserstein_assist_bench
      benchmarks/bottleneck_wasserstein_assist_bench.cpp)
    target_link_libraries(bottleneck_wasserstein_assist_bench PRIVATE bottleneck::core)
    target_compile_features(bottleneck_wasserstein_assist_bench PRIVATE cxx_std_20)

    add_executable(wasserstein_core_bench benchmarks/wasserstein_core_bench.cpp)
    target_link_libraries(wasserstein_core_bench PRIVATE bottleneck::core)
    target_compile_features(wasserstein_core_bench PRIVATE cxx_std_20)

    add_executable(wasserstein_batch_bench benchmarks/wasserstein_batch_bench.cpp)
    target_link_libraries(wasserstein_batch_bench PRIVATE bottleneck::core)
    target_compile_features(wasserstein_batch_bench PRIVATE cxx_std_20)
endif()
