# Developing Topp

[中文](https://proffitteoy.github.io/Topp/DEVELOPMENT.html) · [Contributing](https://github.com/proffitteoy/Topp/blob/main/CONTRIBUTING.en.md)

Regular users consume only the Python API. `include/bottleneck/*`, experimental strategies, and benchmarks support community maintenance and do not promise a stable C++ ABI.

## Python build and tests

Requirements: Python 3.10+, CMake 3.24+, and a C++20 compiler. Visual Studio 2022 Build Tools are recommended on Windows.

```powershell
py -m pip install -v .
py -m pip install pytest
py -m pytest tests/python
```

The oracle suite requires pinned GUDHI 3.13.0:

```powershell
$env:TOPP_RUN_ORACLE_TESTS='1'
py -m pytest tests/python/test_oracles.py
```

## C++ kernel

```powershell
cmd.exe /d /c scripts\build-kernel.cmd
build\manual\bottleneck_core_tests.exe
build\manual\wasserstein_core_tests.exe
```

Disable AVX2 translation units to validate the scalar build:

```powershell
cmake -S . -B build\scalar -G Ninja -DCMAKE_BUILD_TYPE=Release -DBOTTLENECK_ENABLE_AVX2=OFF
cmake --build build\scalar
ctest --test-dir build\scalar --output-on-failure
```

## Benchmarks

See [benchmarks/README.md](https://github.com/proffitteoy/Topp/blob/main/benchmarks/README.md) for commands and fair-comparison rules. Raw output belongs in ignored `benchmarks/results/`; commit only reproducible commands, fixed inputs/seeds, and bounded conclusions.

Experimental strategies must pass reference/oracle differentials before randomized-order median/p95 comparisons. Keep conclusions for failed routes and do not add them to the default dispatcher.

## Layout

- `python/topp/`: public Python layer and type information;
- `python/bindings.cpp`: private pybind11 extension;
- `include/`, `src/`: exact C++20 kernel;
- `tests/python/`, `tests/*.cpp`: Python/C++ correctness;
- `benchmarks/`: reproducible kernel benchmarks;
- `docs/research/`: evidence, historical proposals, and experimental patches.
