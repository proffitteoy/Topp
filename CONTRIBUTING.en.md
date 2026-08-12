# Contributing to Topp

[中文](CONTRIBUTING.md) · [Development guide](docs/DEVELOPMENT.en.md)

Correctness, performance, portability, testing, examples, and documentation contributions are welcome.

1. Keep changes focused and identify whether they affect the public Python API, the stable default kernel, or an experimental strategy.
2. Public behavior changes must update bilingual API/Usage docs, type information, and Python tests.
3. Kernel optimizations must pass the C++ reference, Python contract, and applicable external-oracle differentials before randomized-order median/p95 reporting.
4. Do not commit `build/`, `.obj`, wheels, raw benchmark output, virtual environments, or third-party source trees.
5. Performance claims must record the commit, compiler, CPU, distribution, size, repetitions, and correctness boundary.

Run at least:

```powershell
py -m pytest tests/python
cmd.exe /d /c scripts\build-kernel.cmd
build\manual\bottleneck_core_tests.exe
build\manual\wasserstein_core_tests.exe
```
