# Topp 开发指南

[English](DEVELOPMENT.en.md) · [贡献指南](../CONTRIBUTING.md)

普通用户只使用 Python API。`include/bottleneck/*`、实验策略和 benchmark 面向社区维护者，不承诺稳定 C++ ABI。

## Python 构建与测试

需要 Python 3.10+、CMake 3.24+、C++20 编译器。Windows 推荐 Visual Studio 2022 Build Tools。

```powershell
py -m pip install -v .
py -m pip install pytest
py -m pytest tests/python
```

Oracle 套件需要固定的 GUDHI 3.13.0：

```powershell
$env:TOPP_RUN_ORACLE_TESTS='1'
py -m pytest tests/python/test_oracles.py
```

## C++ 内核

```powershell
cmd.exe /d /c scripts\build-kernel.cmd
build\manual\bottleneck_core_tests.exe
build\manual\wasserstein_core_tests.exe
```

关闭 AVX2 translation units 以验证 scalar 构建：

```powershell
cmake -S . -B build\scalar -G Ninja -DCMAKE_BUILD_TYPE=Release -DBOTTLENECK_ENABLE_AVX2=OFF
cmake --build build\scalar
ctest --test-dir build\scalar --output-on-failure
```

## Benchmark

命令和公平比较规则见 [benchmarks/README.md](../benchmarks/README.md)。原始输出放在已忽略的 `benchmarks/results/`；只提交可复现命令、固定输入/种子和带边界的汇总结论。

实验策略必须先通过 reference/oracle 差分，再比较随机顺序 median/p95。失败路线保留结论，不进入默认 dispatcher。

## 目录

- `python/topp/`：公开 Python 层和类型提示；
- `python/bindings.cpp`：私有 pybind11 扩展；
- `include/`、`src/`：C++20 exact 内核；
- `tests/python/`、`tests/*.cpp`：Python/C++ 正确性；
- `benchmarks/`：可复现内核基准；
- `docs/research/`：证据、历史方案和实验补丁。
