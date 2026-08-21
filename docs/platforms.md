# 平台与支持

[English](https://proffitteoy.github.io/Topp/en/platforms.html)

## 已验证的发布范围

| 组件 | 1.0 版本状态 |
|---|---|
| Python | CPython 3.10–3.14 |
| 预编译 wheels | Windows x64、Linux x86_64（manylinux） |
| 运行时依赖 | NumPy 1.23 或更高版本 |
| 源码构建编译器 | C++20 |
| 构建系统 | CMake 3.24 或更高版本 |

CI 在 Windows Server 2022 和 Ubuntu 上构建并测试 C++ 内核、Python wheels 与 sdist 安装，并覆盖 CPython 3.10–3.14。Linux 本地验证使用 WSL2 Ubuntu 24.04、GCC 13.3 和 Python 3.12。

## Linux 与 macOS

Linux x86_64 是 1.0.1 起的已验证发布平台，发布链生成 manylinux wheels，并在 TestPyPI 阶段逐 Python 版本安装验证。macOS 源码仍以可移植 C++20 为目标，但尚未纳入发布 CI。

## CPU 调度

Windows x64 构建可能包含 AVX2 内核。是否可用会在运行时检测；AVX2 不是必需条件，不支持的 CPU 会使用标量实现。Linux wheel 当前使用可移植标量内核，不承诺与 Windows AVX2 构建相同的性能。

## 公开与实验接口

本站记录的 Python API 属于 `1.x` 公开兼容性范围。C++ 头文件和内部策略开关用于维护及内核实验，不承诺稳定的 C++ ABI。
