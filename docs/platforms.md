# 平台与支持

[English](https://proffitteoy.github.io/Topp/en/platforms.html)

## 已验证的发布范围

| 组件 | 0.1.0 版本状态 |
|---|---|
| Python | CPython 3.10–3.14 |
| 预编译 wheels | Windows x64 |
| 运行时依赖 | NumPy 1.23 或更高版本 |
| 源码构建编译器 | C++20 |
| 构建系统 | CMake 3.24 或更高版本 |

CI 在 Windows Server 2022 上构建并测试 C++ 内核、标量回退、Python wheels 和 sdist 安装。

## Linux 与 macOS

源码设计为可移植的 C++20，但 Linux 和 macOS 尚未纳入 0.1.0 版本的 CI。在这些系统上可能成功完成源码构建，但目前不属于已验证的发布平台。

## CPU 调度

Windows x64 构建可能包含 AVX2 内核。是否可用会在运行时检测；AVX2 不是必需条件，不支持的 CPU 会使用标量实现。

## 公开与实验接口

本站记录的 Python API 是公开兼容性范围。C++ 头文件和内部策略开关用于维护及内核实验；0.1.0 版本不承诺稳定的 C++ ABI。
