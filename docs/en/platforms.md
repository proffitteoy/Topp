# Platforms and support

[中文](https://proffitteoy.github.io/Topp/platforms.html)

## Validated release surface

| Component | Version 1.0 status |
|---|---|
| Python | CPython 3.10–3.14 |
| Prebuilt wheels | Windows x64 and Linux x86_64 (manylinux) |
| Runtime dependency | NumPy 1.23 or newer |
| Compiler for source builds | C++20 |
| Build system | CMake 3.24 or newer |

CI builds and tests the C++ kernel, Python wheels, and sdist installation on Windows Server 2022 and Ubuntu across CPython 3.10–3.14. Local Linux validation uses WSL2 Ubuntu 24.04, GCC 13.3, and Python 3.12.

## Linux and macOS

Linux x86_64 is a validated release platform starting with 1.0.1. The release pipeline builds manylinux wheels and installs each Python version from TestPyPI. macOS remains a portable-C++20 source target but is not part of release CI.

## CPU dispatch

Windows x64 builds may include AVX2 kernels. Availability is checked at runtime; AVX2 is not required. Linux wheels currently use the portable scalar kernel and do not promise the same performance as Windows AVX2 builds.

## Public and experimental interfaces

The Python API documented on this site is the public `1.x` compatibility surface. The C++ headers and internal strategy switches support maintenance and kernel experiments and do not promise a stable C++ ABI.
