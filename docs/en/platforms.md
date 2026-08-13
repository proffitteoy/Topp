# Platforms and support

[中文](https://proffitteoy.github.io/Topp/platforms.html)

## Validated release surface

| Component | Version 0.1.0 status |
|---|---|
| Python | CPython 3.10–3.14 |
| Prebuilt wheels | Windows x64 |
| Runtime dependency | NumPy 1.23 or newer |
| Compiler for source builds | C++20 |
| Build system | CMake 3.24 or newer |

CI builds and tests the C++ kernel, scalar fallback, Python wheels, and sdist installation on Windows Server 2022.

## Linux and macOS

The source is intended to be portable C++20, but Linux and macOS are not part of the version 0.1.0 CI matrix. Source builds on those systems may work, but they are not currently validated release platforms.

## CPU dispatch

Windows x64 builds may include AVX2 kernels. Availability is checked at runtime; AVX2 is not required. Unsupported CPUs use the scalar implementation.

## Public and experimental interfaces

The Python API documented on this site is the public compatibility surface. The C++ headers and internal strategy switches support maintenance and kernel experiments; version 0.1.0 does not promise a stable C++ ABI.
