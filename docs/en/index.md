# Topp

[中文](https://proffitteoy.github.io/Topp/)

```{toctree}
:hidden:
:caption: Getting started

getting-started/installation
getting-started/quickstart
```

```{toctree}
:hidden:
:caption: Distance calculations

guide/bottleneck
guide/wasserstein
guide/prepared-batch
guide/thresholds
guide/input-semantics
guide/threads
```

```{toctree}
:hidden:
:caption: Reference

API
MATHEMATICS
platforms
```

```{toctree}
:hidden:
:caption: Project

DEVELOPMENT
USAGE
Changelog <https://github.com/proffitteoy/Topp/blob/main/CHANGELOG.md>
Source code <https://github.com/proffitteoy/Topp>
```

[![CI](https://github.com/proffitteoy/Topp/actions/workflows/ci.yml/badge.svg)](https://github.com/proffitteoy/Topp/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/topp)](https://pypi.org/project/topp/)
[![Python](https://img.shields.io/pypi/pyversions/topp)](https://pypi.org/project/topp/)

Topp is a focused Python package for **exact distances between persistence diagrams**. It is for users who already have diagrams—from GUDHI, Ripser, or another pipeline—and need strict, repeated comparisons in Python.

Topp currently implements:

- exact Bottleneck distance with internal $L^\infty$;
- exact $W_{1,\infty}$ and $W_{2,2}$ Wasserstein distance;
- immutable prepared diagrams and native one-to-many calls;
- reusable output arrays and exact Bottleneck threshold decisions.

The computational kernel is written in C++20. NumPy is the only runtime dependency.

```{important}
Topp compares persistence diagrams; it does not construct filtrations or compute persistence. It is not a replacement for a complete TDA library.
```

## Install

```console
py -m pip install topp
```

Prebuilt wheels currently support Windows x64 and CPython 3.10–3.14. See [Platforms and support](platforms.md) before attempting a source build elsewhere.

## First calculation

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.3, 0.8]])
y = np.array([[0.0, 1.1], [0.4, 0.9]])

print(topp.bottleneck_distance(x, y))
print(topp.wasserstein_distance(x, y, order=2, internal_p=2))
```

Continue with the [five-minute quickstart](getting-started/quickstart.md), or select a task from the left navigation.

## Python batch-distance performance

Median time in milliseconds for one query diagram against 64 target diagrams. For each size, all timing rounds from five synthetic input families (uniform, near-diagonal, clustered, duplicate-heavy, and separated) are pooled before taking the median; lower is faster.

| Distance | Points per diagram | Topp | GUDHI | giotto-tda | Topp vs GUDHI | Topp vs giotto-tda |
|---|---:|---:|---:|---:|---:|---:|
| Bottleneck | 8 | **1.059** | 10.729 | 14.434 | **10.13×** | **13.62×** |
| Bottleneck | 32 | **4.547** | 43.860 | 90.804 | **9.65×** | **19.97×** |
| Bottleneck | 128 | **16.632** | 235.851 | 558.250 | **14.18×** | **33.56×** |
| Bottleneck | 512 | **236.045** | 1,882.070 | 3,639.137 | **7.97×** | **15.42×** |
| Wasserstein | 8 | **1.104** | 15.344 | 11.758 | **13.90×** | **10.65×** |
| Wasserstein | 32 | **3.246** | 21.727 | 114.656 | **6.69×** | **35.33×** |
| Wasserstein | 128 | **10.101** | 106.639 | 1,000.598 | **10.56×** | **99.06×** |
| Wasserstein | 512 | **12.689** | 2,986.669 | 8,808.310 | **235.37×** | **694.17×** |

Measured on 2026-08-17 on Windows 11 with Python 3.12.13 and one thread, using the current-mainline Topp 0.1.0 MSVC wheel (SHA-256 `dffa357a504121538d63e3fea3675054430f67e04f91c35828f29a18d150f641`), GUDHI 3.13.0, and giotto-tda 0.6.2. Bottleneck uses each library's default call (Topp is exact; GUDHI `e=None` and giotto-tda are approximate). For Wasserstein, Topp and GUDHI compute exact $W_{1,\infty}$, while giotto-tda defaults to approximate $W_2$; that column therefore describes default Python API speed, not an algorithm ranking for the same mathematical task. Runtime varies substantially by input family, so the pooled result does not imply the same speedup for every distribution.
