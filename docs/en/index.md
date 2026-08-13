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

Median time in milliseconds for one query diagram against 64 target diagrams. Each size pools five synthetic input families; lower is faster.

| Distance | Points per diagram | Topp | GUDHI | giotto-tda | Topp vs GUDHI | Topp vs giotto-tda |
|---|---:|---:|---:|---:|---:|---:|
| Bottleneck | 8 | **1.128** | 10.942 | 16.032 | **9.70×** | **14.21×** |
| Bottleneck | 32 | **4.586** | 45.774 | 94.884 | **9.98×** | **20.69×** |
| Bottleneck | 128 | **62.490** | 237.496 | 590.716 | **3.80×** | **9.45×** |
| Bottleneck | 512 | **709.542** | 1,907.244 | 4,297.638 | **2.69×** | **6.06×** |
| Wasserstein | 8 | **1.316** | 17.513 | 13.699 | **13.31×** | **10.41×** |
| Wasserstein | 32 | **3.115** | 22.766 | 118.450 | **7.31×** | **38.03×** |
| Wasserstein | 128 | **10.308** | 108.834 | 1,023.854 | **10.56×** | **99.33×** |
| Wasserstein | 512 | **16.332** | 3,430.789 | 9,456.923 | **210.06×** | **579.03×** |

Measured on Windows 11 with Python 3.12 and one thread, using Topp 0.1.0, GUDHI 3.13.0, and giotto-tda 0.6.2. Bottleneck uses each library's default call (Topp is exact; GUDHI `e=None` and giotto-tda are approximate). For Wasserstein, Topp and GUDHI compute exact $W_{1,\infty}$, while giotto-tda defaults to approximate $W_2$; that column therefore describes default Python API speed, not an algorithm ranking for the same mathematical task.
