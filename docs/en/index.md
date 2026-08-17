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

| Distance | Points per diagram | Topp | GUDHI | Hera | Topp vs GUDHI | Topp vs Hera |
|---|---:|---:|---:|---:|---:|---:|
| Bottleneck | 8 | **1.109** | 11.435 | 17.558 | **10.31×** | **15.83×** |
| Bottleneck | 32 | **4.458** | 49.380 | 132.501 | **11.08×** | **29.72×** |
| Bottleneck | 128 | **11.553** | 265.807 | 765.597 | **23.01×** | **66.27×** |
| Bottleneck | 512 | **53.879** | 1,913.019 | 5,077.575 | **35.51×** | **94.24×** |
| Wasserstein | 8 | **1.161** | 16.813 | 6.976 | **14.48×** | **6.01×** |
| Wasserstein | 32 | **3.154** | 22.152 | 66.302 | **7.02×** | **21.02×** |
| Wasserstein | 128 | **11.987** | 112.666 | 537.364 | **9.40×** | **44.83×** |
| Wasserstein | 512 | **16.400** | 3,112.252 | 4,159.430 | **189.77×** | **253.62×** |

Measured on 2026-08-17 on Windows 11 with Python 3.12.13 and one thread, using an MSVC wheel built from the 1.0 kernel baseline commit `4cf5b4e` (SHA-256 `562161cdc20dda8ab751f611102dda71e5d9d5585c1feed30e0110db1934939c`), GUDHI 3.13.0, and the Hera API bundled with that release. Bottleneck uses each library's default call: Topp is exact, GUDHI `e=None` uses its approximate strategy, and Hera uses approximate `delta=0.01`. For Wasserstein, Topp and GUDHI compute exact $W_{1,\infty}$, while Hera computes approximate $W_{1,\infty}$ (`delta=0.01`); the Hera column therefore describes default Python API speed, not an algorithm ranking at equal precision. Runtime varies substantially by input family, so the pooled result does not imply the same speedup for every distribution.
