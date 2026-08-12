# Topp

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

API.en
MATHEMATICS.en
platforms
```

```{toctree}
:hidden:
:caption: Project

DEVELOPMENT.en
USAGE.en
Chinese documentation <zh/index>
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

## Performance boundary

The repository contains kernel experiments and differential tests. These validate implementations and guide the adaptive default, but they are not a uniform Python-level comparison with other libraries. Until a controlled cross-library suite is published, Topp does not claim an end-to-end speed advantage over GUDHI, Hera, or Persim.
