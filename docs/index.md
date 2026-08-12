---
layout: default
title: Topp
---

# Topp

[![CI](https://github.com/proffitteoy/Topp/actions/workflows/ci.yml/badge.svg)](https://github.com/proffitteoy/Topp/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/topp)](https://pypi.org/project/topp/)
[![Python](https://img.shields.io/pypi/pyversions/topp)](https://pypi.org/project/topp/)

Topp is a lean Python package for exact distances between persistence diagrams. It does one part of a TDA workflow: given diagrams that have already been computed, Topp compares them using exact Bottleneck, `W1-L∞`, or `W2-L2` distance.

The package is intended for applications that repeatedly compare one query diagram with many candidates. Its public API includes immutable prepared diagrams, native one-to-many calls, reusable output arrays, and an exact Bottleneck threshold decision. The computational kernel is implemented in C++20; NumPy is the only runtime dependency.

Topp does **not** compute persistence diagrams or filtrations. It is not a replacement for a complete library such as GUDHI, Ripser, or Dionysus.

## Installation

Prebuilt wheels are currently available for Windows x64 and CPython 3.10–3.14:

```console
py -m pip install topp
```

Source builds require CMake 3.24 or newer and a C++20 compiler. Linux and macOS source builds are not covered by CI in version 0.1.0.

## Example

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.3, 0.8]])
y = np.array([[0.0, 1.1], [0.4, 0.9]])

print(topp.bottleneck_distance(x, y))
print(topp.wasserstein_distance(x, y, order=2, internal_p=2))

query = topp.prepare_diagram(x)
print(topp.bottleneck_distances(query, [y, np.empty((0, 2))]))
print(topp.bottleneck_within(query, y, 0.1))
```

## Supported distances

| Python call | Mathematical convention |
|---|---|
| `bottleneck_distance(X, Y)` | exact Bottleneck with internal `L∞` |
| `wasserstein_distance(X, Y)` | exact `W1-L∞` |
| `wasserstein_distance(X, Y, order=2, internal_p=2)` | exact `W2-L2` |

Other Wasserstein parameter combinations raise `NotImplementedError`. Duplicate points retain multiplicity. Essential points are supported; incompatible essential multiplicities produce an infinite distance. The complete definitions are given in [Mathematical conventions](MATHEMATICS.en.html).

## Prepared and batch calls

Use `prepare_diagram` when the same diagram participates in multiple comparisons:

```python
query = topp.prepare_diagram(x)
targets = [y, [[0.0, 2.0]], np.empty((0, 2))]

out = np.empty(len(targets), dtype=np.float64)
topp.bottleneck_distances(query, targets, out=out)
```

`PreparedDiagram` is immutable and may be read concurrently. Native distance calls release the Python GIL. See [Using Topp](USAGE.en.html) and the [Python API](API.en.html) for the full contract.

## Performance claims

The repository contains internal kernel experiments and differential tests. They establish correctness and guide the adaptive default, but they are not a controlled Python-level comparison with other libraries. Until a uniform cross-library benchmark is published, Topp does not claim an end-to-end speed advantage over GUDHI, Hera, or Persim.

## Documentation

- [Using Topp](USAGE.en.html)
- [Python API](API.en.html)
- [Mathematical conventions](MATHEMATICS.en.html)
- [Development](DEVELOPMENT.en.html)
- [Source repository](https://github.com/proffitteoy/Topp)
- [Changelog](https://github.com/proffitteoy/Topp/blob/main/CHANGELOG.md)

Chinese versions are available for [usage](USAGE.html), the [API](API.html), [mathematical conventions](MATHEMATICS.html), and [development](DEVELOPMENT.html).

## Citing

If Topp is used in research, cite the software version and release tag. Machine-readable citation metadata is available in [`CITATION.cff`](https://github.com/proffitteoy/Topp/blob/main/CITATION.cff).

## License

Topp is distributed under the [MIT License](https://github.com/proffitteoy/Topp/blob/main/LICENSE). GUDHI is used as a test oracle and semantic reference, not as a runtime dependency.
