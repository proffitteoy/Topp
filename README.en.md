<p align="center">
  <img src="https://raw.githubusercontent.com/proffitteoy/Topp/main/assets/topp-mark.svg" width="96" alt="Topp mark">
</p>

<h1 align="center">Topp</h1>

<p align="center">Exact Bottleneck and Wasserstein distances for persistence diagrams.</p>

<p align="center">
  <a href="https://github.com/proffitteoy/Topp/actions/workflows/ci.yml"><img src="https://github.com/proffitteoy/Topp/actions/workflows/ci.yml/badge.svg" alt="CI status"></a>
  <a href="https://pypi.org/project/topp/"><img src="https://img.shields.io/pypi/v/topp" alt="PyPI version"></a>
  <a href="https://pypi.org/project/topp/"><img src="https://img.shields.io/pypi/pyversions/topp?logo=python&logoColor=white" alt="Python versions"></a>
</p>

<p align="center">
  <a href="README.md">中文</a> ·
  <a href="docs/USAGE.en.md">Usage</a> ·
  <a href="docs/API.en.md">API</a> ·
  <a href="docs/DEVELOPMENT.en.md">Development</a> ·
  <a href="CHANGELOG.md">Changelog</a>
</p>

Topp provides a small Python API backed by an adaptive C++20 kernel. Prepare a diagram once and use native batch operations when comparing it with many candidates. NumPy is the only runtime dependency, and inline type information is included.

> **v0.1.0 public preview:** PyPI treats `0.1.0` as a final version number, while this project remains in public testing. The Python API is frozen; the Wasserstein kernel will continue to improve.

## Features

- exact Bottleneck distance with the internal `L∞` metric;
- exact `W1-L∞` and `W2-L2` Wasserstein distances;
- immutable `PreparedDiagram` objects and native one-to-many operations;
- reusable output arrays and exact `bottleneck_within` threshold decisions;
- GIL release during native computation;
- Windows x64 wheels for CPython 3.10–3.14;
- NumPy as the only runtime dependency; AVX2 is selected at runtime and is not required.

## Installation

```powershell
py -m pip install topp
```

The initial wheels target Windows x64. Building the sdist on other platforms requires CMake 3.24+ and a C++20 compiler.

## Quick start

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

See [examples/basic.py](examples/basic.py) for a complete runnable example.

### Batch comparisons with reusable memory

```python
targets = [y, np.empty((0, 2))]
out = np.empty(len(targets), dtype=np.float64)

query = topp.prepare_diagram(x)
topp.wasserstein_distances(
    query, targets, order=2, internal_p=2, out=out
)
```

## Supported metrics

| Function | Meaning | Status |
|---|---|---|
| `bottleneck_distance` | exact Bottleneck, internal `L∞` | supported |
| `wasserstein_distance(..., order=1, internal_p=np.inf)` | exact `W1-L∞` | supported |
| `wasserstein_distance(..., order=2, internal_p=2)` | exact `W2-L2` | supported |
| other Wasserstein parameters | potentially valid mathematically | `NotImplementedError` |

## Input contract

Inputs must convert to a `float64` array of shape `(n, 2)`. Empty diagrams, diagonal points, duplicates, and canonical essential points are valid. NaN, `birth > death`, `birth=+inf`, `death=-inf`, and other invalid infinity forms raise `ValueError` and are never silently repaired.

See the [API reference](docs/API.en.md) for details.

## Experimental features

The C++ source retains experimental candidate, graph, matching, component, and incremental-pricing strategies for maintainers. They are not exposed through the public Python API and are not default performance claims. See [docs/research](docs/research/README.md) for evidence and historical proposals.

## Development

```powershell
py -m pip install -v .
py -m pytest tests/python
cmd.exe /d /c scripts\build-kernel.cmd
```

The existing `include/bottleneck/*` C++ interface supports community maintenance and kernel experiments; it does not promise a stable ABI. See the [development guide](docs/DEVELOPMENT.en.md).

## Project links

| Entry | Contents |
|---|---|
| [Usage guide](docs/USAGE.en.md) | Installation, scalar and batch calls, output arrays, and errors |
| [API reference](docs/API.en.md) | Complete public API and input contract |
| [Development guide](docs/DEVELOPMENT.en.md) | Local builds, tests, and benchmarks |
| [Contributing](CONTRIBUTING.en.md) | Requirements for correctness and performance changes |
| [Research notes](docs/research/README.md) | Kernel experiments, differential evidence, and historical proposals |
| [Changelog](CHANGELOG.md) | Released capabilities and known limitations |

## Citing

For research use, cite the repository version and release tag. Machine-readable metadata is available in [CITATION.cff](CITATION.cff).

## License

Topp is licensed under the [MIT License](LICENSE). GUDHI is used only as a test oracle, semantic reference, and source of a historical patch; it is not a runtime dependency. See [third-party notices](THIRD_PARTY_NOTICES.md).
