# Topp

[中文](README.md) · [Usage](docs/USAGE.en.md) · [API](docs/API.en.md) · [Development](docs/DEVELOPMENT.en.md)

Topp is a lean persistence-diagram distance package. It focuses on exact Bottleneck distance and two exact Wasserstein metrics, backed by an adaptive C++20 kernel.

> **v0.1.0 public preview:** PyPI treats `0.1.0` as a final version number, while this project remains in public testing. The Python API is frozen; the Wasserstein kernel will continue to improve.

## Features

- exact Bottleneck distance with the internal `L∞` metric;
- exact `W1-L∞` and `W2-L2` Wasserstein distances;
- immutable `PreparedDiagram` objects;
- native one-to-many operations and caller-provided output arrays;
- exact threshold decisions with `bottleneck_within`;
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

## Citing

For research use, cite the repository version and release tag. Machine-readable metadata is available in [CITATION.cff](CITATION.cff).

## License

Topp is licensed under the [MIT License](LICENSE). GUDHI is used only as a test oracle, semantic reference, and source of a historical patch; it is not a runtime dependency. See [third-party notices](THIRD_PARTY_NOTICES.md).
