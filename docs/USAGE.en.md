# Using Topp

[中文](USAGE.md) · [README](../README.en.md) · [API](API.en.md) · [Mathematics](MATHEMATICS.en.md)

## Installation

```powershell
py -m pip install topp
```

Binary wheels support Windows x64 and CPython 3.10–3.14. Installing from source requires CMake 3.24+ and a C++20 compiler:

```powershell
py -m pip install -v .
```

Linux and macOS are not yet covered by CI; sdist builds on those platforms are unvalidated.

## Pairwise distances

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.3, 0.8]])
y = np.array([[0.0, 1.1], [0.4, 0.9]])

d = topp.bottleneck_distance(x, y)
w1 = topp.wasserstein_distance(x, y)
w2 = topp.wasserstein_distance(x, y, order=2, internal_p=2)
```

## PreparedDiagram and batches

Prepare a repeatedly used query once. `PreparedDiagram` is immutable and safe for concurrent reads.

```python
query = topp.prepare_diagram(x)
targets = [y, [[0.0, 2.0]], np.empty((0, 2))]

distances = topp.bottleneck_distances(query, targets)
wasserstein = topp.wasserstein_distances(query, targets, order=2, internal_p=2)
```

Reuse an existing output array:

```python
out = np.empty(len(targets), dtype=np.float64)
assert topp.bottleneck_distances(query, targets, out=out) is out
```

`out` must be writable, C-contiguous, `float64`, and shaped `(len(targets),)`.

## Threshold decisions

```python
if topp.bottleneck_within(x, y, 0.1):
    print("close")
```

The threshold must be non-negative and not NaN; `+inf` is valid.

## Input format

- shape `(n, 2)`, with each row `[birth, death]`;
- lists, tuples, NumPy arrays, and non-contiguous views are accepted;
- values are normalized to C-contiguous `float64`;
- use `[]` or `np.empty((0, 2))` for an empty diagram;
- `(finite, +inf)`, `(-inf, finite)`, and `(-inf, +inf)` are supported;
- invalid inputs raise exceptions and are never swapped, removed, or silently repaired.

Run [examples/basic.py](../examples/basic.py) to smoke-test an installation.
