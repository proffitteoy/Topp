# Prepared diagrams and batch calls

[中文](https://proffitteoy.github.io/Topp/guide/prepared-batch.html)

Use the prepared and one-to-many APIs when one diagram is compared repeatedly.

## Prepare once

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.25, 0.75]])
query = topp.prepare_diagram(x)

print(query.n_points)         # all valid input rows
print(query.n_finite_points)  # finite points above the diagonal
```

The prepared object is immutable and owns its input copy:

```python
x[:] = 0.0
assert query.n_finite_points == 2
```

## Bottleneck batch

```python
targets = [
    [[0.0, 1.1], [0.30, 0.80]],
    [[0.0, 2.0]],
    np.empty((0, 2)),
]

distances = topp.bottleneck_distances(query, targets)
assert distances.shape == (3,)
assert distances.dtype == np.float64
```

## Wasserstein batch

```python
w2 = topp.wasserstein_distances(
    query,
    targets,
    order=2,
    internal_p=2,
)
```

The native call reuses a Wasserstein workspace across the target list.

## Caller-provided output

```python
out = np.empty(len(targets), dtype=np.float64)
returned = topp.bottleneck_distances(query, targets, out=out)

assert returned is out
```

`out` must be writable, memory-aligned, C-contiguous, `float64`, and have shape `(len(targets),)`:

```python
bad = np.empty(len(targets), dtype=np.float32)

try:
    topp.bottleneck_distances(query, targets, out=bad)
except TypeError as error:
    print(error)  # out must have dtype float64
```

An empty target list returns an empty `float64` array.
