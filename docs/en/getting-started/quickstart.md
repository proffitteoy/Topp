# Five-minute quickstart

[中文](https://proffitteoy.github.io/Topp/getting-started/quickstart.html)

This page runs through every public operation using two small diagrams.

## Define diagrams

Each row is `[birth, death]`:

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.25, 0.75]])
y = np.array([[0.0, 1.1], [0.30, 0.80]])
```

## Pairwise distances

```python
bottleneck = topp.bottleneck_distance(x, y)
w1 = topp.wasserstein_distance(x, y)
w2 = topp.wasserstein_distance(x, y, order=2, internal_p=2)

print(bottleneck, w1, w2)
```

`wasserstein_distance` defaults to exact $W_{1,\infty}$. Passing `order=2, internal_p=2` selects exact $W_{2,2}$.

## Prepare a repeated query

```python
query = topp.prepare_diagram(x)

print(query.n_points)
print(query.n_finite_points)
print(topp.bottleneck_distance(query, y))
```

`PreparedDiagram` owns a validated copy of the input. Mutating `x` later does not change `query`.

## One query, many targets

```python
targets = [y, [[0.0, 2.0]], np.empty((0, 2))]

bottleneck_batch = topp.bottleneck_distances(query, targets)
wasserstein_batch = topp.wasserstein_distances(
    query, targets, order=2, internal_p=2
)
```

Both results are one-dimensional `float64` arrays with one entry per target.

## Reuse output memory

```python
out = np.empty(len(targets), dtype=np.float64)
returned = topp.bottleneck_distances(query, targets, out=out)

assert returned is out
```

## Exact threshold decision

```python
if topp.bottleneck_within(query, y, 0.1):
    print("distance is at most 0.1")
```

This directly decides $d_B(X,Y) \leq 0.1$; it does not compare against an approximate distance.

## Next steps

- [Bottleneck distance](../guide/bottleneck.md)
- [Wasserstein distance](../guide/wasserstein.md)
- [Prepared diagrams and batches](../guide/prepared-batch.md)
- [Input and essential-point semantics](../guide/input-semantics.md)
