# Inputs, duplicates, and essential points

[中文](https://proffitteoy.github.io/Topp/guide/input-semantics.html)

## Accepted array-like inputs

Every diagram must convert to a `float64` array of shape `(n, 2)`:

```python
import numpy as np
import topp

topp.bottleneck_distance([[0, 1], [1, 3]], np.array([[0.0, 2.0]]))
```

Lists, tuples, NumPy arrays, integer arrays, and non-contiguous views are accepted. Topp validates and copies them into C-contiguous `float64` storage.

## Empty diagrams

Both common spellings are accepted:

```python
empty_a = []
empty_b = np.empty((0, 2))

assert topp.bottleneck_distance(empty_a, empty_b) == 0.0
```

## Diagonal and duplicate points

Finite diagonal points `(a, a)` contribute zero and are ignored by the distance calculation. Duplicate off-diagonal rows retain multiplicity:

```python
duplicates = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]
assert topp.bottleneck_distance(duplicates, duplicates) == 0.0
```

## Essential points

Three forms are supported:

```python
import math

positive = [1.0, math.inf]
negative = [-math.inf, 2.0]
fully = [-math.inf, math.inf]
```

Essential points match only the same type:

```python
a = [[1.0, math.inf], [-math.inf, 2.0], [-math.inf, math.inf]]
b = [[2.0, math.inf], [-math.inf, 4.0], [-math.inf, math.inf]]

assert topp.bottleneck_distance(a, b) == 2.0
assert topp.wasserstein_distance(a, b) == 3.0
```

If the multiplicity of any essential type differs, the result is `inf`.

## Invalid inputs

Topp raises `ValueError` for NaN, `birth > death`, `birth=+inf`, `death=-inf`, and other invalid infinity forms. It never swaps coordinates or silently deletes invalid rows.

```pycon
>>> topp.prepare_diagram([[2.0, 1.0]])
Traceback (most recent call last):
...
ValueError: diagram points must satisfy birth <= death
```

See [Mathematical conventions](../MATHEMATICS.md#essential-points) for the matching rules.
