# Bottleneck distance

`bottleneck_distance(X, Y)` computes the exact Bottleneck distance with $L^\infty$ point costs.

## Pairwise call

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.2, 0.7]])
y = np.array([[0.0, 1.1], [0.3, 0.8]])

distance = topp.bottleneck_distance(x, y)
print(distance)
```

The result is a Python `float`. Empty diagrams, diagonal points, duplicate points, and supported essential points are valid inputs.

## Matching to the diagonal

A finite point $x=(b,d)$ may match the diagonal at cost

$$
d_\infty(x,\Delta)=\frac{d-b}{2}.
$$

For example:

```pycon
>>> topp.bottleneck_distance([[0.0, 2.0]], [])
1.0
```

## Symmetry and prepared input

```python
prepared = topp.prepare_diagram(x)

forward = topp.bottleneck_distance(prepared, y)
reverse = topp.bottleneck_distance(y, prepared)
assert forward == reverse
```

Preparing does not change the mathematical result; it avoids repeating input preprocessing when the diagram is reused.

See [Mathematical conventions](../MATHEMATICS.en.md#bottleneck-distance) for the full definition.
