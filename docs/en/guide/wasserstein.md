# Wasserstein distance

[中文](https://proffitteoy.github.io/Topp/guide/wasserstein.html)

Topp implements exactly two Wasserstein conventions in version 0.1.0.

## $W_{1,\infty}$

This is the default:

```python
import topp

x = [[0.0, 1.0], [0.25, 0.75]]
y = [[0.0, 1.1], [0.30, 0.80]]

w1 = topp.wasserstein_distance(x, y)
```

It is equivalent to spelling out the parameters:

```python
import numpy as np

same_w1 = topp.wasserstein_distance(
    x, y, order=1, internal_p=np.inf
)
assert same_w1 == w1
```

## $W_{2,2}$

```python
w2 = topp.wasserstein_distance(
    x, y, order=2, internal_p=2
)
```

For $W_{2,2}$, a finite point $x=(b,d)$ has Euclidean distance

$$
d_2(x,\Delta)=\frac{d-b}{\sqrt{2}}
$$

to the diagonal.

## Unsupported parameters

Other combinations are rejected rather than silently reinterpreted:

```pycon
>>> topp.wasserstein_distance(x, y, order=3, internal_p=2)
Traceback (most recent call last):
...
NotImplementedError: topp 0.1.0 supports only (order=1, internal_p=inf) and (order=2, internal_p=2)
```

See [Mathematical conventions](../MATHEMATICS.md#wasserstein-distance) for the definition of `order` and `internal_p`.
