# Exact threshold decisions

Use `bottleneck_within` when the downstream question is whether a distance is at most a fixed threshold.

```python
import topp

x = [[0.0, 1.0], [0.25, 0.75]]
y = [[0.0, 1.1], [0.30, 0.80]]

if topp.bottleneck_within(x, y, 0.1):
    print("accept candidate")
```

The function returns `True` exactly when

$$
d_B(X,Y) \leq t.
$$

The boundary is inclusive:

```python
distance = topp.bottleneck_distance(x, y)
assert topp.bottleneck_within(x, y, distance)
```

`+inf` is a valid threshold. Negative values and NaN are rejected:

```pycon
>>> topp.bottleneck_within(x, y, -1)
Traceback (most recent call last):
...
ValueError: threshold must be non-negative and not NaN
```

There is no `wasserstein_within` API in version 0.1.0.
