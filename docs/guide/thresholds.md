# 精确阈值判断

[English](https://proffitteoy.github.io/Topp/en/guide/thresholds.html)

下游问题是距离是否不大于固定阈值时，请使用 `bottleneck_within`。

```python
import topp

x = [[0.0, 1.0], [0.25, 0.75]]
y = [[0.0, 1.1], [0.30, 0.80]]

if topp.bottleneck_within(x, y, 0.1):
    print("接受候选项")
```

当且仅当下式成立时，函数返回 `True`：

$$
d_B(X,Y) \leq t.
$$

边界包含等号：

```python
distance = topp.bottleneck_distance(x, y)
assert topp.bottleneck_within(x, y, distance)
```

`+inf` 是合法阈值；负数和 NaN 会被拒绝：

```pycon
>>> topp.bottleneck_within(x, y, -1)
Traceback (most recent call last):
...
ValueError: threshold must be non-negative and not NaN
```

0.1.0 版本没有 `wasserstein_within` API。
