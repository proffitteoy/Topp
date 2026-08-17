# Wasserstein 距离

[English](https://proffitteoy.github.io/Topp/en/guide/wasserstein.html)

Topp 当前只实现两种 Wasserstein 约定。

## $W_{1,\infty}$

这是默认值：

```python
import topp

x = [[0.0, 1.0], [0.25, 0.75]]
y = [[0.0, 1.1], [0.30, 0.80]]

w1 = topp.wasserstein_distance(x, y)
```

它等价于显式传入参数：

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

对于 $W_{2,2}$，有限点 $x=(b,d)$ 到对角线的欧氏距离为

$$
d_2(x,\Delta)=\frac{d-b}{\sqrt{2}}.
$$

## 不支持的参数

其他组合会明确拒绝，不会被静默解释为其他含义：

```pycon
>>> topp.wasserstein_distance(x, y, order=3, internal_p=2)
Traceback (most recent call last):
...
NotImplementedError: topp supports only (order=1, internal_p=inf) and (order=2, internal_p=2)
```

`order` 和 `internal_p` 的定义见[数学约定](../MATHEMATICS.md#wasserstein-distance)。
