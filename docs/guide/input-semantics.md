# 输入、重复点与 essential points

[English](https://proffitteoy.github.io/Topp/en/guide/input-semantics.html)

## 接受的 array-like 输入

每个 diagram 都必须可以转换为形状为 `(n, 2)` 的 `float64` 数组：

```python
import numpy as np
import topp

topp.bottleneck_distance([[0, 1], [1, 3]], np.array([[0.0, 2.0]]))
```

列表、元组、NumPy 数组、整数数组和非连续视图均可使用。Topp 会验证输入，并复制到按 C 顺序连续的 `float64` 存储中。

## 空 diagrams

两种常见写法都可使用：

```python
empty_a = []
empty_b = np.empty((0, 2))

assert topp.bottleneck_distance(empty_a, empty_b) == 0.0
```

## 对角点与重复点

有限对角点 `(a, a)` 的贡献为零，距离计算会忽略它。重复的非对角行保留重数：

```python
duplicates = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]
assert topp.bottleneck_distance(duplicates, duplicates) == 0.0
```

## Essential points

支持三种形式：

```python
import math

positive = [1.0, math.inf]
negative = [-math.inf, 2.0]
fully = [-math.inf, math.inf]
```

Essential points 只能匹配相同类型：

```python
a = [[1.0, math.inf], [-math.inf, 2.0], [-math.inf, math.inf]]
b = [[2.0, math.inf], [-math.inf, 4.0], [-math.inf, math.inf]]

assert topp.bottleneck_distance(a, b) == 2.0
assert topp.wasserstein_distance(a, b) == 3.0
```

任一 essential 类型的重数不同时，结果为 `inf`。

## 非法输入

NaN、`birth > death`、`birth=+inf`、`death=-inf` 和其他非法无穷形式会触发 `ValueError`。Topp 不会交换坐标，也不会静默删除非法行。

```pycon
>>> topp.prepare_diagram([[2.0, 1.0]])
Traceback (most recent call last):
...
ValueError: diagram points must satisfy birth <= death
```

匹配规则见[数学约定](../MATHEMATICS.md#essential-points)。
