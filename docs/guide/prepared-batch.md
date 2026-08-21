# 预处理 diagrams 与批量调用

[English](https://proffitteoy.github.io/Topp/en/guide/prepared-batch.html)

需要反复比较同一个 diagram 时，请使用预处理和 one-to-many API。

## 只预处理一次

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.25, 0.75]])
query = topp.prepare_diagram(x)

print(query.n_points)         # 所有合法输入行
print(query.n_finite_points)  # 对角线上方的有限点
```

预处理对象不可变，并持有输入副本：

```python
x[:] = 0.0
assert query.n_finite_points == 2
```

## Bottleneck 批量调用

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

## Wasserstein 批量调用

```python
w2 = topp.wasserstein_distances(
    query,
    targets,
    order=2,
    internal_p=2,
)
```

原生调用会在整个目标列表之间复用 Wasserstein workspace。

## 调用方提供输出数组

```python
out = np.empty(len(targets), dtype=np.float64)
returned = topp.bottleneck_distances(query, targets, out=out)

assert returned is out
```

`out` 必须可写、内存对齐、按 C 顺序连续、类型为 `float64`，且形状为 `(len(targets),)`：

```python
bad = np.empty(len(targets), dtype=np.float32)

try:
    topp.bottleneck_distances(query, targets, out=bad)
except TypeError as error:
    print(error)  # out must have dtype float64
```

目标列表为空时，返回空的 `float64` 数组。
