# 五分钟快速上手

[English](https://proffitteoy.github.io/Topp/en/getting-started/quickstart.html)

本页使用两个小型 diagram 演示全部公开操作。

## 定义 diagrams

每一行都是 `[birth, death]`：

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.25, 0.75]])
y = np.array([[0.0, 1.1], [0.30, 0.80]])
```

## 成对距离

```python
bottleneck = topp.bottleneck_distance(x, y)
w1 = topp.wasserstein_distance(x, y)
w2 = topp.wasserstein_distance(x, y, order=2, internal_p=2)

print(bottleneck, w1, w2)
```

`wasserstein_distance` 默认计算精确的 $W_{1,\infty}$。传入 `order=2, internal_p=2` 可选择精确的 $W_{2,2}$。

## 预处理重复使用的查询

```python
query = topp.prepare_diagram(x)

print(query.n_points)
print(query.n_finite_points)
print(topp.bottleneck_distance(query, y))
```

`PreparedDiagram` 持有经过验证的输入副本。之后修改 `x` 不会改变 `query`。

## 一个查询，多个目标

```python
targets = [y, [[0.0, 2.0]], np.empty((0, 2))]

bottleneck_batch = topp.bottleneck_distances(query, targets)
wasserstein_batch = topp.wasserstein_distances(
    query, targets, order=2, internal_p=2
)
```

两个结果都是一维 `float64` 数组，每个目标对应一个元素。

## 复用输出内存

```python
out = np.empty(len(targets), dtype=np.float64)
returned = topp.bottleneck_distances(query, targets, out=out)

assert returned is out
```

## 精确阈值判断

```python
if topp.bottleneck_within(query, y, 0.1):
    print("距离不大于 0.1")
```

该调用直接判断 $d_B(X,Y) \leq 0.1$，不会先计算近似距离再比较。

## 后续阅读

- [Bottleneck 距离](../guide/bottleneck.md)
- [Wasserstein 距离](../guide/wasserstein.md)
- [预处理 diagrams 与批量调用](../guide/prepared-batch.md)
- [输入与 essential points 语义](../guide/input-semantics.md)
