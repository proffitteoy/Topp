# Topp 使用说明

[English](USAGE.en.md) · [README](../README.md) · [API](API.md) · [数学约定](MATHEMATICS.md)

## 安装

```powershell
py -m pip install topp
```

预编译 wheel 支持 Windows x64 与 CPython 3.10–3.14。从源码安装需要 CMake 3.24+ 和 C++20 编译器：

```powershell
py -m pip install -v .
```

Linux 和 macOS 尚未纳入 CI；从 sdist 构建在这些平台上属于未验证路径。

## 单对距离

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.3, 0.8]])
y = np.array([[0.0, 1.1], [0.4, 0.9]])

d = topp.bottleneck_distance(x, y)
w1 = topp.wasserstein_distance(x, y)
w2 = topp.wasserstein_distance(x, y, order=2, internal_p=2)
```

## PreparedDiagram 与批量计算

反复使用同一 query 时先预处理。`PreparedDiagram` 不可变，可供多个线程并发读取。

```python
query = topp.prepare_diagram(x)
targets = [y, [[0.0, 2.0]], np.empty((0, 2))]

distances = topp.bottleneck_distances(query, targets)
wasserstein = topp.wasserstein_distances(query, targets, order=2, internal_p=2)
```

复用已有输出数组：

```python
out = np.empty(len(targets), dtype=np.float64)
assert topp.bottleneck_distances(query, targets, out=out) is out
```

`out` 必须是可写、C-contiguous、`float64`、shape `(len(targets),)`。

## Threshold decision

只需判断距离是否不大于阈值时：

```python
if topp.bottleneck_within(x, y, 0.1):
    print("close")
```

阈值必须非负且不能为 NaN；`+inf` 合法。

## 输入格式

- shape 必须为 `(n, 2)`，每行是 `[birth, death]`；
- list、tuple、NumPy 数组和非连续 view 均可；
- 内部统一为 C-contiguous `float64`；
- 空图可用 `[]` 或 `np.empty((0, 2))`；
- 支持 `(finite, +inf)`、`(-inf, finite)`、`(-inf, +inf)`；
- 非法输入抛异常，不会被交换、删除或静默修正。

可直接运行 [examples/basic.py](../examples/basic.py) 检查安装。
