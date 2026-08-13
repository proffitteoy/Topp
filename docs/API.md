# Topp Python API

[English](https://proffitteoy.github.io/Topp/en/API.html) · [使用说明](USAGE.md) · [数学约定](MATHEMATICS.md)

所有公开对象都从 `topp` 顶层导入。`topp._core` 是私有实现，不受兼容性承诺保护。

## 输入类型

`DiagramLike = ArrayLike | PreparedDiagram`。array-like 必须能转换为 shape `(n, 2)` 的 `float64` 数组。严格输入规则见[使用说明](USAGE.md)。

## `PreparedDiagram`

由 `prepare_diagram` 返回的不可变预处理对象。

- `n_points`：合法原始输入的行数，包括对角点与 essential points；
- `n_finite_points`：严格位于对角线上方的有限点数；
- `len(prepared)`：等于 `n_points`。

## `prepare_diagram(diagram) -> PreparedDiagram`

校验并复制输入，构建批量复用所需缓存。输入之后的修改不会改变 prepared 对象。

## `bottleneck_distance(diagram_a, diagram_b) -> float`

返回 exact Bottleneck Distance，点间使用 `L∞`，对角线具有无限重数。可能返回 `np.inf`。

## `wasserstein_distance(diagram_a, diagram_b, *, order=1, internal_p=np.inf) -> float`

支持：

- `order=1, internal_p=np.inf`：exact `W1-L∞`；
- `order=2, internal_p=2`：exact `W2-L2`。

其他组合抛 `NotImplementedError`。

## `bottleneck_distances(query, diagrams, *, out=None) -> np.ndarray`

一次 native 调用计算 one-to-many exact Bottleneck 距离。返回 shape `(m,)`、dtype `float64`。提供合法 `out` 时写入并返回同一数组。

## `wasserstein_distances(query, diagrams, *, order=1, internal_p=np.inf, out=None) -> np.ndarray`

一次 native 调用计算 one-to-many Wasserstein 距离，并在批内复用 workspace。参数与单对函数相同。

## `bottleneck_within(diagram_a, diagram_b, threshold) -> bool`

exact decision API：仅当 Bottleneck Distance `<= threshold` 时返回 `True`。阈值必须非负且非 NaN；`+inf` 合法。

## 异常

- `TypeError`：不可转换的输入、非法 `out` 类型、非数值 threshold；
- `ValueError`：shape、点语义、threshold 或 `out` 布局非法；
- `NotImplementedError`：暂不支持的 Wasserstein 参数；
- `MemoryError`：原生分配失败。

## 线程与 GIL

原生距离计算和批量循环释放 Python GIL。`PreparedDiagram` 不可变，可在多个线程中并发读取；`out` 的并发写入同步由调用方负责。
