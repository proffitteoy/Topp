# 项目范围：独立 Python 标准库

记录日期：2026-08-10

## 产品目标

本项目交付一个可通过 Python 包管理工具安装、具有稳定公共 API 和预编译 wheel 的 Bottleneck Distance 库。核心求解器可以使用 C++ 等原生实现，但用户面对的是普通 Python/NumPy 接口。

GUDHI 在本项目中承担三种角色：

1. `e=0` exact 结果的正确性 oracle；
2. 单对计算的性能基线；
3. 算法和历史补丁的参考来源。

GUDHI 不是运行时 API 边界。本项目不采用任何特定下游应用的阈值、排序、数据分层或工作流约定。

## v0 能力边界

首阶段聚焦 exact Bottleneck Distance：

- 接受形状为 `(n, 2)`、可转换为 `float64` 的 Python/NumPy 输入；
- 正确定义空图、对角点、重复点、有限点和无穷点行为；
- 单对计算接口；
- 可复用的 prepared diagram，用于避免批量比较时重复预处理；
- one-to-many 批量接口；
- 原生计算期间释放 GIL；
- 结果与固定版本 GUDHI `e=0` oracle 严格差分。

many-to-many、并行调度、通用 threshold-decision 和近似算法属于后续能力，但架构上不能阻塞它们。

## 非目标

- 不复刻完整 GUDHI 包；
- 不绑定任何应用专用的固定阈值、排序或数据维度规则；
- 不要求用户理解 GUDHI 的投影图或底层几何数据结构；
- 不以某一台机器上的专用 wheel 作为标准发行结果；
- 不用内核 microbenchmark 代替 Python API 或批量吞吐结论。

## 公共 API 草案

包名尚未确定，以下只描述形状：

```python
distance = bottleneck_distance(diagram_a, diagram_b)

prepared = prepare_diagram(diagram_a)
distance = bottleneck_distance(prepared, diagram_b)

distances = bottleneck_distances(prepared, diagrams)
```

实验算法、dispatcher 选择和诊断计数器默认隐藏；可通过独立的 benchmark/debug 接口暴露，但不应成为普通用户必须配置的参数。

## 发布约束

目标平台、Python 版本和包名仍待确认。最低发布要求包括：

- 源码发行包；
- 主流 CPython 与 Windows/Linux/macOS wheel；
- 可重复的 Release 构建；
- API、类型、异常和数值语义文档；
- 单元测试、随机差分、边界输入和性能回归基线。

## 仍需冻结的决策

- 发行包名和顶层导入名；
- 支持的 Python/NumPy 版本与平台；
- exact 的“严格一致”是相同数学值还是要求浮点逐位一致；
- 无效输入、`NaN`、反向区间和各类无穷点的异常策略；
- v0 是否同时发布 one-to-many，还是先以内测 API 进入；
- C++ 绑定与构建方案。
