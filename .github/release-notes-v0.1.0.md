# Topp v0.1.0 内测版

Topp 的首个公开 Python 版本，提供 exact persistence-diagram Bottleneck 与 Wasserstein 距离。

## 主要能力

- `bottleneck_distance`、`prepare_diagram`、`bottleneck_distances` 与 `bottleneck_within`；
- exact `W1-L∞` 和 `W2-L2` Wasserstein single/prepared/batch API；
- CPython 3.10–3.14 Windows x64 wheels 和源码包；
- 严格输入校验、类型提示、GIL 释放和不可变 `PreparedDiagram`。

## 验证

- C++ Bottleneck 与 Wasserstein 完整测试；
- Python API、输入契约、prepared/batch/out 一致性测试；
- Bottleneck 对 GUDHI 3.13.0 `e=0` 差分；
- Wasserstein 独立小规模 brute-force 差分；
- AVX2 runtime dispatch 与无 AVX2 scalar 构建；
- 每个 wheel 从 TestPyPI 安装后运行示例。

## 已知限制

- 首发只提供 Windows x64 wheels；其他平台需从 sdist 构建。
- Wasserstein 只支持 `(order=1, internal_p=∞)` 和 `(order=2, internal_p=2)`。
- Wasserstein 内核仍会继续优化；实验策略不属于公开 Python API。
- GitHub 将本版本标记为 Prerelease，但 PyPI 会把 `0.1.0` 视为正式版本号。
