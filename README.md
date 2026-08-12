# Topp

[English](README.en.md) · [使用说明](docs/USAGE.md) · [API](docs/API.md) · [开发指南](docs/DEVELOPMENT.md)

Topp 是一个精简的 persistence diagram 距离包。它只专注于 exact Bottleneck 和两种 exact Wasserstein 距离，并由自适应 C++20 内核完成计算。

> **v0.1.0 内测版：** PyPI 会把 `0.1.0` 视为正式版本号，但本项目仍处于公开内测阶段；Python API 已冻结，Wasserstein 内核仍会继续优化。

## 特性

- exact Bottleneck Distance（点间使用 `L∞`）；
- exact `W1-L∞` 与 `W2-L2` Wasserstein Distance；
- 不可变的 `PreparedDiagram`；
- 原生 one-to-many 批量计算与可复用输出数组；
- exact threshold decision：`bottleneck_within`；
- Windows x64 的 CPython 3.10–3.14 wheels；
- 运行时仅依赖 NumPy；AVX2 路径在运行时检测，不要求所有机器支持 AVX2。

## 安装

```powershell
py -m pip install topp
```

首发 wheel 面向 Windows x64。其他平台可以从 sdist 构建，但需要 CMake 3.24+ 和 C++20 编译器。

## 快速开始

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.3, 0.8]])
y = np.array([[0.0, 1.1], [0.4, 0.9]])

print(topp.bottleneck_distance(x, y))
print(topp.wasserstein_distance(x, y, order=2, internal_p=2))

query = topp.prepare_diagram(x)
print(topp.bottleneck_distances(query, [y, np.empty((0, 2))]))
print(topp.bottleneck_within(query, y, 0.1))
```

完整示例见 [examples/basic.py](examples/basic.py)。

## 支持的度量

| 函数 | 语义 | 状态 |
|---|---|---|
| `bottleneck_distance` | exact Bottleneck，内部 `L∞` | 支持 |
| `wasserstein_distance(..., order=1, internal_p=np.inf)` | exact `W1-L∞` | 支持 |
| `wasserstein_distance(..., order=2, internal_p=2)` | exact `W2-L2` | 支持 |
| 其他 Wasserstein 参数 | 数学上可能合法 | `NotImplementedError` |

## 输入契约

输入必须可转换为 `(n, 2)` 的 `float64` 数组。空图、对角点、重复点和规范 essential points 合法。NaN、`birth > death`、`birth=+inf`、`death=-inf` 及其他非法无穷组合会抛出 `ValueError`，不会被静默修正。

详见 [API 文档](docs/API.md)。

## 实验功能

C++ 源码保留候选生成、图表示、matching、component 和 incremental pricing 等实验策略，供维护者复现与比较。它们不会暴露到普通 Python API，也不代表默认性能承诺。研究证据和历史方案见 [docs/research](docs/research/README.md)。

## 开发

```powershell
py -m pip install -v .
py -m pytest tests/python
cmd.exe /d /c scripts\build-kernel.cmd
```

现有 `include/bottleneck/*` C++ 接口用于社区维护和内核实验，不承诺稳定 ABI。构建、测试和 benchmark 约定见 [开发指南](docs/DEVELOPMENT.md)。

## 引用

研究中使用 Topp 时，请引用仓库版本与发布标签。机器可读元数据见 [CITATION.cff](CITATION.cff)。

## 许可

Topp 使用 [MIT License](LICENSE)。GUDHI 仅作为测试 oracle、语义参考及历史补丁来源，不是运行时依赖；详情见 [第三方声明](THIRD_PARTY_NOTICES.md)。
