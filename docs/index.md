# Topp

```{toctree}
:hidden:
:caption: 快速开始

getting-started/installation
getting-started/quickstart
```

```{toctree}
:hidden:
:caption: 使用指南

guide/bottleneck
guide/wasserstein
guide/prepared-batch
guide/thresholds
guide/input-semantics
guide/threads
USAGE
```

```{toctree}
:hidden:
:caption: 参考

API
MATHEMATICS
platforms
```

```{toctree}
:hidden:
:caption: 项目

DEVELOPMENT
更新日志 <https://github.com/proffitteoy/Topp/blob/main/CHANGELOG.md>
源代码 <https://github.com/proffitteoy/Topp>
English <https://proffitteoy.github.io/Topp/en/>
```

[![CI](https://github.com/proffitteoy/Topp/actions/workflows/ci.yml/badge.svg)](https://github.com/proffitteoy/Topp/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/topp)](https://pypi.org/project/topp/)
[![Python](https://img.shields.io/pypi/pyversions/topp)](https://pypi.org/project/topp/)

[English](https://proffitteoy.github.io/Topp/en/)

Topp 是一个用于计算 persistence diagram 距离的 Python 库。它实现精确 Bottleneck distance、$W_{1,\infty}$ Wasserstein distance 和 $W_{2,2}$ Wasserstein distance，并提供面向重复比较的预处理和批量接口。

计算内核使用 C++20 实现，Python 接口以 NumPy 数组作为主要输入形式。Topp 不构造 filtration，也不计算 persistent homology。

## 安装

```console
pip install topp
```

PyPI 提供 Windows x64 和 Linux x86_64 的 CPython 3.10–3.14 wheels。macOS 源码构建尚未纳入持续集成测试。

详细说明见[安装](getting-started/installation.md)。

## 第一次计算

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.3, 0.8]])
y = np.array([[0.0, 1.1], [0.4, 0.9]])

bottleneck = topp.bottleneck_distance(x, y)
w1 = topp.wasserstein_distance(x, y)
w2 = topp.wasserstein_distance(x, y, order=2, internal_p=2)
```

每个 persistence diagram 使用 shape 为 `(n, 2)` 的数组表示。完整示例见[五分钟快速上手](getting-started/quickstart.md)。

## 支持的距离

| 距离 | Python API |
|---|---|
| Bottleneck, $L^\infty$ | `topp.bottleneck_distance` |
| $W_{1,\infty}$ | `topp.wasserstein_distance(..., order=1, internal_p=np.inf)` |
| $W_{2,2}$ | `topp.wasserstein_distance(..., order=2, internal_p=2)` |

上述距离均按 Topp 文档中的数学约定精确计算。其他 Wasserstein 参数组合目前不实现。

距离定义和对角线匹配规则见[数学约定](MATHEMATICS.md)。

## 重复比较

`PreparedDiagram` 保存经过验证和预处理的 persistence diagram。一个 diagram 需要与多个目标重复比较时，可以复用该对象：

```python
query = topp.prepare_diagram(x)
targets = [y, [[0.0, 2.0]], np.empty((0, 2))]

distances = topp.bottleneck_distances(query, targets)
```

对应的 Wasserstein 批量接口为 `topp.wasserstein_distances`。只需要判断 Bottleneck distance 是否不超过某个阈值时，可以使用 `topp.bottleneck_within`。

相关说明见[预处理与批量计算](guide/prepared-batch.md)和[阈值判断](guide/thresholds.md)。

## 文档结构

文档按用途分为以下几部分：

- **快速开始**：安装和最小可运行示例；
- **使用指南**：各距离函数、批量计算、输入语义和线程行为；
- **参考**：Python API、数学约定和平台支持；
- **项目**：构建、测试、开发和版本记录。

API 的调用契约见 [Python API](API.md)。输入中的对角点、重复点、essential points 和无穷值处理见[输入语义](guide/input-semantics.md)。

## 性能

Topp 针对 persistence-diagram distance 的重复计算实现了专门的 C++ 内核。仓库中的 benchmark 对 Topp、GUDHI 和 Hera 的 Python 接口进行了固定环境下的比较。

完整结果、数据分布、软件版本、近似参数和测试方法见仓库的 [benchmark documentation](https://github.com/proffitteoy/Topp/tree/main/benchmarks)。不同库的默认算法和精度设置并不完全相同，性能数字应与实验配置一起解释。

## 开发与引用

源码构建、测试和 benchmark 说明见[开发指南](DEVELOPMENT.md)。

研究工作中使用 Topp 时，请引用实际使用的软件版本；机器可读的引用信息位于仓库的 `CITATION.cff`。
