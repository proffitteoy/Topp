# Topp

```{toctree}
:hidden:
:caption: 快速开始

getting-started/installation
getting-started/quickstart
```

```{toctree}
:hidden:
:caption: 距离计算

guide/bottleneck
guide/wasserstein
guide/prepared-batch
guide/thresholds
guide/input-semantics
guide/threads
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
USAGE
English <https://proffitteoy.github.io/Topp/en/>
更新日志 <https://github.com/proffitteoy/Topp/blob/main/CHANGELOG.md>
源代码 <https://github.com/proffitteoy/Topp>
```

[![CI](https://github.com/proffitteoy/Topp/actions/workflows/ci.yml/badge.svg)](https://github.com/proffitteoy/Topp/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/topp)](https://pypi.org/project/topp/)
[![Python](https://img.shields.io/pypi/pyversions/topp)](https://pypi.org/project/topp/)

[English](https://proffitteoy.github.io/Topp/en/)

Topp 是一个专注于计算 **persistence diagram 精确距离**的 Python 包。它面向已经从 GUDHI、Ripser 或其他流程获得 diagram，并需要在 Python 中进行严格、重复比较的用户。

Topp 目前提供：

- 点间使用 $L^\infty$ 的精确 Bottleneck 距离；
- 精确的 $W_{1,\infty}$ 与 $W_{2,2}$ Wasserstein 距离；
- 不可变的预处理 diagram 和原生 one-to-many 调用；
- 可复用的输出数组和精确的 Bottleneck 阈值判断。

计算内核使用 C++20 编写，NumPy 是唯一的运行时依赖。

```{important}
Topp 只比较 persistence diagrams，不构造 filtration，也不计算 persistent homology。它不能替代完整的 TDA 库。
```

## 安装

```console
py -m pip install topp
```

预编译 wheel 目前支持 Windows x64 和 CPython 3.10–3.14。在其他平台尝试源码构建前，请阅读[平台与支持](platforms.md)。

## 第一次计算

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.3, 0.8]])
y = np.array([[0.0, 1.1], [0.4, 0.9]])

print(topp.bottleneck_distance(x, y))
print(topp.wasserstein_distance(x, y, order=2, internal_p=2))
```

接下来可以阅读[五分钟快速上手](getting-started/quickstart.md)，或从左侧导航选择具体主题。

## Python 批量距离性能

下表是 1 个查询图与 64 个目标图的批量中位耗时，单位为毫秒；每个规模合并 uniform、near-diagonal、clustered、duplicate-heavy 和 separated 五类合成输入的全部计时轮次后取中位数，数值越低越快。

| 距离 | 每图点数 | Topp | GUDHI | giotto-tda | Topp 相对 GUDHI | Topp 相对 giotto-tda |
|---|---:|---:|---:|---:|---:|---:|
| Bottleneck | 8 | **1.059** | 10.729 | 14.434 | **10.13×** | **13.62×** |
| Bottleneck | 32 | **4.547** | 43.860 | 90.804 | **9.65×** | **19.97×** |
| Bottleneck | 128 | **16.632** | 235.851 | 558.250 | **14.18×** | **33.56×** |
| Bottleneck | 512 | **236.045** | 1,882.070 | 3,639.137 | **7.97×** | **15.42×** |
| Wasserstein | 8 | **1.104** | 15.344 | 11.758 | **13.90×** | **10.65×** |
| Wasserstein | 32 | **3.246** | 21.727 | 114.656 | **6.69×** | **35.33×** |
| Wasserstein | 128 | **10.101** | 106.639 | 1,000.598 | **10.56×** | **99.06×** |
| Wasserstein | 512 | **12.689** | 2,986.669 | 8,808.310 | **235.37×** | **694.17×** |

测试于 2026-08-17，环境为 Windows 11、Python 3.12.13、单线程，Topp 0.1.0 当前主线 MSVC wheel（SHA-256 `dffa357a504121538d63e3fea3675054430f67e04f91c35828f29a18d150f641`）、GUDHI 3.13.0、giotto-tda 0.6.2。Bottleneck 表比较各库默认调用（Topp 为 exact，GUDHI `e=None` 和 giotto-tda 为 approximate）；Wasserstein 中 Topp/GUDHI 为 exact $W_{1,\infty}$，giotto-tda 默认是 approximate $W_2$，因此其数字只表示默认 Python API 的调用速度，不是同一数学任务的算法排名。不同分布的耗时差异可能很大，汇总值不表示每类输入都达到相同倍数。
