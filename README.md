<p align="center">
  <img src="https://raw.githubusercontent.com/proffitteoy/Topp/main/assets/topp-mark.svg" width="96" alt="Topp logo">
</p>

<h1 align="center">Topp</h1>

<p align="center">用于计算 persistence diagram 精确距离的 Python 库。</p>

<p align="center">
  <a href="https://github.com/proffitteoy/Topp/actions/workflows/ci.yml"><img src="https://github.com/proffitteoy/Topp/actions/workflows/ci.yml/badge.svg" alt="CI status"></a>
  <a href="https://pypi.org/project/topp/"><img src="https://img.shields.io/pypi/v/topp" alt="PyPI version"></a>
  <a href="https://pypi.org/project/topp/"><img src="https://img.shields.io/pypi/pyversions/topp?logo=python&logoColor=white" alt="Python versions"></a>
</p>

<p align="center">
  <a href="README.en.md">English</a> ·
  <a href="https://proffitteoy.github.io/Topp/">Documentation</a> ·
  <a href="https://pypi.org/project/topp/">PyPI</a> ·
  <a href="CHANGELOG.md">Changelog</a>
</p>

Topp 计算 persistence diagrams 之间的精确 Bottleneck 与 Wasserstein 距离。Python 接口基于 NumPy，计算内核使用 C++20 实现。

目前支持：

- Bottleneck distance，点间距离为 $L^\infty$；
- $W_{1,\infty}$ Wasserstein distance；
- $W_{2,2}$ Wasserstein distance；
- persistence diagram 预处理与 one-to-many 批量计算；
- Bottleneck distance 的精确阈值判断；
- 原生计算期间释放 Python GIL。

Topp 只计算已经给定的 persistence diagrams 之间的距离，不构造 filtration，也不计算 persistent homology。

## 安装

```bash
pip install topp
```

PyPI 提供 Windows x64 和 Linux x86_64 的 CPython 3.10–3.14 wheels。运行时依赖仅包含 NumPy。

从源码构建需要 CMake 3.24+ 和支持 C++20 的编译器。macOS 尚未纳入 CI。

## 基本用法

Persistence diagram 使用 shape 为 `(n, 2)` 的数组表示，每一行对应一个 `[birth, death]` 点。

```python
import numpy as np
import topp

x = np.array([
    [0.0, 1.0],
    [0.3, 0.8],
])
y = np.array([
    [0.0, 1.1],
    [0.4, 0.9],
])

bottleneck = topp.bottleneck_distance(x, y)
w1 = topp.wasserstein_distance(x, y)
w2 = topp.wasserstein_distance(x, y, order=2, internal_p=2)
```

`wasserstein_distance` 默认计算精确的 $W_{1,\infty}$ 距离。当前公开 API 只支持下列三种距离：

| 距离 | 调用 |
|---|---|
| Bottleneck, $L^\infty$ | `topp.bottleneck_distance(x, y)` |
| $W_{1,\infty}$ | `topp.wasserstein_distance(x, y)` |
| $W_{2,2}$ | `topp.wasserstein_distance(x, y, order=2, internal_p=2)` |

其他 `(order, internal_p)` 组合目前不实现。

## 批量计算

同一个 diagram 需要重复参与计算时，可以先构造 `PreparedDiagram`：

```python
query = topp.prepare_diagram(x)
targets = [
    y,
    np.array([[0.0, 2.0]]),
    np.empty((0, 2)),
]

distances = topp.bottleneck_distances(query, targets)
```

Wasserstein distance 提供对应的 `topp.wasserstein_distances` 批量接口。批量接口可以通过 `out` 参数复用已有的 `float64` 输出数组。

只需要判断 Bottleneck distance 是否不超过给定阈值时，可以直接使用：

```python
close = topp.bottleneck_within(x, y, 0.1)
```

## 输入与数学约定

有限 persistence point 必须满足 `birth <= death`。空 diagram、重复点、对角点和支持的 essential points 都可以作为输入；非法无穷值组合和 NaN 会被拒绝。

距离定义、对角线匹配、重复点和 essential points 的处理见 [Mathematical conventions](docs/MATHEMATICS.md)。完整输入规则见 [Input semantics](docs/guide/input-semantics.md)。

## 性能

Topp 的内核针对 persistence-diagram distance 中的候选边生成、matching、重复点和 repeated-query 场景进行了专门实现。仓库包含与 GUDHI 和 Hera 的可复现实验脚本及结果。

完整的测试环境、调用方式和精度设置见 [benchmark documentation](benchmarks/README.md)。不同库的默认算法和近似参数并不完全相同，因此性能结果应与对应实验设置一起解释。

## 文档

完整文档位于 <https://proffitteoy.github.io/Topp/>。

主要入口：

- [Installation](docs/getting-started/installation.md)
- [Quickstart](docs/getting-started/quickstart.md)
- [Bottleneck distance](docs/guide/bottleneck.md)
- [Wasserstein distance](docs/guide/wasserstein.md)
- [Prepared diagrams and batch computation](docs/guide/prepared-batch.md)
- [Python API](docs/API.md)
- [Mathematical conventions](docs/MATHEMATICS.md)
- [Development](docs/DEVELOPMENT.md)

## 开发

安装本地源码并运行 Python 测试：

```bash
pip install -v .
pytest tests/python
```

C++ 内核的构建、测试和 benchmark 说明见 [Development](docs/DEVELOPMENT.md)。公开 Python API 属于 `1.x` 兼容性范围；内部 C++ 接口和 ABI 不作稳定性承诺。

## 引用

研究工作中使用 Topp 时，请引用实际使用的软件版本。机器可读的引用信息见 [CITATION.cff](CITATION.cff)。

## License

Topp 采用 [MIT License](LICENSE)。第三方代码和测试来源说明见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
