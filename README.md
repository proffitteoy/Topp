<p align="center">
  <img src="https://raw.githubusercontent.com/proffitteoy/Topp/main/assets/topp-mark.svg" width="96" alt="Topp 标志">
</p>

<h1 align="center">Topp</h1>

<p align="center">用于 persistence diagram 的 exact Bottleneck 与 Wasserstein 距离。</p>

<p align="center">
  <a href="https://github.com/proffitteoy/Topp/actions/workflows/ci.yml"><img src="https://github.com/proffitteoy/Topp/actions/workflows/ci.yml/badge.svg" alt="CI status"></a>
  <a href="https://pypi.org/project/topp/"><img src="https://img.shields.io/pypi/v/topp" alt="PyPI version"></a>
  <a href="https://pypi.org/project/topp/"><img src="https://img.shields.io/pypi/pyversions/topp?logo=python&logoColor=white" alt="Python versions"></a>
</p>

<p align="center">
  <a href="README.en.md">English</a> ·
  <a href="https://proffitteoy.github.io/Topp/">网站</a> ·
  <a href="docs/USAGE.md">使用说明</a> ·
  <a href="docs/API.md">API</a> ·
  <a href="docs/MATHEMATICS.md">数学约定</a> ·
  <a href="docs/DEVELOPMENT.md">开发指南</a> ·
  <a href="CHANGELOG.md">更新日志</a>
</p>

Topp 面向**已经拥有 persistence diagrams，需要在 Python 中做严格、重复距离比较**的用户。它提供小型 Python API 和自适应 C++20 内核；对同一个 diagram 执行多次比较时，可预处理一次并直接调用原生批量接口。

> **v1.0.0 稳定版：** 本页记录的 Python API 进入 `1.x` 兼容性范围。C++ 头文件与 ABI 仍是内核维护接口，不属于稳定性承诺。

## 适合什么场景

选择 Topp，如果你：

- 已经从 GUDHI、Ripser 或其他流程获得 persistence diagrams，只需要计算它们之间的距离；
- 需要 exact Bottleneck、`W1-L∞` 或 `W2-L2`，不希望近似参数改变阈值判断；
- 会用一个查询图反复比较许多候选图，希望复用预处理结果、workspace 或输出数组；
- 希望运行时只引入 NumPy，并使用带类型信息的 Python API。

Topp **不负责生成 persistence diagrams**，也不提供任意 `(order, internal_p)`、近似/GPU 距离或完整 TDA 工作流。需要这些能力时，应继续使用覆盖面更广的 TDA 库。公开文档提供了 [Topp、GUDHI 与 Hera 的 Python 批量距离速度表](https://proffitteoy.github.io/Topp/#python-批量距离性能)；结果仅代表表中固定环境和调用方式。

## 提供的能力

- exact Bottleneck Distance（点间使用 `L∞`）；
- exact `W1-L∞` 与 `W2-L2` Wasserstein Distance；
- 不可变的 `PreparedDiagram` 和原生 one-to-many 计算；
- 支持复用输出数组，以及 exact `bottleneck_within` 阈值判断；
- 原生计算期间释放 GIL；
- Windows x64 的 CPython 3.10–3.14 wheels；
- 运行时仅依赖 NumPy；AVX2 在运行时检测，不要求所有机器支持。

## 安装

```powershell
py -m pip install topp
```

预编译 wheel 目前仅面向 Windows x64。其他平台可尝试使用 CMake 3.24+ 和 C++20 编译器从 sdist 构建，但 Linux 和 macOS 尚未纳入 CI，不属于已验证平台。

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

### 批量比较并复用内存

```python
targets = [y, np.empty((0, 2))]
out = np.empty(len(targets), dtype=np.float64)

query = topp.prepare_diagram(x)
topp.wasserstein_distances(
    query, targets, order=2, internal_p=2, out=out
)
```

## 支持的度量

| 函数 | 语义 | 状态 |
|---|---|---|
| `bottleneck_distance` | exact Bottleneck，内部 `L∞` | 支持 |
| `wasserstein_distance(..., order=1, internal_p=np.inf)` | exact `W1-L∞` | 支持 |
| `wasserstein_distance(..., order=2, internal_p=2)` | exact `W2-L2` | 支持 |
| 其他 Wasserstein 参数 | 数学上可能合法 | `NotImplementedError` |

## 输入契约

输入必须可转换为 `(n, 2)` 的 `float64` 数组。空图、对角点、重复点和规范 essential points 合法。NaN、`birth > death`、`birth=+inf`、`death=-inf` 及其他非法无穷组合会抛出 `ValueError`，不会被静默修正。

距离定义、对角线代价、重复点和 essential points 的处理见[数学约定](docs/MATHEMATICS.md)；调用契约见 [API 文档](docs/API.md)。

## 内核维护边界

C++ 源码保留候选生成、图表示、matching、component 和 incremental pricing 等显式策略，用于回归、消融和维护。它们不会暴露到普通 Python API，也不代表默认性能承诺。1.0 默认路径、保留基准和已淘汰路线见[内核最终状态](docs/research/FINAL_STATE.md)。

## 开发

```powershell
py -m pip install -v .
py -m pytest tests/python
cmd.exe /d /c scripts\build-kernel.cmd
```

现有 `include/bottleneck/*` C++ 接口用于社区维护和内核实验，不承诺稳定 ABI。构建、测试和 benchmark 约定见 [开发指南](docs/DEVELOPMENT.md)。

## 项目导航

| 入口 | 内容 |
|---|---|
| [项目网站](https://proffitteoy.github.io/Topp/) | 适用场景、安装、API 与数学语义概览 |
| [使用说明](docs/USAGE.md) | 安装、单次与批量调用、输出数组和异常处理 |
| [API 文档](docs/API.md) | 完整公开 API 与输入契约 |
| [数学约定](docs/MATHEMATICS.md) | 距离定义、对角线、重复点与 essential points |
| [开发指南](docs/DEVELOPMENT.md) | 本地构建、测试与 benchmark |
| [贡献指南](CONTRIBUTING.md) | 正确性和性能修改的提交要求 |
| [研究记录](docs/research/README.md) | 1.0 最终状态、内核实验和差分证据 |
| [更新日志](CHANGELOG.md) | 版本能力与已知限制 |

## 引用

研究中使用 Topp 时，请引用仓库版本与发布标签。机器可读元数据见 [CITATION.cff](CITATION.cff)。

## 许可

Topp 使用 [MIT License](LICENSE)。GUDHI 仅作为测试 oracle、语义参考及历史补丁来源，不是运行时依赖；详情见 [第三方声明](THIRD_PARTY_NOTICES.md)。
