# Persistence Diagram Distance 优化内核

本仓库用于设计高性能 Persistence Diagram 距离内核。当前产品边界仍以独立 Python Bottleneck Distance 库为主；Wasserstein Distance 暂处于纯 C++ 内核优化实验阶段。GUDHI 是正确性 oracle、性能基线和可参考实现，不是本项目的产品边界。

## 当前状态

C++ exact Bottleneck 内核第一阶段已覆盖 Small/Medium-N 算法 zoo、large-N 几何 matching、approximation-assisted exact refinement 和二维 dispatcher。Python 外壳、发行包名、公共 API、wheel 和构建后端暂不处理；没有修改外部 GUDHI 源码树。

Wasserstein 已形成可正交配置的 exact 内核实验矩阵：blocked/AVX2/sweep candidates，多类 weighted graph，dense/sparse primal-dual，component/tiny-component，稳定 matching 重算，duplicate mass compression，以及 reusable batch workspace。当前证据、淘汰项和第二轮缺口见 [docs/WASSERSTEIN_KERNEL.md](docs/WASSERSTEIN_KERNEL.md)。

- GUDHI 源码：`F:\GUDHI\gudhi-devel`
- 当前源码与历史工作的核查结果：[docs/BASELINE.md](docs/BASELINE.md)
- 独立 Python 库的产品边界：[docs/PROJECT_SCOPE.md](docs/PROJECT_SCOPE.md)
- 完整优化路线：[docs/优化方向.md](docs/优化方向.md)
- 内核实验、正确性矩阵与性能结论：[docs/KERNEL_EXPERIMENTS.md](docs/KERNEL_EXPERIMENTS.md)
- Bottleneck 第一阶段报告：[docs/PHASE1_KERNEL_REPORT.md](docs/PHASE1_KERNEL_REPORT.md)
- 历史 small-N 补丁副本：[patches/0001-gudhi-3.13-small-n-neighbors.patch](patches/0001-gudhi-3.13-small-n-neighbors.patch)
- 基准测试约定：[benchmarks/README.md](benchmarks/README.md)

## 仓库结构

```text
benchmarks/  固定输入、基准脚本与结果格式说明
docs/        产品边界、基线、设计和验证记录
include/     C++ 内核接口
src/         exact solver、dispatcher 与可选 AVX2 距离内核
tests/       reference/threshold/GUDHI 差分测试
scripts/     MSVC 普通与 LTO 实验构建脚本
patches/     历史参考及实验补丁；不是最终发布形态
```

## 工作边界

1. 首个算法目标是 exact Bottleneck Distance，并以 GUDHI `e=0` 做严格差分；近似算法作为后续独立能力。
2. 公共 API 面向普通 Python/NumPy 用户，不暴露底层图结构、第三方几何类型或实验 dispatcher 细节。
3. 单对、one-to-many 和 many-to-many 分层设计；批量 API 是正式、通用的产品能力。
4. 所有优化先做结果差分，再做性能测试；内核加速、Python API 加速和批量吞吐分别报告。
5. `F:\GUDHI\gudhi-devel` 当前仅作为外部 oracle/参考源码树，不在未确认实现方案前修改它。
6. wheel、临时构建目录和原始基准结果默认不进入版本控制；可复现实验脚本、输入清单和汇总结论应进入版本控制。

## 下一步

继续补充 Bottleneck 的峰值内存/分配计数、构造型深增广压力测试和跨编译器复测。prepared/native batch 和调用方输出缓冲区已进入内核。内核验证收口前不开始 Python 包装；只有需要公平源码基线时才构建 GUDHI，并保持其外部源码树现有文件不受影响。
"# Tide" 
