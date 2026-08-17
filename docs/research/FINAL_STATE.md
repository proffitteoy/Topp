# Topp 1.0 内核最终状态

- 冻结日期：2026-08-17
- 内核基线：`4cf5b4e`

本文是实验报告与发布文档之间的收口页。它说明哪些能力属于 1.0 默认实现、哪些代码继续作为维护对照，以及哪些一次性材料已从活动树删除。

## 发布契约

Topp 1.0 的稳定范围是 Python API：

- exact Bottleneck Distance，点间度量为 `L∞`；
- exact `W1-L∞` 与 `W2-L2` Wasserstein Distance；
- `PreparedDiagram`、one-to-many 原生批量调用和可复用输出数组；
- exact `bottleneck_within` 阈值判断；
- 空图、对角点、重复点、essential points 和非法输入的明确语义。

`include/bottleneck/*`、内部策略枚举、统计字段和 C++ ABI 用于维护与实验，不属于 `1.x` 稳定兼容性承诺。

## 冻结的默认路径

### Bottleneck

默认 dispatcher 保持 exact 语义，并按廉价、可验证的结构特征依次选择：

1. 完全相同有限 multiset 的 identity certificate；
2. exact duplicate-heavy 输入的 multiplicity capacity flow；
3. 无交叉边时的 no-cross certificate；
4. near-diagonal / cardinality-imbalanced 输入的 mandatory sparse flow；
5. 其余中大型输入的 geometric refinement 或 quickselect；
6. 小输入的 exact quickselect 路径。

component、search、geometry 和 router benchmark 继续保留，因为它们覆盖当前默认路径的 crossover 与负对照，不是一次性测试文件。

### Wasserstein

默认 `adaptive` 组合已验证的跨工作负载赢家：

- DSR positive-saving formulation 与共享 `PreparedDiagram`；
- blocked / AVX2 / sweep / parallel candidate generation；
- dense 或 CSR weighted graph；
- rectangular dense SAP、sparse primal-dual SAP、tiny/component decomposition；
- native batch、caller-provided output 与 workspace reuse。

arena、persistent KD、active-set pricing 和相关 matcher 仍可作为显式 exact 对照，但不进入默认 dispatcher。

## 已淘汰或归档的路线

- 先完整计算 Wasserstein 再辅助 Bottleneck 的 E1–E5 未形成稳定端到端胜区；一次性 assist benchmark 已删除，结论保留在报告中。
- 固定训练集 PGO 与专用 LTO 在部分场景收益、其他场景退化，不作为通用 wheel 构建方式；专用脚本已删除。
- 长篇项目/API/优化方案已经被实现、测试和结果报告取代，活动树不再保留；Git 历史仍可追溯。
- `benchmarks/results/`、对象文件、可执行文件、wheel、虚拟环境和构建目录均为可再生成本地产物，不进入版本控制。

## 保留的验证面

- `tests/bottleneck_core_tests.cpp` 与 `tests/wasserstein_core_tests.cpp`：C++ exact 回归和实验配置差分；
- `tests/python/`：公开 API、输入、prepared/batch/out、brute-force 和可选外部 oracle；
- `tests/test_gudhi_differential.py` 与 C test bridge：手动全配置 GUDHI `e=0` 差分；
- `benchmarks/`：默认路由、批量复用和关键负对照的可复现性能入口；
- 中英文 Sphinx 文档：公开契约与支持边界。

## 1.0 已知边界

- 发布 wheel 只覆盖 Windows x64、CPython 3.10–3.14；Linux/macOS 不属于发布 CI 平台。
- Wasserstein 不支持任意 `(order, internal_p)`。
- 公共速度表是固定机器、合成输入和调用方式下的汇总；不能外推为所有分布都同倍数领先。
- 512 点 uniform/clustered 等 dense Wasserstein 输入仍可能慢于外部实现；这不影响 exact 正确性，但属于后续性能空间。
