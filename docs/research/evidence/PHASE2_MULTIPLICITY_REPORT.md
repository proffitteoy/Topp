# Bottleneck 第二阶段：multiplicity capacity 实验

记录日期：2026-08-13

## 结论

本轮完成提案 E6（multiplicity）的独立实验，并将确认胜出的 exact capacity kernel 接入默认 dispatcher。该结论只覆盖 exact duplicate-heavy 输入；E1–E5 的 Wasserstein 数值界、E7 的 mandatory partial flow、E8/E9 的新消融、E10–E12 的完整搜索与 router 矩阵仍未完成。

核心结果：只合并精确相同的有限点，把 multiplicity 作为容量；候选半径和 threshold feasibility 都在 unique representatives 上计算。压缩不改变 exact Bottleneck 语义。无重复或低重复率时 capacity flow 是负优化，默认 router 不启用。

## Exact 化简

给定 threshold `r`，每个 unique point group 的 multiplicity 是容量：

- `d(x, Δ) > r` 时，该组全部流量必须走 cross edges；否则可送 diagonal；
- cross edge 仅在代表点的精确 `L∞` 距离不超过 `r` 时存在；
- source/group 与 group/sink 使用 lower/upper bounded capacity；
- bounded circulation 经 super-source/super-sink 化为 Dinic feasibility。

这是原 augmented perfect matching 的容量商图。只合并坐标逐位相同的点，不合并 near-duplicates。候选半径同样只需 unique diagonal/cross distances；multiplicity 不会产生新的距离值。

额外加入 exact multiset identity certificate：当高重复两图的 unique points 和 multiplicities 完全相同时，直接返回 0，不进入 flow 或 geometric matching。

## 实现与观测

- `MatcherStrategy::multiplicity_flow`：显式消融开关。
- `SolverStats` 新增 `multiplicity_groups`、`multiplicity_points_removed` 和 `capacity_edges`。
- `bottleneck_multiplicity_bench`：固定 raw N，扫描 duplicate ratio，并在每轮随机 variant 顺序后报告 median/p95。
- distance 与 `bottleneck_within` 使用同一个 adaptive gate。

## Crossover 与 router

环境：本机 MSVC 2022，`/O2`，单线程；每个配置 2 次内层重复、7 轮随机顺序。表中 raw N 是每个 diagram 的点数；router 规则中的 `raw` 是两图有限点总数。以下为固定合成输入的内核 median，不是跨机器 SLA。

| raw N | duplicate ratio | geometric default | capacity | capacity speedup |
|---:|---:|---:|---:|---:|
| 128 | 0.75 | 124.9 us | 163.4 us | 0.76x，负优化 |
| 128 | 0.85 | 201.0 us | 60.6 us | 3.32x |
| 512 | 0.75 | 3021.0 us | 1429.9 us | 2.11x |
| 512 | 0.80 | 3598.5 us | 929.8 us | 3.87x |
| 512 | 0.90 | 2335.0 us | 295.5 us | 7.90x |
| 2048 | 0.75 | 74.0 ms | 18.1 ms | 4.09x |
| 2048 | 0.85 | 88.0 ms | 8.0 ms | 10.97x |
| 2048 | 0.95 | 60.3 ms | 0.94 ms | 64.3x |

因此使用保守分段 gate（`raw` 是两图有限点总数，`unique` 是两图 unique groups 总数）：

- `raw < 128`：禁用；
- `128 <= raw < 512`：`unique <= ceil(raw / 6)`；
- `512 <= raw < 1024`：`unique <= ceil(raw / 5)`；
- `raw >= 1024`：`unique <= ceil(raw / 4)`。

这不是理论最优 crossover；它是当前固定机器和合成输入下的保守手写规则。边界输入仍保留显式 kernel，后续应在真实 diagrams 和 batch 维度上复测。

## 负对照

`bottleneck_large_bench` 的 512×512 uniform、near-diagonal、clustered 和 separated 未满足 duplicate gate，继续走第一阶段路径；该门控只增加 O(1) 的 group-count 判断。

原 repeated workload 是完全相同的重复多重集。若直接路由到 capacity flow，曾从 geometric 的约 0.09 ms 退化到约 0.12 ms。identity certificate 后，所有 matcher 配置都在约 0.1 us 直接返回 0。因此“重复率高”仍不是充分条件，必须先排除可直接证明的 identity。

## 正确性证据

- `build/phase2-final/bottleneck_core_tests.exe`：全部通过。
- 新增 300 组 duplicate-heavy 随机差分；每侧 0–7 个 unique groups、每组 multiplicity 1–8，含不等长与空图；binary/quickselect capacity 均与朴素 reference 完全相等。
- 每组同时验证 exact radius 为 true、前一 ULP threshold 为 false。
- 原全配置随机矩阵已包含 `multiplicity_flow`。
- GUDHI `e=0`：508 cases、5 个几何/默认配置，加 1 个 128 重复点 capacity 专例，共 2,541 次 exact comparison 通过。

本轮早期的 Wasserstein core 回归曾通过；最终验证时，工作区中并行加入的 JV reduction、KD pricing 和 parallel component 实验使 `wasserstein_core_tests` 在新 variant 3 的首个 brute-force case 失败（expected 4.683419，actual 5.423199）。这些文件不属于本轮修改，未回滚或代修。因此当前结论是 Bottleneck 目标测试通过，但全仓 Wasserstein 测试并非全绿。

## 证据边界与下一步

- 这是合成 duplicate-heavy、单机、单线程 C++ 内核证据，不是 Python/wheel 或真实 batch 吞吐结论。
- 尚未测 peak RSS；容量图规模由 stats 间接记录，不能替代进程内存轨迹。
- 当前 Dinic capacity 使用 `int`，极端超过 `INT_MAX` 的 multiplicity 会显式抛出 overflow，而不会静默截断。
- 下一项优先实验应是 E7：把 near-diagonal 比例和 cardinality imbalance 纳入 mandatory-flow crossover；仍需沿用“独立消融 -> 差分 -> negative controls -> router”的顺序。
