# Bottleneck Phase 2 E10：exact threshold search strategy

## 结论

本轮完成提案 E10 的独立消融，比较每次重建 feasibility 的 binary、GUDHI alpha、galloping、quickselect 与 geometric refinement，以及严格按半径递增并复用 matching 的 monotone incremental 和带整块试探/回滚的 blocked incremental。

结论是：当前两种 incremental 实现都不应进入默认 router。

- uniform/clustered 的 unique candidates 多，逐 group 增量导致大量增广搜索；
- near-diagonal 虽然候选窗口很窄，但 E7 mandatory sparse flow 已更直接地消除了 optional/dummy 工作；
- threshold-shell 把未裁剪候选重复率提高到 96.9%–99.6%，且 incremental 只需加入 2 个 group，仍被显式 `N*M` 零权 projection-dummy edges 的 materialization 压倒；
- blocked incremental 将 decision 数压低，但整块试探、完整状态恢复和块内重放没有抵消 dense dummy adjacency 的成本；
- 因而“candidate 重复率高”“decision 少”或“可复用 matching”都不是充分 gate。现有 exact geometric refinement 与 E7 路由保持默认。

本轮没有把 Wasserstein bound 假接入 Bottleneck。E1–E5 尚未给出独立的 bound quality/prepass/anytime/warm-start/continuation 证据；因此 E10 当前回答的是搜索内核在不同候选形状下的成本，而不是宣称 Wasserstein 数值界已经改变默认搜索。

## Exact 与状态复用边界

`incremental_monotone` 严格按 exact `double` 半径 group 递增加边，只从较小半径状态前进，不复用未来 threshold 的 matching。

`incremental_blocked` 先整块加入并测试；若块内首次变为 feasible，则恢复试探前的 adjacency 和 matching，再按 group 单调重放。回滚后不存在来自更大半径的边或 matching state。所有 benchmark 变体均逐 bit 与 geometric exact reference 比较，未发现距离不一致。

这证明当前原型的状态语义安全，但不证明其性能值得默认启用。

## Benchmark

`benchmarks/bottleneck_search_bench.cpp` 在 uniform、clustered、near-diagonal 和 threshold-shell 上记录：

- raw/clipped/unique candidates 与未裁剪候选重复率；
- median/p95、threshold decisions 与 refinement rounds；
- incremental groups、rollbacks、adjacency checks、emitted edges 与 augment searches。

各变体按 round 随机顺序运行。以下结果来自 Windows x64、MSVC `/O2`、3 repetitions × 7 rounds；时间为单次 distance call 的 median。

### 普通与 near-diagonal 分布

| pattern | N | adaptive | quickselect rebuild | monotone incremental | blocked incremental | monotone / adaptive | blocked / adaptive |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| uniform | 64 | 0.330 ms | 0.403 ms | 2.926 ms | 0.762 ms | 8.9x | 2.3x |
| uniform | 128 | 0.934 ms | 2.350 ms | 69.643 ms | 3.978 ms | 74.5x | 4.3x |
| uniform | 256 | 7.177 ms | 10.166 ms | 1855.037 ms | 29.928 ms | 258.5x | 4.2x |
| clustered | 64 | 0.216 ms | 0.436 ms | 4.681 ms | 0.714 ms | 21.7x | 3.3x |
| clustered | 128 | 1.534 ms | 2.460 ms | 31.943 ms | 3.823 ms | 20.8x | 2.5x |
| clustered | 256 | 5.172 ms | 12.939 ms | 138.771 ms | 18.511 ms | 26.8x | 3.6x |
| near-diagonal | 64 | 0.170 ms | 0.105 ms | 2.772 ms | 0.822 ms | 16.3x | 4.8x |
| near-diagonal | 128 | 0.536 ms | 0.550 ms | 5.624 ms | 2.979 ms | 10.5x | 5.6x |
| near-diagonal | 256 | 0.640 ms | 2.579 ms | 16.771 ms | 14.017 ms | 26.2x | 21.9x |

uniform 256 中 monotone 逐个加入 2326 个 groups，累计约 17.75 亿 emitted edges 和 111953 次 augment searches；blocked 虽只做 61 次 threshold decisions，仍累计约 2436 万 emitted edges。对比 adaptive 的约 253 万 emitted edges，matching reuse 没有覆盖维护增量图的成本。

### 高重复候选正向压力样本

threshold-shell 是 incremental 最有利的受控样本：unique candidates 只随 `O(N)` 增长，而 raw candidates 为 `(N+1)^2`；exact distance 在很早的 group 即 feasible。

| N | duplicate ratio | adaptive | galloping rebuild | quickselect rebuild | monotone incremental | blocked incremental |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 64 | 96.92% | 0.085 ms | 0.065 ms | 0.100 ms | 0.566 ms | 0.788 ms |
| 128 | 98.45% | 0.228 ms | 0.286 ms | 0.308 ms | 3.116 ms | 4.680 ms |
| 256 | 99.22% | 0.761 ms | 1.018 ms | 1.325 ms | 14.448 ms | 21.994 ms |
| 512 | 99.61% | 3.287 ms | 5.581 ms | 5.809 ms | 67.505 ms | 106.808 ms |

在 N=512 时 monotone 只加入 2 个 candidate groups，却 materialize 约 4500 万 emitted edges；blocked 加入 35 个 groups、回滚 1 次，并累计约 9000 万 emitted edges。相比之下 refinement 虽做 19 轮，只发出 1024 条边，median 为 3.471 ms。

这排除了“只在高 duplicate ratio 开 incremental”的简单 gate。若未来继续研究，应先设计不显式生成 projection-dummy complete block 的 partial/incremental graph，再与 E7 mandatory sparse flow 单独消融。

## 代码与统计

- `SolverStats::incremental_groups_added` 记录增量搜索累计加入或重放的 exact radius groups；
- `SolverStats::incremental_rollbacks` 记录 blocked search 恢复试探前 adjacency/matching 的次数；
- benchmark 同时报告候选形状、搜索决策和图工作量，避免用 decision count 代替实际成本；
- benchmark 已接入 CMake、MSVC build script 与 `benchmarks/README.md`。

复现实验：

```powershell
build\bn-phase2-e7-dev\bottleneck_search_bench.exe --max-points 256 --repetitions 3 --rounds 7
build\bn-phase2-e7-dev\bottleneck_search_bench.exe --pattern threshold_shell --min-points 512 --max-points 512 --repetitions 3 --rounds 7
```

## Router 决策

- 保留 `incremental` 与 `incremental_blocked` 作为显式 exact 实验配置；
- 不为两者新增 adaptive gate；
- 保留 geometric refinement、multiplicity、no-cross 与 E7 mandatory sparse 的现有优先关系；
- 不从本轮结果推导任何 Wasserstein-bound 默认行为。

## 剩余工作

1. E11：评估 one-to-many batch 是否能复用 prepared diagrams、候选结构或安全的静态索引；
2. E12：汇总 E6–E11 的胜区和负对照，统一重标 adaptive router；
3. E1–E5：分别验证 Wasserstein bound quality、完整/anytime prepass、matching warm-start 与 `q` continuation，不能用本轮搜索原型替代；
4. E8 speculative grid/WSPD 与 E9 component-aware partial flow 仍只在证据显示对应瓶颈时继续。
