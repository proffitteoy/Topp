# Bottleneck Phase 2 E8/E9：geometry crossover 与 component decomposition

## 结论

本轮重新审计了提案 E8/E9。第一阶段已有的 x-sweep CSR、dynamic bitset、exact KD range oracle 与 geometric refinement 都保持 exact；新增的正式 benchmark 把 edge density、组件数与 `s_max` 从分布标签中拆成独立变量。

结论不是“KD 永远开启”或“组件越多越快”：

- 在候选距离高度重复的 threshold-shell 上，`KD matcher + quickselect` 即使在高 edge density 下也胜出，因为只需 1–2 次 threshold decision；
- 在普通 uniform/clustered 数据上，同样的 candidate materialization 会明显变慢，geometric refinement 通过窄化候选环避免了这部分成本；
- 当前完整增广图的 component decomposition 被 projection dummy complete block 限制，所有实验的 `s_max` 至少为 `0.5`，并且稳定大幅变慢；它不能进入默认 router；
- adaptive router 的 x-window 预判从每点两次二分改为单调双指针，复杂度由 `O(n log m)` 降为 `O(n + m)`，计数与路由语义不变。

E8 的 grid/WSPD-inspired speculative expansion 尚未实现；当前结果说明下一步首先需要把 candidate distance multiplicity/shell width 纳入 E10，而不是再增加一个仅按 `rho_E` 判定的空间索引。

## Exact 边界

所有 geometry 配置仍以向外 `nextafter` 扩张查询框，并对返回点重新使用精确 `L_inf` 条件。Geometric refinement 只缩小 infeasible/feasible bracket，最终从 exact `double` candidates 中选择；保留 full-candidate defensive fallback。

本轮没有采用 approximate WSPD，也没有从一次 sparse-graph infeasible 推导完整图 infeasible。新增 benchmark 中每个变体均与 geometric exact reference 做逐 bit 结果比较。

## Benchmark

### E8 density crossover

`benchmarks/bottleneck_geometry_bench.cpp` 用 block size 精确控制 exact-threshold cross density：

`rho_E = block_size / N`。

各配置共享 x-window candidate generation，比较 adaptive、geometric refinement、KD matcher、on-demand、x-sweep CSR 与 dynamic bitset。配置逐轮随机顺序，以下为 Windows x64、MSVC `/O2`、3 repetitions × 7 rounds 的 median。

受控 threshold-shell 的 1024×1024 结果：

| `rho_E` | KD matcher | refinement | x-sweep CSR | dynamic bitset |
| ---: | ---: | ---: | ---: | ---: |
| 0.00098 | 0.375 ms | 1.829 ms | 12.378 ms | 12.742 ms |
| 0.01563 | 0.356 ms | 1.507 ms | 7.092 ms | 7.485 ms |
| 0.25000 | 1.735 ms | 3.205 ms | 10.235 ms | 9.713 ms |
| 1.00000 | 11.355 ms | 17.420 ms | 29.943 ms | 19.843 ms |

这里 KD 的优势不能外推到普通数据：现有 large benchmark 的 uniform 512×512（exact density 约 0.040）中，geometric candidate path 约 68.7 ms，而 refinement 约 21.3 ms；clustered 512×512（density 约 0.019）中分别约 17.9/3.9 ms。edge density 相近而赢家相反，差别来自候选半径数量、重复度和 refinement decisions。

因此 `rho_E < rho* -> KD` 不是充分 gate。E10 应至少记录 retained unique candidates、annulus width 与 decision count，再研究 quickselect/refinement crossover。

### E9 component decomposition

`benchmarks/bottleneck_component_bench.cpp` 构造完全分离的 cross blocks，并在 exact threshold 单独统计组件。完整增广图始终额外包含一个由所有 projection dummy 构成的 complete component；对称 `N×N` 时该组件包含 `2N` 个顶点，而整个二分增广图共有 `4N` 个顶点，所以 `s_max = 0.5`。

256×256 的结果：

| cross block size | exact components | `s_max` | augmented component | adaptive | best non-component |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 257 | 0.500 | 4.555 ms | 0.171 ms | 0.173 ms mandatory sparse |
| 4 | 65 | 0.500 | 8.966 ms | 0.355 ms | 0.341 ms mandatory sparse |
| 16 | 17 | 0.500 | 12.917 ms | 0.610 ms | 0.621 ms refinement |
| 64 | 5 | 0.500 | 22.288 ms | 1.061 ms | 1.051 ms refinement |
| 256 | 2 | 0.500 | 39.803 ms | 2.284 ms | 2.273 ms refinement |

组件越少时完整 component matcher 反而继续变慢，因为它先 materialize 正反 bit graph，并在 dummy dense component 上执行大量 DFS。block size 1 时每次 distance call 累计 emitted edges 约 296 万；block size 256 时约 2329 万。

E7 mandatory sparse flow 已经避免显式 dummy block，在小组件输入上达到同等或更好的效果。现有证据不支持再把完整 augmented `component_kuhn` 接入 adaptive。若未来实现 component-aware partial flow，必须与当前全局 Dinic 做独立消融；不能用本轮“原始 cross components 很多”作为收益证明。

## 代码与统计

- `SolverStats::largest_component_vertices` 记录单次统计生命周期内见到的最大二分组件顶点数；
- component benchmark 的 `exact_components` 与 `s_max` 来自单独的 exact-radius `bottleneck_within`，不与多次 threshold search 累计值混用；
- `SolverStats::router_x_window_pairs` 记录 adaptive distance 预判的精确窗口 pair 数；专项测试将线性 sweep 与 brute-force pair count 比较；
- 线性 x-window count 使用两组已排序 births 和向外 `nextafter` 边界，lower/upper 指针均单调前进。

## Router 决策

本轮只合入无语义变化的 x-window count 复杂度优化：

- 保留现有 geometric refinement/KD/bitset 选择；
- 不启用 `component_kuhn`；
- 不按单一 `rho_E` 新增 KD gate；
- E7 near-diagonal 与 multiplicity/no-cross 的已验证优先路径保持不变。

## 剩余工作

1. E10：正式扫描 candidate multiplicity、annulus width、binary/quickselect/incremental/refinement decisions；
2. E8：只有在 E10 仍显示 range-query 构建为主要瓶颈时，再比较 exact grid 与 speculative radius expansion；
3. E9：仅在真实数据出现去 dummy 后 `k >> 1` 且全局 mandatory sparse flow 成为瓶颈时，实验 component-aware partial flow；
4. 最终 E12 必须使用随机轮序、median/p95 和真实/合成负对照统一重标 router。
