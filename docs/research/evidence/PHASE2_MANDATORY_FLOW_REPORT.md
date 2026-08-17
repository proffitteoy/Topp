# Bottleneck Phase 2 E7：mandatory partial flow 实验报告

## 结论

提案 E7 已形成 exact-safe 的 `mandatory_sparse_flow`，并在保守 near-diagonal gate 内接入默认 dispatcher。该路径只保留至少一端不能在当前阈值送往 diagonal 的 cross edge，并用 birth-coordinate x-window 枚举可能边；optional–optional cross edge 不影响可行性，因此可以精确删除。

E7 的两个变量已独立消融：

- `x_sweep_clipped` 只优化 candidate generation，单独使用并不稳定，不能无条件进入默认路径；
- `mandatory_sparse_flow` 是主要收益来源；
- x-sweep 在已经命中 E7 gate 的输入上继续带来正收益，因此只在该 gate 内与 sparse flow 组合。

这不是整个第二阶段的完成结论。E1–E5、E8–E12 仍需继续按提案逐项审计和实验。

## Exact 语义

给定阈值 `r`，点 `p` 在 `d_inf(p, diagonal) <= r` 时是 optional，否则是 mandatory。完整增广图中的 optional–optional cross match 可以被两条独立 diagonal match 替代，所以删除该类 cross edge不会改变 threshold feasibility。

对剩余边，`d_inf(p, q) <= r` 必然蕴含 `abs(birth(p) - birth(q)) <= r`。实现先以向外 `nextafter` 扩张 birth window，再计算完整的精确 `L_inf` 距离；x-window 只负责不漏候选，不能自行接受边。

候选半径同样使用可行的 all-to-diagonal upper bound。upper bound 之外的 cross distance 不可能改善当前解，因此可以不生成；窗口内仍计算精确距离并参与原有 exact threshold search。

## 实现与 router

新增公开实验策略：

- `CandidateStrategy::x_sweep_clipped`；
- `MatcherStrategy::mandatory_sparse_flow`；
- C API candidate `4`、matcher `12`；
- `x_window_candidates`、`mandatory_vertices`、`optional_pairs_pruned` 统计。

默认 distance router 只在总有限点数至少为 384 且满足以下 near-diagonal 结构之一时进入组合路径：

1. 至少一半有限点的 diagonal distance 不超过本对图最大 diagonal distance 的 5%；
2. 最大 diagonal distance 不超过合并 birth span 的 2%。

进入后使用 `x_sweep_clipped + quickselect + recompute_soa + mandatory_sparse_flow`。`bottleneck_within` 和其他仍使用 adaptive matcher 的决策路径，则在实际阈值下满足 sampled cross density `<= 3%` 或 optional–optional pair fraction `>= 75%` 时选择 sparse flow。已有 multiplicity route 和 no-cross geometric certificate 的优先级保持不变。

这些阈值来自本轮胜区/负对照，不是理论最优常数。后续 E12 必须在更完整 phase space 上重新校准。

## Benchmark 设计

正式入口为 `benchmarks/bottleneck_mandatory_bench.cpp`。它覆盖：

- 小侧 `N = 128, 256, 512, 1024, 2048`；
- cardinality ratio `1:1, 1:2, 1:4, 1:8, 1:16`，受 `--max-points` 控制；
- near-diagonal target optional ratio `0.50, 0.75, 0.90, 1.00`；
- uniform、clustered、separated 负对照；
- current adaptive、legacy full scan、x-sweep only、dense mandatory、sparse mandatory、E7 combined 六个变体。

配置逐轮随机顺序，先预热，再报告 7 轮、每轮 3 次调用的 median/p95。所有变体必须与 geometric exact reference 得到逐 bit 相同的 `double` 结果，否则 benchmark 失败。以下数据来自 Windows x64、MSVC `/O2`、AVX2 构建 `build/bn-phase2-e7-dev`；它们是单机 C++ 内核证据，不代表 Python API 或其他机器。

## 结果

### Near-diagonal 胜区

下表使用 target optional ratio `0.90`。`adaptive-before` 是接入 E7 前同一 benchmark 中的默认路径，`E7 combined` 是新组合路径；speedup 按 median 计算。

| small × large | exact optional ratio | adaptive-before median | E7 combined median | speedup |
| ---: | ---: | ---: | ---: | ---: |
| 256 × 256 | 0.898 | 1.693 ms | 0.694 ms | 2.44x |
| 512 × 512 | 0.900 | 7.713 ms | 2.441 ms | 3.16x |
| 1024 × 1024 | 0.900 | 28.120 ms | 6.376 ms | 4.41x |
| 2048 × 2048 | 0.900 | 100.924 ms | 19.734 ms | 5.11x |
| 128 × 2048 | 0.988 | 5.324 ms | 2.388 ms | 2.23x |
| 256 × 2048 | 0.977 | 11.754 ms | 3.405 ms | 3.45x |
| 1024 × 2048 | 0.934 | 56.569 ms | 12.587 ms | 4.49x |

较低 optional ratio 仍有结构性收益：512×512 的 target `0.50` 从 6.922 ms 降至 3.399 ms（2.04x），target `0.75` 从 6.805 ms 降至 3.219 ms（2.11x）。128×128 总点数低于 gate，继续走原路径。

最终 router 构建的独立烟测中，512×512、target `0.90` 的 current adaptive median/p95 为 2.376/3.462 ms，并记录 2,352 个累计 mandatory vertices 与 2,894,330 个累计 pruned optional pairs；同轮 uniform 512×512 的 adaptive median/p95 为 12.926/13.053 ms，三项 E7 route 统计均为零，确认负对照未误分派。

以 2048×2048、target `0.90` 为例，完整 cross product 为 4,194,304 对；组合路径每次 distance call 的 x-window 统计为 880,024（含 candidate generation 和各 threshold decision 的 matcher window），最终 flow capacity edge 累计仅 21,729。`optional_pairs_pruned` 是跨 17 次 threshold decision 的累计值，不能当作单图唯一 pair 数。

### 负对照

显式 E7 变体在负对照上会变慢，因此 gate 是必要条件。

| pattern / size | exact optional ratio | adaptive median | E7 combined median | 判断 |
| --- | ---: | ---: | ---: | --- |
| uniform 256×256 | 0.463 | 4.129 ms | 5.192 ms | E7 禁用 |
| uniform 512×512 | 0.387 | 12.404 ms | 15.856 ms | E7 禁用 |
| clustered 256×256 | 0.262 | 4.507 ms | 11.165 ms | E7 禁用 |
| clustered 512×512 | 0.269 | 11.573 ms | 36.097 ms | E7 禁用 |
| uniform 128×512 | 0.814 overall | 4.525 ms | 4.508 ms | 仅 overall rho 不足以判路由 |

separated 输入继续由已有 no-cross certificate 返回：128、256、512 对称规模的 adaptive median 分别约 0.003、0.006、0.012 ms。E7 不应覆盖该更便宜的路径。

### 独立消融

`x_sweep_only` 在 near-diagonal 上有时降低 candidate generation 成本，但在 256+ 规模经常与 legacy full scan 持平或更慢；它不是独立默认赢家。`mandatory_dense` 虽避免完整 dummy matching，仍扫描每次决策的 `n*m` cross pairs；在 512×512、target `0.90` 上约 13.566 ms，而 sparse mandatory 为 3.221 ms、组合为 2.441 ms。收益来自删除 optional–optional pairs和避免全 cross scan，不只是换用 flow。

## Correctness 与验证边界

专项 C++ 测试覆盖随机 near-diagonal、`1:1–1:16`、正反方向、exact radius、前一 ULP、全 optional、单侧 mandatory、birth-window 边界、death-coordinate 精确复核和统计触发。全配置回归及 GUDHI differential 的最终命令与结果应与代码变更一起保留；benchmark 的 exact equality 只证明被执行的合成矩阵，不替代独立 oracle。

当前证据仍缺少 Python API/轮子端到端计时、非 Windows 编译器的稳定性能数据、真实数据集 phase-space 覆盖，以及 threshold-shell/adversarial dense 专项 crossover。E12 之前不能把本轮手写 gate 描述为普适最优。

## 下一步

1. 在 E8/E9 中复用本轮 x-window/optional 统计，比较 geometric sparsification 与 components；
2. 在 E10 中衡量 candidate window 是否应改变 threshold search，而不是只改变 candidate materialization；
3. E12 汇总完整 `N × distribution × cardinality ratio × batch size` 矩阵，重新拟合 router，并保留本轮负对照；
4. E1–E5 的 Wasserstein 数值信息实验仍需单独结论，不能由 E7 的 transport-structure 复用替代。
