# Bottleneck Phase 2 E12：adaptive router 统一消融

## 结论

本轮将 E6–E11 的胜区和负对照放入统一 exact benchmark，并用实际 route counters 核对默认 dispatcher。最终只调整一条已有、已验证的 E7 性能门槛：`mandatory_sparse_flow` 的 distance route 从总有限点数至少 384 下调到至少 256。

该调整覆盖 128×128 以及总规模 256+ 的高 cardinality-ratio near-diagonal 输入；不改变 exact 条件，也不放宽 near-diagonal 几何 gate。统一结果支持：

- uniform/clustered 继续走 geometric refinement，adaptive 通常在显式 refinement 的约 1.00–1.09x 内；
- duplicate-heavy 继续优先走 multiplicity capacity flow，64–256 点的 adaptive regret 约 1.00–1.03x；
- separated 继续由 no-cross certificate 返回，比任何 matcher 快两个数量级以上；
- 128×128 near-diagonal 从 quickselect 切到 mandatory sparse，代表 median 由约 0.43 ms 降到约 0.29 ms；32×256 从约 0.31 ms 降到约 0.24 ms；
- E9 component、E10 incremental/blocked 与按 B 改 single-pair kernel 仍不进入 router；
- threshold-shell 上 geometric-candidate path 仍是已知 oracle winner，adaptive regret 约 2.5–3.1x，但现有 cheap point-level features 无法安全识别 candidate-distance multiplicity；本轮不为单一合成形状加入过拟合 gate。

因此 E12 的结果不是“每个样本都命中事后最快配置”，而是一个 exact-first、负对照受控的保守 router：只默认启用有稳定胜区且能由廉价结构特征可靠区分的路径。

## 当前路由顺序

distance dispatcher 的有效顺序为：

1. duplicate-heavy 且两图有限 multiset 完全相同时，identity certificate 直接返回 0；
2. duplicate-heavy 命中保守压缩率 gate 时，使用 multiplicity flow；
3. 总点数至少 128 时计算 exact-safe x-window pair count；若为 0，使用 no-cross certificate；
4. 总点数至少 256 且 near-diagonal gate 命中时，使用 x-sweep candidates + quickselect + mandatory sparse flow；
5. 其余总点数至少 128 的输入，按 x-window fraction 在 geometric refinement 与 quickselect 之间选择；
6. 更小输入保持 quickselect。

优先顺序很重要：multiplicity 和 no-cross 不应被较重的 E7 matcher 覆盖。batch size 不进入该 single-pair 顺序；batch 调用方通过持有 `PreparedDiagram` 显式摊销 preparation。

## 门槛重标证据

原 E7 专项 benchmark 在 128×128、3 repetitions × 7 rounds 上复测：

| target optional ratio | old adaptive | E7 combined | speedup |
| ---: | ---: | ---: | ---: |
| 0.50 | 0.740 ms | 0.413 ms | 1.79x |
| 0.75 | 0.504 ms | 0.267 ms | 1.89x |
| 0.90 | 0.451 ms | 0.304 ms | 1.48x |

optional ratio 1.00 时现有 no-cross shortcut约 0.002 ms，继续先于 E7 返回。相同规模的 uniform/clustered 虽然显式 E7 分别约 1.60/3.04 ms、慢于 adaptive 的约 0.84/0.99 ms，但它们不满足 point-level near-diagonal gate；门槛下调不会误分派这两个负对照。

门槛下调后的统一复测：

| workload | old route / median | new route / median | best explicit in new run |
| --- | ---: | ---: | ---: |
| near-diagonal 128×128 | quickselect / 0.434 ms | mandatory / 0.292 ms | geometric candidates / 0.279 ms |
| near-diagonal 256×256 | mandatory / 0.526 ms | mandatory / 0.438 ms | adaptive / 0.438 ms |
| near-diagonal 32×256 | quickselect / 0.306 ms | mandatory / 0.244 ms | mandatory / 0.235 ms |
| near-diagonal 256×32 | quickselect / 0.297 ms | mandatory / 0.244 ms | mandatory / 0.160 ms |

最后一档 p95 较宽，且不同方向使用不同随机 diagrams，不能把 0.160/0.244 当作稳定的方向性差异；能确认的是 route 已切换且不再承担 quickselect 的 full dummy matching 工作。

## 统一 benchmark

新增 `benchmarks/bottleneck_router_bench.cpp`，覆盖：

- uniform、clustered、near-diagonal、duplicate-heavy、separated、threshold-shell；
- 64/128/256/512 对称规模；
- uniform/clustered/near-diagonal 的 32×256、256×32、64×512、512×64；
- adaptive、refinement、quickselect、geometric candidates、multiplicity、mandatory sparse；
- duplicate fraction、near-diagonal fraction、exact cross density、实际 route、median/p95、best-measured regret 和图工作量。

变体逐轮随机执行，并逐 bit 与 geometric refinement exact reference 比较。`best` 只表示本轮显式候选配置中的最快 median，不是所有可能算法的理论 oracle。

### 正/负对照摘要

Windows x64、MSVC `/O2`、3 repetitions × 7 rounds 的 64–256 点代表结果：

| pattern / size | adaptive route | adaptive median | best measured | adaptive regret |
| --- | --- | ---: | ---: | ---: |
| uniform 64 | refinement | 0.267 ms | 0.260 ms refinement | 1.03x |
| uniform 128 | refinement | 0.746 ms | 0.714 ms refinement | 1.04x |
| uniform 256 | refinement | 3.812 ms | 3.692 ms refinement | 1.03x |
| clustered 64 | refinement | 0.310 ms | 0.310 ms refinement | 1.00x |
| clustered 128 | refinement | 1.143 ms | 1.143 ms adaptive | 1.00x |
| clustered 256 | refinement | 4.763 ms | 4.763 ms adaptive | 1.00x |
| duplicate-heavy 64 | multiplicity | 0.012 ms | 0.012 ms multiplicity | 1.03x |
| duplicate-heavy 128 | multiplicity | 0.013 ms | 0.013 ms adaptive | 1.00x |
| duplicate-heavy 256 | multiplicity | 0.012 ms | 0.012 ms multiplicity | 1.03x |
| separated 64 | no-cross | 0.001 ms | 0.001 ms adaptive | 1.00x |
| separated 128 | no-cross | 0.002 ms | 0.002 ms adaptive | 1.00x |
| separated 256 | no-cross | 0.004 ms | 0.004 ms adaptive | 1.00x |

显式 negative optimizations 的代价也被保留：uniform 256 上 multiplicity/mandatory/geometric-candidates 分别约为 adaptive 的 3.2x/1.6x/4.0x；clustered 256 上约为 3.2x/2.5x/3.3x。这是它们不能无条件开启的直接证据。

## 已知 router miss：threshold-shell

| N | adaptive refinement | geometric candidates | regret |
| ---: | ---: | ---: | ---: |
| 64 | 0.125 ms | 0.041 ms | 3.08x |
| 128 | 0.356 ms | 0.125 ms | 2.85x |
| 256 | 1.121 ms | 0.444 ms | 2.53x |

该样本 point duplicate fraction 为 0、near-diagonal fraction 为 0、x-window fraction 为 1、exact density 为 1；真正让 geometric candidates 获胜的是 cross candidate distances 高度重复且只需一次 threshold decision。计算精确 candidate duplicate ratio 本身需要接近 `O(NM)` materialization，会吞掉 router 收益；仅用 x-window fraction 或 bbox 宽度又会误伤随机 vertical-dense 输入。

因此当前把它记录为已知 miss，而不是加入未经负对照证明的规则。后续若研究 cheap sketch，必须同时加入 birth/death 窄簇但 candidate distances 连续的反例，并把 sketch 自身成本计入 regret。

## 代码与可观测性

`SolverStats` 新增六个只读计数器：

- `router_identity_shortcuts`；
- `router_no_cross_shortcuts`；
- `router_multiplicity_routes`；
- `router_mandatory_routes`；
- `router_refinement_routes`；
- `router_quickselect_routes`。

这些字段不参与决策，只让测试、benchmark 和调用方区分“结果快”与“确实命中目标路径”。专项测试覆盖 identity、no-cross、small-N quickselect、multiplicity 和 128×128/256×256 mandatory 路由。

## Router 决策

- 合入：E6 multiplicity、E7 mandatory sparse、no-cross、geometric refinement 与本轮 256 总点数 E7 门槛；
- 保留显式但不默认：component decomposition、incremental、blocked incremental、完整 mandatory flow、单独 x-sweep、按 density 单变量切 KD；
- batch：不影响 single-pair kernel，只复用 prepared diagrams/output buffer；
- Wasserstein numerical bounds：E1–E5 未完成前不进入 router。

## 剩余工作

1. E1–E5 必须单独验证 Wasserstein bound quality、完整/anytime prepass、matching warm-start 与 q-continuation；
2. threshold-shell cheap candidate-multiplicity sketch 仅在有连续-distance 负对照后继续；
3. 最终发布验证需重跑 C++ 全配置、GUDHI e=0 differential、Wasserstein 回归、两套编译和 diff hygiene。
