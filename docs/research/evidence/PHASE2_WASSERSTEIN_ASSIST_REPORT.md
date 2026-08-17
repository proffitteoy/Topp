# Bottleneck Phase 2 E1–E5：Wasserstein 数值信息与 matching screening

## 结论

本轮完成提案中“真正先计算 Wasserstein”的 E1–E5 隔离实验。结果表明，Wasserstein 数值和 matching 与 Bottleneck 高度相关，但在当前已经过 E6/E7/E12 优化的 exact router 上，完整 prepass 仍没有稳定的端到端胜区。

- E1：随 q 增大，`Lq` 普遍变紧，最优 Wasserstein matching 的 `Uq` 常在 q=4–8 达到 exact Bottleneck；候选削减率可达 90%–99% 以上；
- E2：生产 optimized `W1/L∞` 虽在 near-diagonal 很快，但 `W1 + lower-assisted Bottleneck` 的对称 32/64/128 点 median 分别只有 baseline 的 0.83x/0.95x/0.85x；不进入默认 router；
- E3：完整 greedy matching 有时给出好 upper bound，但只处理 5%–10% positive edges 在 duplicate-heavy/threshold-shell 上仍可能差 31x/1024x；生产 solver 又没有可独立验证的 anytime global dual，因此不能把 partial state 当作 certified lower bound；
- E4：Wasserstein matching 在 exact Bottleneck threshold 下的 surviving fraction 多数超过 0.9，说明结构相关性真实存在；但当前 winning routes 是 multiplicity flow、mandatory partial flow、no-cross 或 geometric refinement，不消费 augmented Kuhn seed，prepass 成本也未回本；暂不接 warm start；
- E5：相邻 q 的 matching 在 `1→2` 时经常明显变化，而 q≥4 后 upper bound 才趋于 exact；dense q-reference 在 128 点通常已比 Bottleneck 慢 6x–27x，continuation 没有工程收益。

结论不是“Wasserstein bounds 无信息”，而是“信息价值没有覆盖获取和验证信息的成本”。这两件事必须分开汇报。

## 实验实现与 exact 边界

新增 `benchmarks/bottleneck_wasserstein_assist_bench.cpp`：

- ground metric 对所有 q 统一为 `L∞`；
- 用 long-double dense Hungarian reference 求 `q=1,2,4,8,16` 的最优 positive-saving matching；
- unmatched point 继续使用 diagonal cost，得到完整 feasible Wasserstein matching；
- q=1 value 与生产 optimized `W1/L∞` 逐 case 数值交叉验证；
- 每个 matching 计算 `Wq`、`Lq=Wq/K^(1/q)`、`Uq=max edge`、候选保留数、exact-threshold survival 和相邻 q cross-edge overlap；
- 另测 global greedy matching，以及只扫描 5%/10% positive edges 后的 feasible upper bound；
- 所有 Bottleneck 结果继续由默认 exact solver计算，研究 benchmark 不修改默认 router。

这里 `K` 是该 matching 的实际非零/零 cost summand 数；它不大于 `n+m`。端到端 W1 prototype 为保持生产接口可用，只采用更保守的安全上界 `N=n+m`。

## E1：bound quality

### Near-diagonal

代表 128×128：

| q | `Lq/dB` | `Uq/dB` | candidate reduction | matching survival | previous-q overlap |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.506 | 1.437 | 98.74% | 0.962 | 1.000 |
| 2 | 0.562 | 1.019 | 99.08% | 0.996 | 0.674 |
| 4 | 0.649 | 1.019 | 99.25% | 0.996 | 0.933 |
| 8 | 0.736 | 1.000 | 99.50% | 1.000 | 0.857 |
| 16 | 0.812 | 1.000 | 99.62% | 1.000 | 0.958 |

W1 已经大幅缩小 candidate window，但 upper 仍可能比 dB 松 44%；q=8 才在该样本给 exact upper。64×64 中 q1 upper ratio 1.495，q2 为 1.009，q4 起 exact。

### Uniform / clustered

uniform 128×128 中，q1 的 lower/upper ratio 为 0.355/1.092，candidate reduction 92.68%；q2 upper 已 exact，q16 lower ratio 提升到 0.738、reduction 97.07%。

clustered 128×128 更难：q1 lower/upper ratio 为 0.401/2.100，q4 upper 仍为 1.064，q8 才 exact。q1 survival 仍有 0.945，说明“多数 W matching edges 可保留”并不等价于“最大 edge upper 已经紧”。

### Duplicate-heavy / separated / threshold-shell

- duplicate-heavy 的所有 q 都有 `Lq=Uq=dB`，候选削减约 96.6%–100%；
- separated 的 optimal matching 全部走 diagonal，所有 q 的 `Uq=dB`，candidate reduction 约 97%–99.7%；
- threshold-shell 的所有 q 同样得到 exact lower/upper，candidate reduction 98.5%–99.6%。

这三个结果从信息论角度很强，但后两类已有更便宜的 Bottleneck certificate；duplicate-heavy 已由 E6 capacity compression 直接解决。

## E2：完整 W prepass 与 lower-assisted prototype

为了不把“候选少”误当成“总时间快”，本轮给 `SolverConfig` 增加显式 `lower_bound_hint` 研究入口：

- hint 大于内建 lower bound 时，先在 hint 的前一 ULP 运行一次 exact feasibility；
- probe infeasible 才接受并裁掉更小 candidates；
- probe feasible 则安全忽略错误 hint；
- 默认值为 0，现有 router 不使用该字段；
- stats 记录 tested/accepted，测试覆盖正确 hint 与故意过大的错误 hint。

端到端 prototype 计算 optimized W1，再用 `W1/(n+m)` 作为 lower hint。Windows x64、MSVC `/O2`、1 repetition × 7 rounds 的代表 median：

| pattern / size | Bottleneck | optimized W1 | W1 + assisted B | speedup |
| --- | ---: | ---: | ---: | ---: |
| near-diagonal 32×32 | 0.028 ms | 0.007 ms | 0.034 ms | 0.83x |
| near-diagonal 64×64 | 0.303 ms | 0.026 ms | 0.318 ms | 0.95x |
| near-diagonal 128×128 | 0.387 ms | 0.066 ms | 0.455 ms | 0.85x |
| near-diagonal 32×128 | 0.223 ms | 0.023 ms | 0.224 ms | 0.99x |
| near-diagonal 128×32 | 0.289 ms | 0.028 ms | 0.252 ms | 1.15x |
| duplicate-heavy 128×128 | 0.025 ms | 0.020 ms | 0.035 ms | 0.72x |
| uniform 128×128 | 1.115 ms | 0.944 ms | 2.339 ms | 0.48x |
| clustered 128×128 | 0.703 ms | 0.696 ms | 1.418 ms | 0.50x |
| separated 128×128 | 0.004 ms | 0.034 ms | 0.008 ms | 0.52x |
| threshold-shell 128×128 | 0.386 ms | 2.562 ms | 2.696 ms | 0.14x |

near-diagonal 128×32 是唯一超过 1 的代表样本，但反方向 32×128 仅持平，而且两者不是同一对 diagrams 的转置。单一随机样本不足以建立 cardinality-direction gate；默认禁止。

注意 uniform/clustered 的 adaptive route 走 geometric refinement，它不 materialize candidate vector，因此 lower hint 不会缩短其内部 refinement；prototype 仍支付 W1。separated 在总点数至少 128 时先命中 no-cross，hint 不测试；duplicate-heavy 先命中 E6。现有优先路径已经吃掉了最有希望的结构。

## E3：greedy 与 anytime screening

完整 global greedy matching 的 upper 在 near-diagonal、duplicate-heavy、separated、threshold-shell 通常 exact；但这不能外推到“只做 5%–10% 就够”。

代表 worst cases：

| pattern | full greedy `U/dB` | 5% positive edges | 10% positive edges |
| --- | ---: | ---: | ---: |
| uniform 128 | 1.75 | 1.72 | 1.67 |
| clustered 128 | 2.99 | 2.20 | 2.20 |
| duplicate-heavy | 1.00 | 31.25 | 31.25 |
| threshold-shell | 1.00 | about 1024.4 | about 1024.4 |

partial greedy 尚未覆盖足够 cross pairs 时，大量点仍走 diagonal，upper 被最大 persistence 控制。near-diagonal/separated 的 5% bound 看似 exact，是因为 diagonal upper 本来就等于 dB，不是 partial Wasserstein 找到了额外信息。

生产 priced sparse solver 的 Dijkstra potentials 已在 Wasserstein 阶段被审计为 solver-internal potentials，不能自动当作全局 LP dual。没有独立 certificate 前，E3 不输出伪 `Dq(t)`，也不把 partial pricing state接入 Bottleneck。

## E4：matching warm-start screening

optimal W matching 在 dB 下的 surviving fraction：

- near-diagonal q1：约 0.956–1.000，q≥4 常接近 1；
- uniform q1：约 0.947–0.977；
- clustered q1：约 0.906–0.984；
- duplicate/separated/threshold-shell：1.000。

这证明 W matching 与 Bottleneck matching 有强相关性。但 warm-start 必须与当前 kernel 对齐：

- E6 使用 compressed capacity flow，原 point-level seed不直接对应压缩变量；
- E7 使用 mandatory partial circulation，已删除 optional–optional edges；
- no-cross 不执行 matching；
- geometric refinement 使用 range oracle 和自己的增广状态。

为一个高 survival fraction 重建完整 augmented graph，会丢掉当前 router 的主要收益。加上 E2 的 prepass 未回本，本轮不实现 B3/B4 production warm start。若未来某个 dense-Kuhn workload 成为真实默认瓶颈，应先在同一 matcher 内比较 empty/greedy/W-surviving seed，而不能用本轮 survival 单独宣称加速。

## E5：q-continuation

相邻 q cross-edge overlap 显示：

- near-diagonal 128 的 `q1→q2` overlap 约 0.674，之后为 0.857–0.958；
- uniform 128 的 `q1→q2` 约 0.702，之后约 0.846–0.927；
- clustered 128 的 `q1→q2` 约 0.555，q4→q8 仍只有约 0.748；
- duplicate/separated matching 稳定；threshold-shell 存在大量等价最优 matching，edge identity overlap 很低但 upper 始终 exact。

因此 q1 matching 不是普遍稳定的 q2/q4 warm start。更关键的是 128 点 dense q-reference 的时间：near-diagonal q2–q16 约为 Bottleneck 的 6.0x–6.8x，duplicate-heavy 约 20x–27x，clustered q8/q16 约 1.7x；即使 continuation 完全消除最后的 Bottleneck，也无法回本。

本轮不实现 `W1→W4→dB` continuation，也不把 q 当作默认 router feature。

## Router 决策

- 不默认运行完整 W1/Wq prepass；
- 不默认使用 partial greedy/anytime state；
- 不接 Wasserstein matching warm start；
- 不接 q-continuation；
- 保留 exact-safe `lower_bound_hint` 作为显式实验/外部 certified-bound 入口，但默认 0，错误 hint 会被前一 ULP feasibility 安全拒绝；
- 继续使用 E6–E12 中无需先计算 Wasserstein 的结构迁移赢家。

## 验证与复现

```powershell
build\bn-phase2-e7-dev\bottleneck_wasserstein_assist_bench.exe --pattern near_diagonal --max-points 128 --repetitions 1 --rounds 7
build\bn-phase2-e7-dev\bottleneck_wasserstein_assist_bench.exe --pattern duplicate_heavy --max-points 128 --repetitions 1 --rounds 7
build\bn-phase2-e7-dev\bottleneck_wasserstein_assist_bench.exe --pattern uniform --max-points 128 --repetitions 1 --rounds 7
```

最终发布级验证仍需覆盖 C++ 全配置、GUDHI exact differential、Wasserstein 回归、MSVC/CMake 两套构建和 diff hygiene。
