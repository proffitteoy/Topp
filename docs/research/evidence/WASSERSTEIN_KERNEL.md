# Wasserstein 内核实验记录

记录日期：2026-08-11

## 范围

本轮按 `../proposals/wasserstein第二阶段.md` 的方法重做 Wasserstein 优化：每个 candidate、graph representation、matcher 和 component 分支独立配置、独立差分、随机顺序 benchmark，再组合默认 adaptive kernel。

当前仍然只处理 C++ exact 内核，不处理 Python 外壳、绑定、wheel 或 approximate auction。已覆盖 `W1-L∞` 和 `W2-L2`。

## 已实现实验矩阵

### Candidate

- `dense_scalar`：A0/A1 reference，原始 birth/death all-pairs DSR；
- `dense_blocked`：32×32 rotated DSR block scan；
- `dense_avx2`：MSVC x64 AVX2，一次计算四条 saving；无 AVX2 时回退 blocked；
- `dense_parallel`：按 row 无锁分片，最多 8 个线程，每个线程只写自己的 candidate rows；
- `sweep_binary`：按 `u` 排序并对每行做 binary-search window；
- `sweep_two_pointer`：W1 merged-event active sweep；W2 回退 binary window；
- `topk_pricing_full_scan`：每行仅物化 top-k positive saving，restricted solve 后以全对扫描做 exact pricing；
- `topk_pricing_sweep`：相同 active-set 流程，但 seed/pricing 只枚举 rotated geometric window；
- `topk_pricing_kdtree`：按 rotated `u` 建 exact branch-and-bound tree，节点缓存 `u` 区间、最大 `v` 和最小列 dual；
- `topk_pricing_sweep_incremental`：保留上一轮 matching，以重建 residual network 后的 exact negative-cycle cancellation 做 D4 matching-incremental 对照；
- `topk_pricing_sweep_persistent`：跨 pricing round 持久复用 residual network、matching 和 Bellman–Ford scratch，只向网络追加新候选边；
- `adaptive`：小图 scalar；高 window density 的至少 262144 个候选对使用 parallel，较小 dense 输入用 AVX2/blocked；其余用 binary sweep。

### Weighted graph

- `dense_matrix`；
- `csr`：`row_offsets + column_indices + weights`；
- `row_vectors`；
- `fixed_degree_4/8/16/32`：每行 inline 固定槽，超过容量的整行进入 compact overflow；
- `bitmask_lazy`：只保存 positive-edge mask，访问边时从 `PreparedDiagram` 原始坐标重算 saving；
- `block_sparse_16/32`：只实例化含正边的 dense tile；
- `adaptive`：根据 `n,m,E,ρ` 选择 dense 或 CSR；矩形 SAP 完成后，中等以上密度的非对称图可以安全选择 dense。

### Matcher 和 component

- `dense_hungarian`：完整方阵 baseline；
- `dense_sap`：独立的矩形 shortest augmenting path；删除无正边的行列、以较短侧为 augmentation 侧、不补成方阵，并跨 augmentation 复用 scan buffers；
- `dense_sap_row_reduction`：独立 exact 实验；先做 row dual reduction，并把无冲突的 row-minimum 零 reduced-cost 边作为初始 partial matching，再只增广冲突行；
- `dense_sap_jv_reduction`：column reduction、reduction transfer 和两轮 augmenting-row reduction 的直接矩形迁移；返回前以 matching LP dual certificate 验证，失败则回退 clean SAP；
- `sparse_sap`：CSR/row-list residual graph 上的 exact primal-dual Dijkstra；
- `priced restricted SAP`：`k=2/4/8/16/32`，component-aware restricted solve，返回可验证的 matching LP dual；
- `component + dense`；
- `component + sparse`；
- `tiny component`：`1×k/k×1` 闭式，最多 8 个顶点时 stack-local exhaustive；
- `row_max/global_descending` greedy certificate：只有达到 row/column upper bound 时直接返回，否则继续 exact solver；
- `adaptive`：同时查看 `n,m,E,ρ`、非对称性和 component 规模。

原先的 SPFA residual search 已淘汰。它在 GCC extended `long double` 下通过，但 MSVC 的 `long double==double` 会在某个随机 residual graph 上出现浮点负零循环。替换为可行势函数 + non-negative reduced cost 的 primal-dual Dijkstra 后，同一回归从超过 5 分钟卡住变为约 0.6 秒通过。

## 数学不变量

有限部分始终解同一个 positive-saving matching：

```text
D = Σ ai + Σ bj
sij = ai + bj - cij
Wq,p(X,Y)^q = D - max_M Σ sij
```

只保留 `sij>0`。finite、empty、diagonal、duplicate 和 essential points 继续沿用 `PreparedDiagram` 的拆分；essential 类型或数量不一致返回 infinity。

所有 candidate/graph/matcher/component 组合都必须与 `dense_scalar + dense_matrix + dense_hungarian` reference 对齐。

## 正确性证据

测试文件：`tests/wasserstein_core_tests.cpp`。

- 手工 empty、identity、shift、diagonal、zero-saving boundary、essential 和 disconnected component；
- 每个 metric 150 组 `n,m<=5` exhaustive partial-matching oracle；
- 每个 metric 500 组随机 diagram，52 条实验配置逐项差分；
- duplicate、equal cost/tie、near-diagonal、extreme coordinates、`n≠m`、finite+essential 混合；
- candidate/graph/matcher/component/adaptive 交换 diagram 后的对称性；
- GCC 和 MSVC/AVX2 两套构建。

当前结果：

```text
GCC wasserstein_core_tests: all checks passed
MSVC wasserstein_core_tests: all checks passed
```

普通随机、exhaustive 和性能矩阵跨路径门限均为 `2e-12`。所有 solver 现在返回 matching，最终值从原始 cross/diagonal cost 稳定重算，不再使用易消去的 `D-S*` 浮点减法；专门的 128 点 near-identical dense W2 回归使用 `2e-14` 门限。

## Benchmark 方法

```powershell
build\manual\wasserstein_core_bench.exe `
  --repetitions 5 --rounds 7 `
  --min-size 128 --max-size 128 `
  --pattern near_diagonal
```

每轮随机打乱 56 个配置的执行顺序，报告 median/p95，并拆分：

```text
T_prepare + T_candidate + T_graph + T_component + T_solver + T_pricing
```

同时输出 candidate/edge density、average/max degree、component count、largest component、augmentations、pricing rounds、priced/violated/materialized edges、greedy certificates、fallback count、累计与 peak graph bytes。

覆盖：uniform、near-diagonal、clustered、separated、duplicate-heavy、imbalanced、adversarial-dense、adversarial-sparse；规模清单已经包含 8–8192。超过 512 的 dense/普通分布只运行显式选择的 priced/adaptive 路径，不再误跑 O(N³) dense baseline。

## 第一轮 MSVC 结果

128×128，5 repetitions × 7 randomized rounds：

| 分布 | metric | winner median | adaptive median | adaptive/winner | edge density |
| --- | --- | ---: | ---: | ---: | ---: |
| uniform | W1-L∞ | AVX2 dense 864.6 μs | 909.1 μs | 1.05× | 21.1% |
| uniform | W2-L2 | sweep dense 1155.9 μs | 1170.9 μs | 1.01× | 20.5% |
| near-diagonal | W1-L∞ | tiny component 28.7 μs | 30.4 μs | 1.06× | 0.4% |
| near-diagonal | W2-L2 | tiny component 27.2 μs | 27.5 μs | 1.01× | 0.4% |
| clustered | W1-L∞ | dense Hungarian 1257.4 μs | 1309.1 μs | 1.04× | 99.9% |
| clustered | W2-L2 | component dense 4129.8 μs | 4911.9 μs | 1.19× | 100% |
| separated | W1-L∞ | sweep CSR 2.3 μs | 2.8 μs | 1.19× | 0% |
| separated | W2-L2 | sweep CSR 2.6 μs | 2.8 μs | 1.09× | 0% |
| duplicate-heavy | W1-L∞ | adaptive 827.7 μs | 827.7 μs | 1.00× | 50% |
| duplicate-heavy | W2-L2 | component dense 784.6 μs | 848.1 μs | 1.08× | 50% |
| adversarial-dense | W1-L∞ | dense Hungarian 1415.1 μs | 1437.6 μs | 1.02× | 100% |
| adversarial-dense | W2-L2 | component dense 5088.0 μs | 5386.9 μs | 1.06× | 100% |
| adversarial-sparse | W1-L∞ | adaptive 40.9 μs | 40.9 μs | 1.00× | 0.8% |
| adversarial-sparse | W2-L2 | greedy tiny sparse 38.6 μs | 40.4 μs | 1.05× | 0.8% |

32×128 imbalanced：adaptive 分别为 winner 的 1.08×（W1）和 1.07×（W2）。

8192×8192 large sparse：

| 分布 | metric | winner | adaptive/winner | 主要耗时 |
| --- | --- | ---: | ---: | --- |
| separated | W1-L∞ | 317 μs | 1.00× | candidate 206 μs |
| separated | W2-L2 | 315 μs | 1.01× | candidate 197 μs |
| adversarial-sparse | W1-L∞ | 2878 μs | 1.01× | component discovery 1926 μs |
| adversarial-sparse | W2-L2 | 2793 μs | 1.02× | component discovery 1913 μs |

`1×k/k×1` 闭式分支把 8192 个 `1×1` component 的耗时从约 48–51 ms 降到约 2.8–2.9 ms。

## 第二轮：矩形 dense SAP

原先的 `dense_hungarian` 和 `dense_sap` 实际共享同一套方阵算法，主要区别只有 active compaction，不能算独立 dense matcher。现在 `dense_sap` 已替换为独立矩形 shortest-augmenting-path：若 `n>m` 就转置局部访问，以 `min(n,m)` 次 augmentation 扫描长边，不再把 `128×512` 补成 `512×512`。

新增测试覆盖 `1×64`、`64×1`、`7×64`、`64×7`、`31×128`、`128×31`，含 duplicate，并与完整方阵 baseline 双 metric 差分；MSVC 与 GCC 全配置回归继续通过。

MSVC，`128×512` imbalanced，5 repetitions × 9 randomized rounds：

| metric | full-square Hungarian | rectangular SAP | adaptive before | adaptive after | after / SAP |
| --- | ---: | ---: | ---: | ---: | ---: |
| W1-L∞ | 127.82 ms | 1.18 ms | 47.51 ms | 0.88 ms | 0.74× |
| W2-L2 | 155.43 ms | 3.39 ms | 57.30 ms | 3.51 ms | 1.03× |

`adaptive after` 与显式 SAP 的小幅倒挂来自 candidate 分支和轮间频率波动，不代表 matcher 能低于自身；重要结论是数量级一致。旧 dispatcher 曾在四倍非对称且 density `<0.30` 时强制 sparse，这条规则已被实验推翻：现在 `ρ>=0.15` 选择 dense matrix + rectangular SAP，低密度非对称图仍走 CSR + sparse SAP。

方形 `128×128` 上，新 SAP 的收益依赖权重结构：本轮 W2 uniform/clustered 相对完整方阵 baseline 有收益，W1 uniform/adversarial-dense 未稳定占优。默认仍以结构分派，不删除 Hungarian oracle；完整 Jonker-Volgenant reduction transfer 等初始化技巧仍作为后续独立实验，不把矩形 SAP 冒充完整 JV。

## Native batch / caller-buffer 实验

已加入 prepared one-to-many 的 caller-buffer 和 allocated-vector overload。benchmark 改为预热、逐轮随机模式顺序、median/p95；旧的固定顺序单次测量已废弃。

MSVC，64 点 query 对 128 个 target，5 repetitions × 21 randomized rounds；workspace 已扩展为复用 candidate row capacities、graph storage 和 dense SAP scratch：

| mode | median per distance |
| --- | ---: |
| one-shot | 262.4 μs |
| prepared loop | 249.6 μs |
| native caller-buffer | 224.9 μs |
| native workspace | 221.7 μs |
| native allocated-vector | 224.1 μs |

native batch 现在自动创建一个本地 workspace，在同批 128 个 target 间复用；相对手写 prepared loop 的 median 改善约 9.9%。显式 `WassersteinWorkspace` 还能跨多次 batch 保留 capacity，本轮再改善约 1.4%。allocated-vector 复用同一 native batch 内核，因此分配输出 vector 没有抹掉主体收益。当前 workspace 尚未复用 sparse residual network，但 dense repeated-query 路径已有稳定端到端收益。

## Parallel candidate 实验

`dense_parallel` 按 candidate row 分片，线程之间不共享写入位置，只在结束时归并 candidate/positive edge 计数；少于 64 行时回退 scalar。显式实验配置继续固定为 dense matrix + rectangular SAP，并关闭 duplicate compression，避免把 graph/matcher 或 mass solver 收益计入 candidate 对照。

MSVC，`512×512`，5 repetitions × 11 randomized rounds；通过 `--experiments avx2_dense,parallel_dense,adaptive` 在同一进程内配对：

| 分布 | metric | AVX2 total | parallel total | total ratio | AVX2 candidate | parallel candidate | candidate ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| uniform | W1-L∞ | 23.82 ms | 24.13 ms | 1.013× | 1.300 ms | 1.299 ms | 1.000× |
| uniform | W2-L2 | 51.83 ms | 49.47 ms | 0.955× | 1.284 ms | 1.300 ms | 1.012× |
| clustered | W1-L∞ | 35.96 ms | 33.88 ms | 0.942× | 3.803 ms | 2.068 ms | 0.544× |
| clustered | W2-L2 | 240.90 ms | 240.23 ms | 0.997× | 3.283 ms | 1.991 ms | 0.607× |
| adversarial-dense | W1-L∞ | 42.63 ms | 40.18 ms | 0.943× | 3.836 ms | 2.117 ms | 0.552× |
| adversarial-dense | W2-L2 | 264.69 ms | 262.39 ms | 0.991× | 3.543 ms | 2.015 ms | 0.569× |

uniform 的 sampled window density 约 `0.31–0.32`，并不属于 dense candidate 分支；其端到端差异主要来自占绝对多数的 solver 波动，parallel candidate 本身没有收益。clustered/adversarial-dense 的 density 约为 `1.0`，candidate 稳定缩短约 39%–46%，W1 端到端改善约 5.7%–5.8%，W2 改善约 0.3%–0.9%。因此 adaptive 只在 `n*m>=262144 && sampled_density>=0.75` 时启用 parallel；128/256 和中低 density 输入继续保留原路径。

## Dense SAP row-reduction 实验

这不是完整 Jonker–Volgenant。当前实验只增加 exact row reduction：每个短侧 row 先计算最小 cost 作为 dual potential，把列不冲突的 row minimum 直接加入 partial matching；由于这些边的 reduced cost 为零且 dual 可行，剩余冲突 row 继续用原 SAP 增广即可保持全局 exact。完整 column reduction transfer、augmenting-row reduction 和 JV 的 free-row 两轮策略仍未实现。

MSVC，5 repetitions × 7–9 randomized rounds，`dense_sap` 与 `row_reduced_dense_sap` 同进程配对：

| 场景 | metric | baseline | row reduction | ratio | baseline augmentations | reduced augmentations |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| uniform 128 | W1-L∞ | 0.962 ms | 1.013 ms | 1.053× | 128 | 65 |
| clustered 128 | W1-L∞ | 1.412 ms | 1.292 ms | 0.915× | 128 | 63 |
| adversarial-dense 128 | W1-L∞ | 1.559 ms | 1.470 ms | 0.943× | 128 | 67 |
| uniform 512 | W1-L∞ | 25.42 ms | 24.41 ms | 0.961× | 511 | 276 |
| clustered 512 | W1-L∞ | 36.06 ms | 34.57 ms | 0.959× | 512 | 276 |
| adversarial-dense 512 | W1-L∞ | 42.55 ms | 41.99 ms | 0.987× | 512 | 272 |
| imbalanced 128×512 | W1-L∞ | 0.983 ms | 0.890 ms | 0.906× | 128 | 31 |
| uniform 512 | W2-L2 | 52.84 ms | 53.30 ms | 1.009× | 511 | 469 |
| clustered 512 | W2-L2 | 242.45 ms | 241.81 ms | 0.997× | 512 | 503 |
| adversarial-dense 512 | W2-L2 | 268.88 ms | 270.23 ms | 1.005× | 512 | 509 |

W1 的 best-column 冲突较少，512 档增广约减少 46%，四类 512/非对称输入端到端均有 1.3%–9.4% 收益；128 uniform 的额外 reduction scan 反而使总时间退化 5.3%。W2 的冲突率高，512 档只减少 1%–8% 增广，收益不稳定。因此 adaptive 仅在 `metric=W1-L∞ && max(n,m)>=512` 且原本已选择 dense matcher 时启用 row reduction；W2 和较小 W1 保留原矩形 SAP。

### 完整 JV 初始化直接迁移

`dense_sap_jv_reduction` 增加了 column reduction、reduction transfer 和两轮 augmenting-row reduction，再把剩余 free rows 交给现有矩形 SAP。由于经典 JV 初始化依赖方阵/虚拟指派约定，直接迁移到“短侧全指派、长侧可空”的 optional rectangular matching 后不能仅凭内部势函数假定 exact；因此每次先从返回 matching 构造并验证 LP dual certificate，失败时计入 `jv_fallbacks` 并回退 clean SAP。

MSVC 同进程随机轮序摘要：

| 场景 | metric | clean SAP | row reduction | JV direct | JV/clean | clean aug. | JV aug. | fallback |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| uniform 128 | W1 | 1.035 ms | 0.983 ms | 2.574 ms | 2.49× | 2639 | 3563 | 21/21 |
| uniform 128 | W2 | 1.453 ms | 1.402 ms | 3.643 ms | 2.51× | 2681 | 2730 | 7/21 |
| clustered 128 | W1 | 1.378 ms | 1.297 ms | 3.429 ms | 2.49× | 5760 | 1539 | 0/45 |
| clustered 128 | W2 | 4.459 ms | 4.810 ms | 11.366 ms | 2.55× | 5760 | 7749 | 27/45 |
| clustered 512 | W1 | 32.35 ms | 31.29 ms | 117.12 ms | 3.62× | 10752 | 3213 | 0/21 |
| clustered 512 | W2 | 241.62 ms | 237.26 ms | 414.60 ms | 1.72× | 10752 | 8848 | 0/21 |
| imbalanced 128×512 | W1 | 1.053 ms | 1.057 ms | 3.196 ms | 3.04× | 1920 | 330 | 15/15 |
| imbalanced 128×512 | W2 | 2.133 ms | 2.362 ms | 9.303 ms | 4.36× | 1920 | 1340 | 15/15 |

clustered W1 的 JV 初始化能稳定通过 exact certificate，并把增广次数减少约 70%；但多次全矩阵 reduction scan 和 dual certificate 构造仍远超节省的 SAP 工作。uniform/非对称图还暴露了经典方阵 reduction transfer 不能无条件作为矩形 optional matching 初态的问题。结论：完整 JV 初始化的直接矩形迁移实验完成，保留 exact-guarded 显式负面对照，不进入 adaptive；已冻结的基础 row reduction 仍是 W1 large dense 的赢家。

## Large-N top-k + exact pricing 实验

新增两条独立 exact active-set 路线，均测试 `k=2/4/8/16/32`：

```text
top-k seed graph
→ component-aware restricted sparse solve
→ 从 matching 构造 α/β LP dual
→ pricing: 检查 s(i,j) - α(i) - β(j) > 0
→ 每行加入 top-k violated edges
→ 无 violation 时以 dual certificate 结束
```

`topk_pricing_full_scan` 是 E1 baseline，每轮检查所有尚未物化的 pair；`topk_pricing_sweep` 是 E2，利用 W1/W2 rotated window 只枚举可能为正的 saving。restricted graph 先分 component，`<=8` 顶点走 exhaustive，其余走 sparse primal-dual。dual 会同时验证非负性、所有已物化边约束以及 primal/dual objective equality；证书无效时保守加入全部遗漏正边，不能以近似答案提前结束。

MSVC 同进程随机轮序摘要：

| 场景 | metric | full-scan priced | sweep priced | adaptive | priced/adaptive | peak materialized | full positive |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| near-diagonal 512 | W1 | 2.873 ms (`k=8`) | 1.348 ms (`k=8`) | 0.544 ms | 2.48× | 1127 | 约 1132 |
| near-diagonal 512 | W2 | 2.700 ms (`k=8`) | 1.108 ms (`k=8`) | 0.476 ms | 2.33× | 1052 | 约 1054 |
| adversarial-sparse 512 | W1 | 1.919 ms | 0.603 ms | 0.191 ms | 3.16× | 512 | 512 |
| adversarial-sparse 512 | W2 | 2.094 ms | 0.629 ms | 0.202 ms | 3.11× | 512 | 512 |
| uniform 512 | W1 | 500.9 ms (`k=32`) | 489.2 ms (`k=32`) | 18.47 ms | 26.5× | 15796 | 约 55800 |
| near-diagonal 1024 | W1 | — | 15.25 ms (`k=32`) | 13.85 ms | 1.10× | 4661 | 约 4661 |
| near-diagonal 1024 | W2 | — | 9.93 ms (`k=32`) | 8.17 ms | 1.22× | 4343 | 约 4343 |
| near-diagonal 2048 | W1 | — | 1.408 s (`k=32`) | 1.414 s | 1.00× | 18239 | 18239 |

full scan 的 512 near-diagonal pricing 从约 549 万次 pair 检查降到 sweep 的约 1.1 万次；说明 rotated pricing oracle 有效。uniform 512 的 peak graph 从 dense 2 MiB 降到约 189 KiB、物化边降约 72%，但 repeated sparse solve 令总时间慢 26.5×。1024/2048 near-diagonal 的 `k=32` 已覆盖全部正边，没有产生 active-set 边数收益；接近持平来自相同 solver 主体，不能算 pricing 胜出。

### D4 matching-incremental baseline

`topk_pricing_sweep_incremental` 第一轮仍使用 component-aware restricted solve；后续 pricing 加边后保留上一轮 matching，重建 residual network，以该 matching 初始化 flow，并通过 `sink→source` 零成本边上的 Bellman–Ford 式 negative-cycle cancellation 恢复 exact optimum。每轮结束仍重建并验证 matching LP dual，因此该分支复用了 matching，但尚未持久复用 residual network、dual 和 workspace。

当前 39 条配置在 MinGW GCC 13.2 和 MSVC/AVX2 下均通过 exhaustive/random/tie/duplicate/非对称与 near-identical 回归。MSVC 同进程随机轮序结果：

| 场景 | metric | sweep rebuild | matching incremental | incremental/rebuild | adaptive | pricing violations |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| uniform 128, `k=8` | W1 | 24.64 ms | 194.50 ms | 7.89× | 1.38 ms | 4725 |
| uniform 128, `k=8` | W2 | 22.88 ms | 853.60 ms | 37.31× | 2.17 ms | 23639 |
| uniform 128, `k=32` | W1 | 16.62 ms | 23.58 ms | 1.42× | 1.38 ms | 42 |
| uniform 128, `k=32` | W2 | 12.77 ms | 10.27 ms | 0.80× | 2.17 ms | 0 |

512 uniform 的 `k=32`、两种 metric、1 repetition × 5 randomized rounds 配对在 10 分钟上限内未完成；该组没有用超时前的部分轮次生成 median。cycle cancellation 只有在几乎没有 violation 时可能局部微赢，一旦 active set 需要实质扩张，Bellman–Ford 式修复就会远慢于重解；即使 128 W2 `k=32` 局部快约 20%，仍比 adaptive 慢约 4.7×。因此这条 exact matching-incremental baseline 保留为显式负面对照，不进入 adaptive，也不再扩大规模矩阵。

结论：D1 top-k、D2 restricted solve、D3 full/sweep exact pricing，以及一条 D4 matching-incremental baseline 已形成可复现实验。当前 active-set 路线没有击败有代表性的现有 adaptive；下一步 D4 只继续研究真正持久的 residual/dual/workspace 增量求解。KD-tree pricing 保留到增量 solver 的主体成本被解决之后。

### D4 persistent residual baseline

`topk_pricing_sweep_persistent` 在首轮 restricted solve 后构造一次 residual network，并保留 source/row/column/sink 边、当前 flow matching、每行已知边索引、Bellman–Ford distance/predecessor scratch；pricing 加边后只追加对应 residual edge，不再重建 residual network。候选 CSR 目前仍为 matching LP dual certificate 重建，因此这条分支完成了 residual/matching/scratch 持久化，但还没有消除 graph-layer certificate view。

MSVC 同进程随机轮序摘要：

| 场景 | metric/k | sweep rebuild | matching incremental | persistent residual | persistent / matching incremental | adaptive |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| uniform 64 | W1, `k=8` | 2.669 ms | 3.144 ms | 3.141 ms | 1.00× | 0.211 ms |
| uniform 64 | W2, `k=8` | 3.699 ms | 11.014 ms | 10.944 ms | 0.99× | 0.241 ms |
| uniform 64 | W1, `k=32` | 1.827 ms | 1.673 ms | 1.819 ms | 1.09× | 0.211 ms |
| uniform 128 | W1, `k=8` | 15.13 ms | 112.35 ms | 113.90 ms | 1.01× | 0.974 ms |
| uniform 128 | W2, `k=8` | 21.32 ms | 616.11 ms | 638.72 ms | 1.04× | 1.251 ms |
| uniform 128 | W2, `k=32` | 10.02 ms | 9.64 ms | 9.72 ms | 1.01× | 1.251 ms |

residual 容器与 scratch 的持久复用没有产生可测收益，差异都在约 ±9% 的轮间范围内；小 `k` 的灾难性退化完全保留。热点因此被进一步定位为每轮全 residual graph 的 Bellman–Ford negative-cycle detection，而不是 network allocation、matching materialization 或 scratch 初始化。结论：persistent exact baseline 完成并保留为显式负面对照，不进入 adaptive；继续优化这条路线必须更换增量最短路/负环算法，而不是继续做容器复用。

### E3 KD-tree pricing

`topk_pricing_kdtree` 与 sweep 使用相同 top-k seed 和 restricted solver，只替换 pricing oracle。树按 second diagram 的 rotated `u` 排序；每个节点缓存 `min/max u`、`max v` 与当前 `min β`，由此给出 W1 的 `2 min(v_i,max_v)-min|du|` 或 W2 的 `4 v_i max_v-2 min|du|²` saving 上界。若上界不足以超过 `α_i+min β` 则整棵子树 exact 剪枝；节点摘要的 binary64 舍入使用向外容差，dual 无效时回退 sweep 全枚举。

MSVC 同进程随机轮序摘要：

| 场景 | metric/k | sweep priced edges | tree priced edges | sweep pricing | tree pricing | sweep total | tree total |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| uniform 128 | W1, `k=8` | 226863 | 4718 | 0.577 ms | 0.735 ms | 18.64 ms | 20.47 ms |
| uniform 128 | W2, `k=8` | 242473 | 23639 | 0.760 ms | 1.399 ms | 28.84 ms | 28.76 ms |
| near-diagonal 128 | W1, `k=32` | 700 | 0 | 0.013 ms | 0.185 ms | 0.145 ms | 0.276 ms |
| uniform 512 | W1, `k=32` | 975330 | 1280 | 9.002 ms | 6.511 ms | 588.8 ms | 623.5 ms |
| uniform 512 | W2, `k=32` | 1093470 | 75850 | 17.590 ms | 17.049 ms | 1088.8 ms | 1051.1 ms |

树在 128 uniform 将叶子 exact pricing 访问减少 90%–98%，在 512 W1 减少 99.9%，说明 E3 上界和剪枝有效；但 128 的建树/遍历开销高于简单连续 sweep，near-diagonal 即使零叶子访问仍慢一个数量级。512 的 pricing 局部最多缩短约 28%，总耗时由 repeated sparse solve 主导，W1 反而慢约 5.9%，W2 约快 3.5% 但仍比 adaptive 慢约 17×。结论：E3 exact 实验完成并保留显式配置，不进入 adaptive；在 active-set solver 主体未解决前不实现 E4 adaptive pricing。

## B3 fixed-small-degree adjacency 实验

新增 `fixed_degree_4/8/16/32` 四档 exact representation。每行 degree 不超过容量时直接放入连续 inline 槽；超出容量时整行进入共享 compact overflow，因此不存在截断或近似。四档均加入全配置随机差分，MinGW GCC 13.2 与 MSVC/AVX2 通过。

MSVC 随机轮序结果摘要：

| 场景 | metric | CSR | fixed 最好值 | fixed/CSR | CSR bytes | fixed bytes |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| near-diagonal 128 | W1 | 560.8 μs | 353.4 μs (`K=16`) | 0.63× | 1.8 KiB | 25.1 KiB |
| near-diagonal 128 | W2 | 503.6 μs | 498.1 μs (`K=8`) | 0.99× | 1.7 KiB | 13.1 KiB |
| near-diagonal 512 | W1 | 34.47 ms | 35.88 ms (`K=32`) | 1.04× | 16.9 KiB | 196.5 KiB |
| near-diagonal 512 | W2 | 60.83 ms | 58.08 ms (`K=32`) | 0.95× | 16.3 KiB | 196.5 KiB |
| adversarial-sparse 512 | W1 | 28.59 ms | 26.55 ms (`K=8`) | 0.93× | 10.0 KiB | 52.5 KiB |
| adversarial-sparse 512 | W2 | 41.75 ms | 35.11 ms (`K=4`) | 0.84× | 10.0 KiB | 28.5 KiB |

这些是 `components=none + sparse SAP` 的表示层对照；在同一数据上，component/tiny adaptive 只需约 `0.18–0.84 ms`。因此 fixed-degree 的局部 scan 收益没有落到真实端到端 winner 上，而且容量越大，零填充、初始化和内存放大越明显。结论：四档保留为显式可复现实验，不进入 adaptive；后续 bitmask/lazy-weight 应直接与 CSR + component winner 比较，而不是只击败无 component 的 sparse solver。

## B4 bitmask + lazy-weight 实验

`bitmask_lazy` 只存每行 positive-saving bit mask；遍历命中位时，从两个 `PreparedDiagram` 的原始 finite point 与 half-persistence 重算 saving。mask 只是候选存在性结构，权重没有缓存。直接 sparse solve 和 `tiny_sparse` component 两条配置均加入 exact 差分；MinGW GCC 13.2 与 MSVC/AVX2 通过。

MSVC，5 repetitions × 7 randomized rounds，端到端 `tiny_sparse` 对照：

| 场景 | metric | N | CSR tiny | bitmask lazy tiny | bitmask/CSR | CSR bytes | mask bytes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| near-diagonal | W1 | 64 | 21.6 μs | 20.5 μs | 0.95× | 0.71 KiB | 0.50 KiB |
| near-diagonal | W1 | 128 | 33.8 μs | 38.5 μs | 1.14× | 1.80 KiB | 2.00 KiB |
| near-diagonal | W1 | 256 | 160.9 μs | 171.1 μs | 1.06× | 5.28 KiB | 8.00 KiB |
| near-diagonal | W2 | 128 | 38.3 μs | 66.7 μs | 1.74× | 1.73 KiB | 2.00 KiB |
| adversarial-sparse | W1 | 128 | 38.5 μs | 53.0 μs | 1.38× | 2.51 KiB | 2.00 KiB |
| adversarial-sparse | W2 | 256 | 89.7 μs | 94.7 μs | 1.06× | 5.01 KiB | 8.00 KiB |

64 点的少量胜出没有延续到 128/256；component discovery 会多次遍历边，lazy 原始权重重算逐渐超过少一次 weight load 的收益。mask 内存又按 `nm/8` 固定增长，在当前低 degree 图上从 128/256 点开始高于 CSR。结论：保留为 `N≤256` 的显式 exact 实验，不进入 adaptive。后续如果 matcher 能直接消费 mask、避免 component materialization，可重新评估；当前不能用单次 graph-build 降低冒充端到端收益。

## B5 block-sparse adjacency 实验

新增 `block_sparse_16/32`。只为含正边的 tile 分配 dense weight block，按 block-row 保存 tile columns；最后一个不完整 tile 仍按实际列数遍历。两档均与 sparse SAP 和 `tiny_sparse` component 正交组合，并加入 exact 差分；MinGW GCC 13.2 与 MSVC/AVX2 通过。

MSVC，3 repetitions × 5 randomized rounds，端到端 `tiny_sparse` 摘要：

| 场景 | metric | N | CSR tiny | block16 tiny | block32 tiny | CSR bytes | block bytes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| near-diagonal | W1 | 128 | 84.5 μs | 246.0 μs | 186.1 μs | 1.8 KiB | 80–125 KiB |
| near-diagonal | W2 | 128 | 79.5 μs | 209.9 μs | 155.6 μs | 1.7 KiB | 74–123 KiB |
| near-diagonal | W1 | 256 | 252.1 μs | 538.7 μs | 646.6 μs | 5.1 KiB | 327–504 KiB |
| adversarial-sparse | W1 | 128 | 48.4 μs | 91.5 μs | 128.7 μs | 2.5 KiB | 16–32 KiB |
| adversarial-sparse | W2 | 256 | 169.4 μs | 208.2 μs | 352.3 μs | 5.0 KiB | 32–64 KiB |

在 clustered 128 W1 的无 component sparse matcher 中，block32 曾把 CSR 的 `35.24 ms` 降到 `23.03 ms`，同时把 CSR 197 KiB 降到接近 dense matrix 的 128 KiB；但同一输入的 adaptive dense kernel 只需 `1.93 ms`。W2 clustered 没有复现该收益。局部 tile 压缩只在本来就应选择 dense matcher 的高密度图有一次局部胜出，在真正低 degree 图上则为少量边分配整个 tile，内存和清零成本显著放大。结论：两档保留为显式 exact 实验，不进入 adaptive。

## Matching 返回与稳定原始 cost 重算

`dense_hungarian`、rectangular SAP、sparse primal-dual、greedy certificate、tiny exhaustive 和 component closed-form 现在统一返回 `row_to_column` matching。component 局部下标会映射回全局行列；zero-saving dummy assignment 不进入最终 matching。

最终 powered distance 直接累加：

```text
essential cost
+ matched pair 的原始 cross cost
+ unmatched X/Y 点的原始 diagonal cost
```

因此 solver 仍可在 positive-saving 图上优化，但最终数值不再执行两个约为同一大数的 `D-S*`。新增 128 点 dense W2 near-identical 回归，真实距离约 `3.4e-8`；52 个 candidate/graph/matcher/component/duplicate/pricing 配置在 MinGW GCC 13.2 与 MSVC/AVX2 下均以 `2e-14` 对齐。性能矩阵的跨路径门限已从临时的 `2e-10` 收紧回 `2e-12`。

MSVC 复测中，128 点 uniform adaptive 约 `1.04 ms`（W1）/`1.02 ms`（W2），near-diagonal adaptive 约 `31.6 μs`/`29.9 μs`，没有出现与 matching materialization 相关的数量级回退。

## Exact duplicate / mass compression

`PreparedDiagram` 现在缓存 exact duplicate representatives 和 integer multiplicities。Wasserstein duplicate 分支把完全相同的 finite points 压成 group，在 source→row、row→column、column→sink 上使用整数容量，以 primal-dual Dijkstra 求 exact positive-profit flow；最终按 group flow multiplicity 从原始 cross/diagonal cost 稳定重算。

配置包含 `none/exact/adaptive`：

- `exact` 强制走 group flow，用于 oracle 和负面对照；
- `adaptive` 只有当 `unique_pair_count <= original_pair_count/16` 时启用 mass solver；
- 轻重复或无重复只读取 `PreparedDiagram` 已缓存的 group count，直接回到普通 adaptive kernel，不再每次排序或分配 group vectors。

MSVC，128×128，5 repetitions × 7–9 randomized rounds：

| 场景 | metric | adaptive no mass | adaptive mass | ratio | groups | removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| duplicate-heavy | W1 | 1392.1 μs | 19.3 μs | 0.014× | 8 | 248 |
| duplicate-heavy | W2 | 1197.2 μs | 12.1 μs | 0.010× | 8 | 248 |
| uniform | W1 | 1424.9 μs | 1441.1 μs | 1.011× | 256 | 0 |
| uniform | W2 | 1571.4 μs | 1496.7 μs | 0.952× | 256 | 0 |

强制 `exact` 在无重复 uniform 128 上需约 `15–18 ms`，证明“只要出现一个重复就压缩”会是严重负优化；因此保留 1/16 pair-ratio 门限。duplicate-heavy 512 的 group 数仍为 8，exact 路径约 `25–40 μs`，而未压缩 adaptive 约 `71 ms`。默认 `WassersteinConfig` 已启用 duplicate `adaptive`，显式 candidate/graph/matcher 实验仍关闭 duplicate 分支，保证正交对照。

## 已淘汰或限制使用的路径

- SPFA residual solver：MSVC 数值循环，淘汰；
- sparse solver on clustered/dense：通常慢一个数量级，只保留显式实验配置；
- unconditional component split：dense 图纯 overhead，adaptive 在 `ρ>=0.05` 时跳过；
- unconditional global greedy：中高密度排序成本大，只在 tiny/低密度图尝试 certificate；
- full dense Hungarian on sparse graph：不删除 isolates，near-diagonal 明显负优化；
- W1 two-pointer row-list：没有稳定击败 binary sweep + CSR，不进入默认；
- fixed-degree、bitmask/lazy-weight、block-sparse：均完成独立实验，但没有击败对应端到端 adaptive winner，不进入默认；
- unconditional duplicate flow：无重复 uniform 会慢约一个数量级；只保留 exact 对照，默认使用 1/16 compression-ratio dispatcher；
- 固定只看 `n*m` 的 dispatcher：已删除。

## Parallel component 实验

新增 `parallel_sparse` 和 `parallel_dense` 两条 exact component strategy。component discovery 仍为串行；分解完成后，以最多 8 个 worker 通过原子索引调度独立 component，每个 worker 只写自己的 matching/stats，主线程按原 component 顺序合并，避免浮点归并顺序和 matching 冲突。小于 256 个总顶点或少于两个 component 时回退对应串行分支。

为避免只在 `1×1` tiny 图上测线程开销，benchmark 新增 `multi_component`：每个 component 固定为独立的 `32×32` dense positive-saving 子图。MSVC 同进程随机轮序摘要：

| N | metric | component dense | parallel dense | component sparse | parallel sparse | parallel dense / serial dense |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 256 | W1 | 0.878 ms | 1.300 ms | 7.77 ms | 3.09 ms | 1.48× |
| 256 | W2 | 0.694 ms | 1.354 ms | 7.82 ms | 3.18 ms | 1.95× |
| 512 | W1 | 1.506 ms | 2.212 ms | 17.00 ms | 4.87 ms | 1.47× |
| 512 | W2 | 1.187 ms | 1.678 ms | 15.77 ms | 4.47 ms | 1.41× |
| 1024 | W1 | 2.592 ms | 3.731 ms | — | 8.15 ms | 1.44× |
| 1024 | W2 | 2.223 ms | 2.469 ms | — | 6.78 ms | 1.11× |
| 2048 | W1 | 5.415 ms | 5.509 ms | — | 14.91 ms | 1.02× |
| 2048 | W2 | 5.653 ms | 6.408 ms | — | 13.66 ms | 1.13× |

parallel sparse 相对串行 sparse 在多个中型分量上可快约 1.9–3.5×，但这个局部对照不是端到端 winner：串行 component dense 又比 parallel sparse 快 2.4–4.6×。parallel dense 直到 2048 W1 才接近持平，所有稳定配置仍未胜出。512 个 `1×1` adversarial-sparse component 上，线程化还会把约 0.18–0.20 ms 放大到 0.72–0.78 ms。结论：parallel component 工程实验完成，两条配置保留为显式 exact 对照，不进入 adaptive；当前赢家是串行 component+dense。

该实验同时暴露了原 adaptive 的结构盲点：它在全局 density `>=0.05` 时无条件跳过 component，导致 512 点、16 个独立 `32×32` dense component 的输入走全局 sparse/dense matcher，耗时约 `274–307 ms`。component skip 阈值现已与 dense matcher 的 `ρ=0.15` 对齐；density `0.062–0.125` 会先分解，再按局部 density 选择 dense SAP，且局部 dense component 关闭无收益的 global-greedy 排序。MSVC 512 点复测：

| metric | adaptive before | adaptive after | explicit component dense | after / explicit |
| --- | ---: | ---: | ---: | ---: |
| W1 | 307.4 ms | 1.792 ms | 1.699 ms | 1.05× |
| W2 | 296.3 ms | 1.390 ms | 1.232 ms | 1.13× |

同一 dispatcher 修改不影响 density `0.212–0.220` 的 uniform 512，仍跳过 component 并走原 dense matcher。该结构修正进入 adaptive。

## 尚未完成的第二轮

以下仍是 `../proposals/wasserstein第二阶段.md` 的未完成项，不能把当前状态称为完整 Wasserstein 优化结束：

1. 更高效的增量最短路/负环算法与 E4 adaptive pricing（matching、residual 和 scratch 持久化 baseline 已完成，确认容器复用不是瓶颈；active-set 各路径因端到端负优化不进入 adaptive）；
2. GUDHI/POT exact 和 Hera exact 四层外部基线；
3. real persistence diagrams。

这些完成并重新跑完整 8–8192 median/p95 矩阵后，才能按第二阶段文档的标准宣告整体完成。
