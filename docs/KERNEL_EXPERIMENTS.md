# 内核实验登记

本文件只记录 C++ exact Bottleneck Distance 内核。Python 绑定、wheel、公共 API 与任何下游工作流不在当前性能结论内。

所有实验项必须满足：独立开关、与 reference 对拍、可单独 benchmark、负优化可保留为基线但不得进入默认路径。

## Reference contract

- 语义目标：有限点、essential points、空图、重复点与对角点均对齐 GUDHI `bottleneck_distance(..., e=0)`。
- reference variant：`sort_all + binary + dense_aos + on_demand + Kuhn + natural order`。
- threshold contract：exact radius 必须返回 true；`nextafter(exact, 0)` 必须返回 false。
- 默认优化配置：`sort_unique_clipped + quickselect + dense_aos + adaptive adjacency + adaptive matcher + degree ascending`。
- AVX2 是运行时检测后的可选距离内核，不是二进制的强制 CPU 前提。

## 已接入实验开关

| 维度 | 变体 | 当前结论 |
|---|---|---|
| candidate | `sort_all` | reference；保留 |
| candidate | `sort_unique` | 正确；重复值多时有效 |
| candidate | `sort_unique_clipped` | 保留；合法对角上界通常显著裁剪 cross candidates |
| candidate | `sort_unique_greedy_clipped` | 条件保留；重复点数据可胜出，普通随机图常因额外排序而变慢 |
| threshold | binary | reference 搜索；保留 |
| threshold | GUDHI alpha-style | 正确；小规模仍有竞争力 |
| threshold | exponential | 已验证负优化；当前分布 decision 数明显增加 |
| threshold | quickselect | 保留；避免全量候选排序，多数中大规模、稀疏和非对称数据胜出 |
| threshold | incremental | 已验证负优化；逐 unique weight 激活事件过多 |
| threshold | blocked incremental | 条件保留；按约 `sqrt(K)` 分块、命中块回滚重放，64 点事件数由 57,885 降至 3,556 |
| distance | dense AoS | 保留；当前稳定基线 |
| distance | dense SoA | 正确；未稳定胜过 AoS |
| distance | dense SoA AVX2 | 条件保留；运行时分派，较大或重复/稀疏输入有收益，小图可能被调用开销抵消 |
| distance | recompute AoS/SoA | 正确；小图偶有竞争力，中大规模通常不如 dense |
| adjacency | on-demand | reference；不分配邻接矩阵 |
| adjacency | dense byte | 正确；当前没有稳定优势 |
| adjacency | sparse CSR | 正确；部分较大稀疏图有竞争力 |
| adjacency | fixed 64/128/256-bit | 保留；按扩展图单侧大小选择固定字宽 |
| adjacency | dynamic bitset | 保留；固定字宽超限后的通用实现 |
| adjacency | bitset auto | 保留；自动选择 64/128/256 或 dynamic bitset |
| adjacency | x-sweep CSR | 条件保留；复用 birth 排序，在 x 窗口内再检查 y；总点数较大时获益 |
| adjacency | adaptive | 保留；小图 fixed bitset，较大图 x-sweep，并将较小 diagram 放到扫描侧 |
| matcher | Kuhn | reference；保留 |
| matcher | greedy + Kuhn | 保留；16–256 点区间整体强势 |
| matcher | reusable greedy Kuhn | 条件保留；thread-local buffer 减少重复分配，收益依分布变化 |
| matcher | fixed-array greedy Kuhn | 条件保留；扩展图单侧不超过 256 时消除 matching 数组分配，收益不稳定 |
| matcher | constraint + Kuhn | 已优化为 bitset、增量 degree、队列驱动；仍是负优化，保留实验基线 |
| matcher | component + Kuhn | 正确；当前合成分布构图开销大于收益，保留负基线 |
| matcher | Hopcroft–Karp | 正确；当前 small/medium-N 通常不如 greedy Kuhn |
| matcher | greedy + Hopcroft–Karp | 正确；比冷 HK 好，但通常仍不如 greedy Kuhn |
| matcher | mandatory flow | 条件保留；不显式建立投影 dummy 块，256 点极稀疏图胜出 |
| matcher | adaptive | 保留；大规模阈值图先固定成本采样，极稀疏时用 mandatory flow，否则用 reusable greedy Kuhn |
| vertex order | natural | reference；极小图可避免排序成本 |
| vertex order | degree ascending | 保留；中等规模通常更稳定 |
| batch | prepared scalar loop | reference batch 路径 |
| batch | native allocated output | 正确；消除未来语言边界的逐项调用，但 C++ 内部与 prepared loop 接近 |
| batch | native caller-provided output | 保留；避免结果 vector 分配，当前 C++ 单线程收益受噪声支配 |

## 已进入默认路径的无损优化

1. `PreparedDiagram` 一次保存 AoS、SoA、对角距离、最大对角距离和 birth 排序。
2. exact search 使用最大对角距离差作为合法 lower bound。
3. 使用“全部匹配对角线”的合法 upper bound 裁剪 cross candidates。
4. 物化扩展图时只访问可能存在的 cross/diagonal 边；数学上恒成立的 dummy 块直接写入，不再调用通用 edge predicate。
5. quickselect 对候选区间做三路划分，避免为二分搜索先完整排序。
6. adjacency 与 matcher 根据扩展图大小和阈值处采样密度分派。

## 正确性验证

- C++ 确定性边界：空图、单点、相同图、对角点、正/负 essential points、essential 数量不匹配。
- C++ 随机差分：固定种子，250 组、每侧 0–12 点，24,000 种配置与 reference 对拍。
- threshold decision：所有配置均验证 exact radius 成功、前一个可表示 `double` 失败。
- GUDHI oracle：507 个输入 × 24,000 种配置，共 **12,168,000 次** `e=0` exact 比较，全部通过。
- x-sweep 原型曾被 threshold contract 抓到端点窗口多接纳一条边；加入 birth 精确复核后重新通过完整 C++ 矩阵与 GUDHI 差分。

## 性能证据边界

`bottleneck_grid_bench` 覆盖：

- 对称 `8–256` 点；
- 非对称 `8×64`、`16×128`、`32×256` 及反向；
- uniform、near-diagonal、clustered、repeated、separated；
- exact 阈值处的原始 cross-edge density。

最近一次 20 对合成图网格中的 best-of-zoo 示例：

| workload | best mean | reference ratio | exact cross density | 主要胜出路线 |
|---|---:|---:|---:|---|
| uniform 8×8 | 4.53 μs | 0.254 | 0.205 | quickselect + bitset |
| uniform 64×64 | 373.43 μs | 0.052 | 0.123 | quickselect + fixed bitset + greedy Kuhn |
| uniform 128×128 | 1.82 ms | 0.035 | 0.096 | quickselect + adaptive adjacency |
| uniform 256×256 | 10.46 ms | 0.020 | 0.056 | quickselect；本轮 mandatory flow 胜出 |
| near-diagonal 256×256 | 2.01 ms | 0.007 | 0.008 | adaptive matcher -> mandatory flow |
| separated 256×256 | 1.03 ms | 0.009 | 0.000 | adaptive matcher -> mandatory flow |

这些数字只比较本仓库 C++ 变体，且是合成输入上的探索性 best-of-zoo，不是 GUDHI Release 公平性能对照，也不能代表 Python 单次调用或批量吞吐。后续已增加多轮交错顺序的 median/p95 输出；AVX2、workspace 与相近配置仍需在更多机器复测。

另一次 200 对 uniform 单规模 microbenchmark 中，当前默认 dispatcher 相对本仓库 reference 的结果为：4 点 `3.20x`、8 点 `3.34x`、16 点 `7.86x`、32 点 `11.61x`、64 点 `19.93x`。这组倍数只说明默认配置相对朴素 reference 的内核改进，不是相对 GUDHI 的加速。

额外工程实验：

- MSVC `/O2 /GL /LTCG` 对 145 个均匀 microbenchmark 行的中位加速为约 `0.96x`，即整体略慢且方差很大；LTO 不进入默认构建。
- native batch 在纯 C++ 中本质仍是同一求解循环。prepared diagram 对 16 点、20 项批次的一轮结果约减少 `1.77x` 的重复预处理成本；到 256 点后求解本身占主导，allocated/caller-buffer/native loop 基本接近。该结果只证明内核批量边界和预处理复用，不代表未来 Python 加速倍数。

## 下一批仍属内核的实验

1. 在 128–1024 点和更多密度档位校准 adaptive matcher 的规模与密度阈值。
2. 增加内存分配计数和峰值内存；继续压缩 ThresholdGraph 的 per-decision allocation。
3. 比较 Clang/GCC、平台原生 SIMD 多版本与 PGO；编译器结果和算法结果分开记录。
4. 继续评估 pair/triple Hall、bucket/grid 邻域和更紧合法 bounds；负优化同样登记。
5. 将 caller-provided solver workspace 扩展到 adjacency/candidate 缓冲区；native batch 与 caller-provided output 已进入内核。
