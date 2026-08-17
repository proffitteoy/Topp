# 基准测试约定

本目录保留覆盖 1.0 默认路由、批量接口和关键负对照的 C++ exact 内核 benchmark。一次性实验程序与原始结果不属于发布树。

## 三层证据

1. **内核基准**：直接测量 C++ Bottleneck 算法，定位数据结构和算法开销。
2. **Python API 基准**：测量 `gudhi.bottleneck_distance`，包含数组转换和扩展模块调用。
3. **批量 API 基准**：测量 one-to-many、many-to-many 的吞吐、延迟、内存与并行扩展性。

三层结果分别报告，不能用第一层的倍数代替第三层结论。

## 公平比较最低要求

- 固定 GUDHI `e=0` oracle 版本；性能对照使用相同编译器、构建类型、依赖版本和编译选项
- 相同输入顺序和调用参数；exact 与未来的 approximate 接口分开报告
- 固定输入清单和随机种子，覆盖空图、单点、小图、阈值附近图和较大图
- 预热后重复测量，报告样本数、中位数及尾部延迟
- 先完成基线/优化版结果差分，再汇总性能
- 保存机器、CPU、线程设置和 wheel/二进制摘要

## 当前命令

```powershell
cmd.exe /d /c scripts\build-kernel.cmd
build\manual\bottleneck_core_bench.exe --repetitions 100
build\manual\bottleneck_grid_bench.exe --repetitions 20 --rounds 5
build\manual\bottleneck_batch_bench.exe --repetitions 3 --rounds 7
build\manual\bottleneck_large_bench.exe --repetitions 3
build\manual\bottleneck_multiplicity_bench.exe --repetitions 3 --rounds 7
build\manual\bottleneck_mandatory_bench.exe --repetitions 3 --rounds 5
build\manual\bottleneck_component_bench.exe --repetitions 3 --rounds 7
build\manual\bottleneck_geometry_bench.exe --repetitions 3 --rounds 7
build\manual\bottleneck_search_bench.exe --repetitions 3 --rounds 7
build\manual\bottleneck_router_bench.exe --repetitions 3 --rounds 7 --max-points 512
build\manual\wasserstein_core_bench.exe --repetitions 5 --rounds 7 --max-size 128
build\manual\wasserstein_batch_bench.exe --repetitions 5 --rounds 21
```

`bottleneck_grid_bench` 覆盖 uniform、near-diagonal、clustered、repeated、separated，以及对称/非对称 `8–256` 点输入，并输出 exact cross-edge density。多轮模式会打乱配置执行顺序并报告 median/p95；单轮输出只用于探索，不作为稳定回归阈值。

`bottleneck_large_bench` 覆盖 `32–4096` 点、稠密/稀疏/重复/分离分布和非对称输入，对比 geometric refinement、全候选 geometric matcher、默认 dispatcher 与旧 adaptive 路径。`--max-points N` 可限制交叉区实验；当前输出是固定种子多输入的均值，仍属于内核探索数据。

`bottleneck_multiplicity_bench` 固定 raw N 并扫描 exact duplicate ratio，对比 geometric default、显式 multiplicity capacity flow、当前 adaptive router 和旧 quickselect。配置逐轮随机顺序，输出 median/p95、unique group 数、移除点数和 capacity edge 数；可用 `--min-duplicate-ratio` 与 `--max-points` 缩小 crossover 扫描。

`bottleneck_mandatory_bench` 扫描 `128–2048` 点、`1:1–1:16` cardinality ratio、near-diagonal optional ratio，以及 uniform、clustered、separated 负对照。它独立消融 x-window candidate generation、dense mandatory flow、sparse mandatory flow 和二者组合，逐轮随机顺序并输出 median/p95、实际 `rho_delta`、候选/窗口/mandatory vertex/optional pair/flow edge/decision 统计。可用 `--min-points`、`--max-points`、`--ratio`、`--optional-ratio` 和 `--pattern` 缩小扫描。

`bottleneck_component_bench` 构造 `1–2048` 点的独立 cross components，扫描总规模、block size、组件数和 `s_max`，对比 geometric refinement、reusable Kuhn、完整增广图 component decomposition、mandatory sparse flow 与默认 router。结构统计在 exact threshold 上单独采集，计时配置逐轮随机顺序并报告 median/p95；可用 `--min-points`、`--max-points` 和 `--block-size` 缩小扫描。

`bottleneck_geometry_bench` 用 block size 精确控制 threshold cross-edge density（`block_size / N`），比较 KD matcher、geometric refinement、on-demand、x-sweep CSR、dynamic bitset 与默认 router。所有配置共享 x-window candidate generation，逐轮随机执行并报告 median/p95、KD visits、edge checks 和 augment searches；可用 `--min-points`、`--max-points` 和 `--block-size` 缩小扫描。

`bottleneck_search_bench` 比较 binary、GUDHI alpha、galloping、quickselect、严格单调 incremental、blocked incremental、geometric refinement 与默认 router。它在 uniform、clustered、near-diagonal、threshold-shell 上记录候选 raw/clipped/unique 数、重复率、decision/refinement rounds、增量 group additions/rollbacks 和 median/p95；可用 `--min-points`、`--max-points` 与 `--pattern` 缩小扫描。

`bottleneck_batch_bench` 扫描 `B=1,4,16,64,256,1024`，把 query/targets preparation、逐对 raw、只复用 prepared query、逐对全 prepared、batch 自分配输出和 batch 调用方复用输出缓冲区拆开统计。各调用变体逐轮随机顺序并报告 median/p95；可用 `--min-points`、`--max-points`、`--min-batch` 与 `--max-batch` 缩小扫描。

`bottleneck_router_bench` 是 E12 统一消融：在六类分布、对称与 1:8 cardinality ratio 上比较 adaptive、refinement、quickselect、geometric candidates、multiplicity 和 mandatory sparse，输出实际路由、cheap features、median/p95 与相对本轮最佳配置的 regret。配置逐轮随机执行；可用 `--pattern` 和 `--max-points` 缩小矩阵。

`wasserstein_core_bench` 覆盖 uniform、near-diagonal、clustered、separated、duplicate-heavy、imbalanced、adversarial-dense、adversarial-sparse 和 32×32 多分量专用输入的 `8–8192` 规模清单。配置按轮随机执行并报告 median/p95，同时拆分 prepare/candidate/graph/component/solver/pricing 时间，并记录 pricing rounds、priced/violated/materialized edges、各 oracle round 数、KD build/update 数、sparse arena build/scratch reuse/heap growth、dynamic inserted edges/Dijkstra runs/batch groups、max degree、peak arena bytes 与 peak graph bytes。可用 `--min-size`、`--max-size`、`--pattern`、`--metric`、`--repetitions` 和 `--rounds` 缩小实验矩阵；`--experiments sweep_csr_sparse,arena_sparse` 可随机轮序配对 clean sparse SAP 与 contiguous residual arena，`--experiments priced_dynamic_topk8,priced_dynamic_batched_topk8` 可配对逐边与按目标 column 批量化的 reduced-cost Dijkstra，`--experiments priced_topk8,priced_simd_topk8,priced_sweep_topk8,priced_kdtree_topk8,priced_kdtree_persistent_topk8,priced_adaptive_topk8,priced_incremental_topk8,priced_persistent_topk8,adaptive` 可配对 scalar full-scan、exact-guarded AVX2 full-scan、sweep、KD-tree、persistent KD geometry/scratch、E4 oracle router、matching-incremental、persistent-residual cycle-cancel 与默认配置，五档 top-k 名称均为 `2/4/8/16/32`。component 并行对照名为 `parallel_component_dense` / `parallel_component_sparse`；旧的单值 `--experiment` 仍兼容。超过 512 的普通分布只有在显式选择 priced 实验时才会开放，`multi_component` 只开放 component 专用配置，其他 dense baseline 继续跳过。

`wasserstein_batch_bench` 比较 one-shot、手写 prepared loop、native caller-buffer、显式 reusable workspace 和 native allocated-vector。五种模式先预热，再逐轮随机执行并报告 median/p95；固定顺序的单次数字不能用于判断 batch/workspace 收益。

## 结果归档

- 原始逐轮计时、临时诊断数据、探针和本地辅助脚本统一放入 `benchmarks/results/`；该目录默认忽略，不纳入版本控制。
- 可长期复现的正式 benchmark 应整理为 `benchmarks/` 下的源码；面向用户的最终汇总放入正式文档，并只保留理解表格所需的测试口径。
- 不要在仓库根目录或文档目录散放 JSON、NPZ、对象文件、可执行文件及其他实验中间产物。
