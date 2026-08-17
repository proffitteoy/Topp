# Bottleneck Phase 2 E11：one-to-many preparation 与 batch reuse

## 结论

本轮完成提案 E11 的独立 batch 消融，将 batch size 扩展为 `B=1,4,16,64,256,1024`，并把 query preparation、targets preparation、逐对 exact solve、结果向量分配和调用方输出缓冲区拆开计时。

结论是：现有 `PreparedDiagram` 值得由重复调用方显式复用，但 batch size 本身不需要改变 Bottleneck 的默认 pair router。

- `PreparedDiagram` 已缓存排序后的 births/midpoints、SoA 字段、diagonal distances、duplicate groups 和无限点序列；复用它能安全消除重复预处理；
- 对 16 点小图，预处理占比可见，全 prepared 相对 raw loop 的 median 收益约 19%–28%；
- 对 64/256 点图，exact candidate/matching 成为主成本，prepared 收益多在约 0%–6%，且小 batch 的差异容易落入系统噪声；
- `bottleneck_distances` 当前有意等价于 prepared pair loop；allocated-vector 与 caller-buffer 没有不同的 matching kernel，后者只避免结果向量分配；
- B 增大没有产生新的 algorithmic reuse：不同 target 的 candidates、threshold bracket、adjacency 和 matching state 都不安全共享。

因此本轮保留现有 batch API，不增加基于 B 的隐藏 router，也不引入跨 target matching state。若调用方会重复使用同一批 diagrams，应提前构造并持有 `PreparedDiagram`；一次性调用无需为“batch”额外建立更重的空间索引。

## Exact 边界

所有调用路径最终执行同一个 `bottleneck_distance(const PreparedDiagram&, const PreparedDiagram&)` exact kernel。benchmark 对 raw loop、prepared query、prepared loop、allocated batch 和 caller-buffer batch 的每个输出做逐 bit 比较。

跨 target 只共享 immutable prepared query，不共享候选半径、feasibility、matching 或 solver workspace。因此没有 E10 所审计的未来 threshold state，也没有一个 target 的图状态污染另一个 target。

## Benchmark

`benchmarks/bottleneck_batch_bench.cpp` 使用固定数据种子；五种调用变体在每轮随机排序，并报告 median/p95：

1. `raw_loop`：每个 pair 都从两份 raw diagrams 构造 prepared state；
2. `prepared_query`：只复用 query，每次重新 prepare target；
3. `prepared_loop`：query 与 targets 全部预处理后手写循环；
4. `batch_allocated`：原生 batch API 每次分配结果向量；
5. `batch_caller_buffer`：调用方复用结果缓冲区。

以下为 Windows x64、MSVC `/O2` 的代表 median。16 点使用 3 repetitions × 7 rounds；64 点使用 1 × 21 rounds；256 点 B≤256 使用 1 × 7 rounds，B=1024 使用 1 × 3 rounds。较重档位减少重复次数，但每个 timed block 已包含数百到数千次 exact pair solve。

### N=16

| B | query prepare | targets prepare | raw loop | prepared loop | batch caller-buffer | prepared speedup |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.003 ms | 0.003 ms | 0.032 ms | 0.026 ms | 0.026 ms | 1.26x |
| 4 | 0.002 ms | 0.012 ms | 0.124 ms | 0.102 ms | 0.102 ms | 1.23x |
| 16 | 0.004 ms | 0.051 ms | 0.671 ms | 0.528 ms | 0.519 ms | 1.27x |
| 64 | 0.005 ms | 0.317 ms | 3.078 ms | 2.336 ms | 2.482 ms | 1.32x |
| 256 | 0.005 ms | 1.153 ms | 11.687 ms | 9.827 ms | 9.198 ms | 1.19x |
| 1024 | 0.003 ms | 3.181 ms | 27.557 ms | 21.575 ms | 24.896 ms | 1.28x |

小图中 target preparation 是可测成本；若同一 targets 只使用一次，把 preparation 加回 prepared solve 后通常不会优于 raw one-shot。收益来自 prepared objects 的跨调用复用，而不是仅仅改用 batch 函数名。

### N=64

21 轮复测中，B=1/64/1024 的 raw → prepared loop median 分别为：

| B | query prepare | targets prepare | raw loop | prepared loop | caller-buffer |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.011 ms | 0.009 ms | 0.563 ms | 0.541 ms | 0.545 ms |
| 64 | 0.014 ms | 0.546 ms | 17.467 ms | 16.615 ms | 18.267 ms |
| 1024 | 0.006 ms | 8.192 ms | 324.539 ms | 319.075 ms | 318.474 ms |

B=4/16/256 的五种调用路径也都在约 ±7% 内交错，p95 明显宽于 median。没有证据支持声称 batch wrapper 或 caller buffer 能稳定加速 solver；能确认的是它们没有改变 exact 结果，且输出分配相对 pair cost 很小。

### N=256

| B | query prepare | targets prepare | raw loop | prepared loop | caller-buffer |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.041 ms | 0.038 ms | 6.314 ms | 6.903 ms | 6.043 ms |
| 16 | 0.037 ms | 0.939 ms | 111.315 ms | 97.275 ms | 93.415 ms |
| 64 | 0.043 ms | 2.072 ms | 321.741 ms | 329.838 ms | 314.488 ms |
| 256 | 0.040 ms | 8.213 ms | 1162.863 ms | 1141.963 ms | 1136.873 ms |
| 1024 | 0.022 ms | 37.593 ms | 5292.557 ms | 5003.919 ms | 5029.386 ms |

B=1024 时完整 target preparation 只占 prepared solve 的约 0.75%，query preparation 更可忽略；prepared loop 相对 raw 的 1.058x 收益没有随 B 形成新的数量级变化。大 N 的优化重点仍是单 pair candidate/matching kernel，而不是 batch dispatch。

## 实现与复现

本轮将正式 batch benchmark 改为：

- 支持 `--repetitions`、`--rounds`、`--min-points`、`--max-points`、`--min-batch`、`--max-batch`；
- 固定扫描提案要求的六档 B；
- 单独输出 query/targets preparation median；
- 五种调用变体随机轮序，输出 batch median/p95、per-pair median 和 raw-relative speedup；
- 所有变体与同一 expected vector 做逐 bit 检查。

复现实验：

```powershell
build\bn-phase2-e7-dev\bottleneck_batch_bench.exe --min-points 16 --max-points 16 --max-batch 1024 --repetitions 3 --rounds 7
build\bn-phase2-e7-dev\bottleneck_batch_bench.exe --min-points 64 --max-points 64 --max-batch 1024 --repetitions 1 --rounds 21
build\bn-phase2-e7-dev\bottleneck_batch_bench.exe --min-points 256 --max-points 256 --max-batch 256 --repetitions 1 --rounds 7
build\bn-phase2-e7-dev\bottleneck_batch_bench.exe --min-points 256 --max-points 256 --min-batch 1024 --max-batch 1024 --repetitions 1 --rounds 3
```

## Router 决策

- 不把 B 加入 single-pair adaptive router；
- batch API 继续接受 prepared query 和 prepared targets，并提供 allocated/caller-buffer 两种结果所有权方式；
- 不自动构建 Wasserstein state、KD index 或跨 target matching workspace；本轮没有证据证明其准备成本可被 B 稳定摊薄；
- 若未来出现同一 targets 多 query 或固定 geometry 的真实 workload，应另建 benchmark，不能从 one-query/many-target 结果外推。

## 剩余工作

1. E12：汇总 E6–E11 的正/负对照，统一审计 adaptive router；
2. E1–E5：独立验证 Wasserstein 数值信息，不把 prepared Bottleneck fields 当作 Wasserstein prepass；
3. 如需 batch workspace，必须先证明 pair 内 allocation 是主要瓶颈，并与当前 caller-buffer 路径随机轮序对照。
