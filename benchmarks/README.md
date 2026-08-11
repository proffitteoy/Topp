# 基准测试约定

本目录同时包含单一均匀规模 microbenchmark、分布/非对称 Bottleneck 网格 benchmark，以及 Wasserstein Phase 1/2 内核 benchmark。当前阶段只测 C++ exact 内核。

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
build\manual\bottleneck_batch_bench.exe --repetitions 5
build\manual\bottleneck_large_bench.exe --repetitions 3
build\manual\wasserstein_core_bench.exe --repetitions 5 --rounds 7 --max-size 128
build\manual\wasserstein_batch_bench.exe --repetitions 5 --rounds 21
```

`bottleneck_grid_bench` 覆盖 uniform、near-diagonal、clustered、repeated、separated，以及对称/非对称 `8–256` 点输入，并输出 exact cross-edge density。多轮模式会打乱配置执行顺序并报告 median/p95；单轮输出只用于探索，不作为稳定回归阈值。

`bottleneck_large_bench` 覆盖 `32–4096` 点、稠密/稀疏/重复/分离分布和非对称输入，对比 geometric refinement、全候选 geometric matcher、默认 dispatcher 与旧 adaptive 路径。`--max-points N` 可限制交叉区实验；当前输出是固定种子多输入的均值，仍属于内核探索数据。

`wasserstein_core_bench` 覆盖 uniform、near-diagonal、clustered、separated、duplicate-heavy、imbalanced、adversarial-dense、adversarial-sparse 的 `8–8192` 规模清单。配置按轮随机执行并报告 median/p95，同时拆分 prepare/candidate/graph/component/solver/pricing 时间，并记录 pricing rounds、priced/violated/materialized edges、max degree 与 peak graph bytes。可用 `--min-size`、`--max-size`、`--pattern`、`--metric`、`--repetitions` 和 `--rounds` 缩小实验矩阵；`--experiments priced_topk8,priced_sweep_topk8,adaptive` 可在同一进程中随机轮序配对多个指定配置，旧的单值 `--experiment` 仍兼容。超过 512 的普通分布只有在显式选择 priced 实验时才会开放，其他 dense baseline 继续跳过。

`wasserstein_batch_bench` 比较 one-shot、手写 prepared loop、native caller-buffer、显式 reusable workspace 和 native allocated-vector。五种模式先预热，再逐轮随机执行并报告 median/p95；固定顺序的单次数字不能用于判断 batch/workspace 收益。

生成的原始结果放入 `benchmarks/results/`，该目录默认忽略。可复现脚本、脱敏固定输入或输入清单以及最终汇总应纳入版本控制。
