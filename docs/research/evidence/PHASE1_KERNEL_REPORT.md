# Bottleneck exact 内核第一阶段实验报告

记录日期：2026-08-11

## 结论

第一阶段已形成一个保持 GUDHI `e=0` exact 语义的自适应 C++ Bottleneck Distance 内核。Small/Medium-N 继续使用 fixed/dynamic bitset、x-sweep、quickselect、Kuhn/mandatory-flow 等组合；较大或候选较密的输入进入本地 Hera-style 几何后端；近似只用于缩小候选区间，最终结果仍从精确 `double` 候选中选出。

该阶段没有 Python 外壳、wheel 或任何下游应用特化，也没有修改 `F:\GUDHI\gudhi-devel`。本文数字只代表本仓库 C++ 内核，不代表 Python 调用或相对 GUDHI Release 构建的端到端加速。

## 已完成的内核路线

### Small/Medium-N

- 预处理：有限/essential/对角点一次分类，AoS/SoA、对角距离、最大对角距离和 birth 排序进入 `PreparedDiagram`。
- candidate：全量、去重、合法上界裁剪、greedy 上界、quickselect、incremental 和 blocked incremental。
- adjacency：on-demand、dense byte、CSR、x-sweep、64/128/256-bit 和 dynamic bitset。
- matching：Kuhn、greedy/reusable/fixed Kuhn、constraint、component、mandatory-flow、HK 及 density-aware adaptive。
- 距离核：dense/recompute、AoS/SoA 和运行时检测的 AVX2。
- 工程路径：prepared diagram、native one-to-many、调用方提供输出缓冲区、普通和 LTO 构建。

### 几何后端

- 连续节点数组的二维 KD-tree，节点保存 bounding box。
- range query 后重新计算精确 `L∞` 距离；查询框使用向外 `nextafter`，避免浮点端点漏边或多边。
- greedy 初始化使用 subtree active count 和惰性删除，直接取得尚未占用的右侧点。
- 投影 dummy 块隐式枚举，不物化 expanded adjacency。
- 同一个 exact search 内保留 matching；阈值减小时裁掉失效边，再继续 augment。
- 几何 matcher 和 geometric refinement 都是显式实验开关；外部 Hera 不是运行时依赖。

### Approximation-assisted exact refinement

1. 以合法 lower/upper bound 启动几何 threshold oracle。
2. 得到一侧 infeasible、另一侧 feasible 的相对窄区间。
3. 只收集区间 `(lower, upper]` 内的对角距离和 cross `L∞` 精确候选。
4. cross candidate 通过已排序 birth 的 x-window 提取，再做 death 和精确距离复核。
5. 对精确候选二分并调用同一个 warm oracle，返回最小 feasible candidate。
6. 保留全候选 defensive fallback；任何未来空间索引回归都不能静默返回近似 upper endpoint。

公开默认语义仍是 exact。内部 refinement 的 1% 区间宽度不构成返回误差。

## 默认 dispatcher

默认 threshold 策略现在为 `adaptive`，判断只依赖通用的图规模和几何候选特征：

1. 总有限点数小于 128：保留 Small-N quickselect 路径。
2. 总点数至少 128 时，先计算“全部匹配对角线”上界下的 x-window pair 数。
3. x-window pair 为 0：cross match 在该上界内不可能存在，直接 exact 返回对角上界。
4. 较小 diagram 少于 32 点：保留旧路径，避免 8×128、16×128 一类窄侧输入的几何固定成本。
5. x-window pair 比例低于 1%：保留 quickselect/legacy matcher；该区间覆盖 near-diagonal 稀候选输入。
6. 其余输入：使用 geometric refinement。

这不是按单一 `n=512` 硬切换；显式 `quickselect`、`geometric_refinement` 和 `geometric_hopcroft_karp` 仍可独立消融。

## 正确性证据

- 原有 C++ 矩阵：250 组随机输入、24,000 种旧配置，以及空图、重复点、对角点和 essential points 边界，全部通过。
- 原有 GUDHI oracle：507 个输入 × 24,000 种旧配置，共 12,168,000 次 `e=0` exact 比较通过。
- 新增几何差分：300 组、每侧 0–40 点，包含重复点；4 个显式几何配置和默认 dispatcher 均与朴素 reference 相等。
- 新增 threshold contract：exact radius 返回 true，`nextafter(exact, 0)` 返回 false。
- 新增数值边界：`1e12` 坐标附近的一个 ULP 位移通过 exact/previous-threshold 检查。
- 新增 dispatcher 差分：64×64、32×256 和无 cross 的大图捷径与 reference 对拍。
- 新增 GUDHI oracle：507 个输入 × 5 个新增配置，共 2,535 次 `e=0` exact 比较通过。
- large benchmark 在可运行交叉区内同时计算 geometric refinement 和旧 exact 路径；所有被比较的 512–2048 点结果完全相等。

## 性能结果

环境为本机 MSVC 2022、`/O2`、单线程。下表是固定种子、每个 workload 3 个输入的探索性均值；不是跨机器稳定回归阈值。

| workload | adaptive | 旧 adaptive | 旧/新 | exact cross density |
|---|---:|---:|---:|---:|
| uniform 512×512 | 16.61 ms | 71.85 ms | 4.33× | 0.040 |
| uniform 1024×1024 | 74.89 ms | 491.99 ms | 6.57× | 0.019 |
| uniform 2048×2048 | 389.09 ms | 1.074 s | 2.76× | 0.012 |
| clustered 512×512 | 1.91 ms | 18.31 ms | 9.56× | 0.019 |
| clustered 2048×2048 | 30.27 ms | 450.64 ms | 14.89× | 0.013 |
| repeated 512×512 | 0.153 ms | 10.76 ms | 70.3× | 0.024 |
| repeated 2048×2048 | 0.512 ms | 343.69 ms | 671.9× | 0.024 |
| separated 512×512 | 0.015 ms | 3.84 ms | 255.8× | 0.000 |
| separated 2048×2048 | 0.057 ms | 96.91 ms | 1709× | 0.000 |
| uniform 512×2048 | 224.73 ms | 1.902 s | 8.46× | 0.167 |
| uniform 2048×512 | 227.27 ms | 1.791 s | 7.88× | 0.159 |

补充观察：

- uniform 4096×4096 的 adaptive 单次 exact 均值约 1.22 s；旧全候选路径未在该档运行。
- 4096×4096 的 3 个输入理论 raw candidates 合计约 5035 万，refinement 最终只保留 4,473 个精确候选，并避免 dense `n×m` 距离矩阵。
- near-diagonal 512/2048 点落回旧路径，和显式 legacy 处于同一性能区间；没有为追求 large-N 数字牺牲稀候选交叉区。
- 32×32 uniform 仍留在旧内核；64×64 uniform/clustered 开始进入几何路径。8×128、16×128 因较小侧限制继续走旧路径，32×256 可进入几何路径。

## 正优化与负优化

确认保留：

- warm matching + approximate interval + exact annulus refinement；
- KD range query 的向外边界与精确复核；
- active-count/lazy-delete greedy。相邻的单样本消融中，clustered 512 从约 3.53 ms 降到 1.55 ms，repeated 2048 从约 1.73 ms 降到 0.41 ms；
- 二维 dispatcher 和无-cross exact shortcut。

确认不进入默认：

- 强制标准 HK “找到 free-right 后只保留最短层”。单样本 uniform 4096 从约 1.9 s 退化到 3.4 s，原因是 BFS phase 和重复 KD 查询增加；已回滚。
- geometric matcher 配合全量 candidate/quickselect。512–2048 点通常比 refinement 慢 3–6×，只保留消融开关。
- LTO。旧 small-N 网格中位约 0.96×；本轮 large-N 单次结果同样分布依赖且噪声大，不改变默认构建。

## 证据边界和后续风险

- large-N 性能表是合成输入、单机、单线程 C++ 内核数据；尚无公平的 GUDHI 源码同编译器 large-N 对照。
- 4096 点仅对 geometric/adaptive 路径计时；旧路径因时间和内存成本未全量运行。
- 尚未加入峰值内存/分配计数，当前内存结论来自数据结构和 candidate 计数，不是进程 RSS 轨迹。
- 递归 augment 在极端构造图上的栈深仍需压力测试。
- Python 数组转换、绑定、GIL、wheel 和批量并行均未开始，也不能从本文倍数外推。

第一阶段的下一项内核工作应是分配/峰值内存仪表、构造型深增广压力测试和跨编译器复测；在这些完成前不进入 Python 外壳。
