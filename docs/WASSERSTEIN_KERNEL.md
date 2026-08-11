# Wasserstein 内核实验

记录日期：2026-08-11

## 当前边界

本阶段只实现并验证 `docs/wasserstein优化方案.md` 的 Phase 1 与 Phase 2，不处理 Python 外壳、绑定、wheel、近似 auction、large-N column generation 或通用 `(order, internal_p)`。

当前 exact 内核覆盖：

- `W1-L∞`；
- `W2-L2`；
- finite、empty、diagonal、duplicate 与 essential points；
- dense DSR reference generator；
- rotated one-dimensional sweep generator；
- exact dense primal-dual weighted assignment；
- DSR 后没有正 saving 边的行列剔除；
- 小矩阵 dense、较大矩阵 sweep 的实验 dispatcher。

公共内核声明位于 `include/bottleneck/wasserstein.hpp`，实现位于 `src/wasserstein.cpp`。它们不会改变现有 Bottleneck Distance 的默认配置或调用路径。

## 数学约化

有限点分别匹配到 diagonal 的总成本为

```text
D = Σ ai + Σ bj
```

交叉边 `(i, j)` 的 saving 为

```text
sij = ai + bj - cij
```

所以 exact 目标等价于

```text
Wq,p(X, Y)^q = D - max_M Σ(i,j in M) sij
```

`sij <= 0` 的边可以严格删除。未出现在正 saving 图中的点直接匹配 diagonal，不进入加权求解器。

使用 `u=(birth+death)/2`、`v=(death-birth)/2` 后，两条专用路径为：

```text
W1-L∞: sij = 2 min(vi, vj) - |ui-uj|
W2-L2: sij = 4 vi vj - 2(ui-uj)^2
```

两者都按 `u` 排序并用一维窗口生成候选；dense 路径直接从原始 birth/death 计算 cost，作为独立差分 reference。

## 正确性验证

构建与测试：

```powershell
cmd.exe /d /c scripts\build-kernel.cmd
build\manual\wasserstein_core_tests.exe
build\manual\bottleneck_core_tests.exe
```

当前测试包括：

- 手工可计算的 empty、identity、shift、diagonal、zero-saving boundary 与 essential cases；
- tiny diagram 对穷举 partial matching 的随机差分；
- dense DSR、rotated sweep、adaptive 三路径的 1,000 组随机差分；
- 交换两个 diagram 的对称性；
- 原有 Bottleneck 全量内核回归，确认共享 `PreparedDiagram` 扩展没有改变结果。

2026-08-11 本机结果：

```text
wasserstein_core_tests: all checks passed
bottleneck_core_tests: all checks passed
```

## 初步性能结果

命令：

```powershell
build\manual\wasserstein_core_bench.exe --repetitions 30
```

以下是同一进程内 C++ kernel microbenchmark 的探索性结果。`speedup` 为 dense/sweep；它不是 Python API 或端到端结论。

| 分布 | 规模 | metric | sweep 候选密度 | 正 saving 密度 | dense/sweep |
| --- | ---: | --- | ---: | ---: | ---: |
| uniform | 128 | W1-L∞ | 31.7% | 22.1% | 1.00× |
| uniform | 128 | W2-L2 | 30.8% | 21.2% | 1.10× |
| uniform | 256 | W1-L∞ | 32.5% | 22.4% | 1.02× |
| uniform | 256 | W2-L2 | 31.4% | 21.4% | 0.98× |
| near-diagonal | 128 | W1-L∞ | 0.6% | 0.4% | 1.11× |
| near-diagonal | 128 | W2-L2 | 0.6% | 0.4% | 1.17× |
| separated | 128 | W1-L∞ | 0.0% | 0.0% | 2.03× |
| separated | 128 | W2-L2 | 0.0% | 0.0% | 2.21× |
| separated | 256 | W1-L∞ | 0.0% | 0.0% | 2.02× |
| separated | 256 | W2-L2 | 0.0% | 0.0% | 2.60× |

clustered 输入的候选和正 saving 密度接近 100%，说明 sweep 不具备稀疏化收益。uniform 输入虽然删除约 78% 的边，但总耗时仍主要由 dense Hungarian 求解支配，因此总体加速接近 1×。

## 当前限制与下一轮实验

1. Savings matrix 仍按 `n*m` 分配；rotated sweep 减少 cost 计算和活跃点，但没有消除完整矩阵存储。
2. 活跃点很多但正边很少时，dense Hungarian 仍然是主要瓶颈。
3. 当前 dispatcher 的 `1024` pair 阈值只来自本机初步校准，不是稳定跨平台结论。
4. 尚未加入与固定版本 GUDHI/POT 的外部 oracle 差分；当前证据边界是手工值、穷举 oracle 和内部独立路径差分。
5. 尚未覆盖 generic `(q,p)`、integer mass、one-to-many workspace 和 approximate 模式。

下一轮应进入 Phase 3，但仍保持纯 C++ 内核实验：将 positive-saving edges 改为 CSR/edge list，比较 dense Hungarian、sparse primal-dual 和 component-wise weighted matching，并在相同四类分布上按 `(n, density)` 校准 dispatcher。
