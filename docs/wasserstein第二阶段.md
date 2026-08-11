重新按你 Bottleneck 那轮的标准来设计。之前的错误在于把“DSR + 两个 sweep”当成两个完整 Phase；实际上它们只能算 **两个独立实验分支**。

Wasserstein 应该和 Bottleneck 一样走：

[
\boxed{
\text{大量独立实验}
\rightarrow
\text{逐项 benchmark}
\rightarrow
\text{淘汰负优化}
\rightarrow
\text{组合赢家}
\rightarrow
\text{adaptive default kernel}
}
]

而不是“实现一个理论点就进入下一阶段”。

---

# 一、先重新定义最终目标

Bottleneck 现在的默认内核大致是：

```text
clipped candidates
+ quickselect
+ adaptive adjacency
+ adaptive matcher
+ large-N special branch
+ prepared diagram / workspace / batch
```

Wasserstein 最后也应该形成同级别的东西：

```text
DSR / geometric candidates
+ adaptive weighted graph representation
+ adaptive weighted matcher
+ large-N active-set/pricing
+ component decomposition
+ prepared diagram / workspace / batch
```

数学核心统一写成

[
W_{q,p}(X,Y)^q
==============

D-S^\star,
]

其中

[
D=\sum_i a_i+\sum_j b_j,
]

[
s_{ij}=a_i+b_j-c_{ij},
]

[
S^\star
=======

\max_M
\sum_{(i,j)\in M}s_{ij}.
]

只需要考虑

[
E^+={(i,j):s_{ij}>0}.
]

所以 Wasserstein 最终是在解：

[
\boxed{\text{positive-saving weighted bipartite matching}}
]

---

# 二、第一轮完整优化：Candidate 实验组

你目前完成的工作只算这里面的两个实验。

## A0. Dense baseline

必须保留一个干净基线：

```text
all-pairs cost
→ dense saving matrix
→ dense Hungarian
```

以后所有实验相对这个比较。

---

## A1. DSR block scan ✅ 已做一部分

直接计算

[
s_{ij}=a_i+b_j-c_{ij}
]

并过滤

[
s_{ij}\le0.
]

但是这里还应该拆：

```text
scalar DSR
blocked DSR
branchless DSR
SIMD DSR
```

而不能“DSR 写了”就结束。

---

## A2. W1-(L^\infty) rotated sweep ✅ 原型已有

利用

[
s_{ij}
======

2\min(v_i,v_j)-|u_i-u_j|.
]

候选条件：

[
|u_i-u_j|
<
2\min(v_i,v_j).
]

这里至少应该继续实验：

```text
two-pointer sweep
binary-search window
blocked sweep
branchless sweep
sorted-index reuse
small-N brute / large-N sweep crossover
```

---

## A3. W2-(L^2) rotated sweep ✅ 原型已有

[
s_{ij}
======

2v_iv_j-(u_i-u_j)^2.
]

候选：

[
(u_i-u_j)^2<2v_iv_j.
]

同样不能只保留一种 sweep。

---

## A4. SIMD block scan

这应该对应 Bottleneck 的 AVX2 实验。

例如 W2-L2 热循环实际上只有：

[
du=u_i-u_j,
]

[
s=2v_iv_j-du^2.
]

非常适合一次算 4/8 个 (j)。

比较：

```text
scalar
AVX2
AVX-512（如果环境值得）
```

尤其 clustered 数据很重要。

因为 clustered 几乎全稠密：

[
\rho\approx1,
]

这时候 sweep 可能是负优化，反而 SIMD dense scan 很可能是正确路线。

---

# 三、第二组：Weighted adjacency 数据结构实验

这一步对应你 Bottleneck 的：

* 64/128/256-bit adjacency
* adaptive adjacency
* blocked representation

Wasserstein 虽然有权值，但仍然应该做同级别实验。

## B0. Dense matrix

基线。

---

## B1. CSR

```cpp
row_offsets
col_indices
weights
```

目标场景：

[
E\ll nm.
]

这是当前最应该马上实现的。

---

## B2. Edge-list per row

每行：

```cpp
struct Edge {
    uint32_t col;
    double saving;
};
```

比较 CSR 和 row-vector。

不要先假定 CSR 一定更快。

weighted solver 的访问模式可能更适合连续 row-edge list。

---

## B3. Fixed small-degree adjacency

对应 Bottleneck fixed 64/128/256-bit 的思想，但这里按 degree 做：

```text
inline 4 edges
inline 8 edges
inline 16 edges
inline 32 edges
overflow → heap/CSR
```

near-diagonal 数据只有 0.6% density 时非常值得测。

比如 (N=256)：

[
0.006\times256\approx1.5
]

平均每个 row 可能只有一两条边。

为了两条边建完整 CSR 不一定划算。

---

## B4. Bitset mask + lazy weight

对于：

[
N\le64,;128,;256
]

保存：

```text
positive-saving bitset
```

但不保存所有 weight。

solver 真正访问某条边时：

[
s_{ij}
]

即时重算。

这可能成为 Wasserstein 对应 Bottleneck fixed-bitset 的方案。

因为

[
s_{ij}
]

在 W1-L∞ / W2-L2 下非常便宜。

需要实测：

[
\text{load weight from memory}
\quad vs\quad
\text{recompute weight}.
]

---

## B5. Block sparse

例如 (16\times16) 或 (32\times32) tile。

一个 tile 内如果存在候选，再存 tile。

适合候选具有几何局部性的 diagram。

这和普通 CSR 的 cache 行为可能完全不同。

---

# 四、第三组才是真正的核心：Weighted Matcher 大实验

这是现在最大的缺口。

不能只“写一个 sparse primal-dual”。

至少要同时做下面这些独立实验。

## C0. 当前 dense Hungarian

保留 oracle / baseline。

---

## C1. Dense shortest augmenting path / JV

Hungarian 未必是 small-N 的最快算法。

你应该直接比较：

```text
Hungarian
Jonker-Volgenant style
dense shortest augmenting path
```

这可能直接改善 clustered 和 uniform。

换句话说：

> 即使 DSR 完全失效，Wasserstein 也必须有一个比当前 dense baseline 更强的 dense kernel。

这是之前方案漏掉的关键。

---

## C2. Sparse shortest augmenting path

复杂度真正开始依赖：

[
E
]

而不是 (nm)。

这是 near-diagonal 0.6% density 最应该吃到收益的地方。

---

## C3. Sparse primal-dual

维护：

[
\alpha_i,\beta_j
]

和 reduced cost。

和 C2 单独 benchmark。

不要理论上认为谁更优。

---

## C4. Positive-profit flow

因为我们的边权全部满足：

[
s_{ij}>0.
]

可以构造：

```text
source → X → Y → sink
```

每次只寻找存在正收益的 augmenting path。

当最优剩余增广收益

[
\le0
]

直接结束。

这可以作为 Wasserstein 的一个 large sparse solver 分支。

它某种意义上对应你 Bottleneck 的 `mandatory-flow` 独立实验。

---

# 五、Component 实验必须独立做

对应 Bottleneck 的 component branch。

令

[
G=(X\cup Y,E^+).
]

如果

[
G=C_1\cup\cdots\cup C_k,
]

则

[
S^\star
=======

\sum_{r=1}^kS^\star(C_r).
]

所以 component decomposition 是严格 exact。

但要单独测试：

```text
no component split
component + dense
component + sparse
component + adaptive
```

而不是默认打开。

因为 clustered 上 DFS/BFS 本身可能就是纯 overhead。

---

# 六、Tiny component 专用内核

这一步很可能非常有价值。

DSR 后可能出现大量：

```text
1×1
1×k
2×2
2×3
3×3
```

component。

这些根本没必要进通用 matcher。

例如 (1\times k)：

[
S^\star=\max_j s_{1j}.
]

(1\times1)：

[
S^\star=s_{11}.
]

小 component 可以：

```text
closed form
stack-local exhaustive
tiny dense solver
```

建议分别测：

```text
size <= 2
size <= 4
size <= 8
size <= 16
```

找到 crossover。

这非常像你 Bottleneck 的“小图直接 bitset + greedy Kuhn”。

---

# 七、Greedy 不用于求最终答案，而用于 warm start

weighted matching 不能像 Bottleneck Kuhn 那样简单。

但可以先做：

```text
row-max greedy
global descending-saving greedy
degree-first greedy
```

得到一个 initial matching。

然后 sparse primal-dual 从这个 matching 开始。

测试：

```text
cold solver
row-greedy warm start
global-greedy warm start
```

有可能减少 augmentations。

如果没有收益就淘汰。

---

# 八、Large-N 要有一个真正独立的路线

这才对应你 Bottleneck large sparse graph 的 `mandatory-flow`。

不要让 (N=4096/8192) 还必须：

[
\text{生成全部 }E^+
]

之后才求解。

## D1. Top-k seed graph

每个 row 先只加入 top-(k) saving：

[
k=2,4,8,16,32.
]

得到 restricted graph。

注意这本身**不是 exact**。

---

## D2. Restricted solve

求出 primal + dual：

[
\alpha_i,\beta_j.
]

---

## D3. Pricing

检查遗漏边是否存在：

[
s_{ij}-\alpha_i-\beta_j>0.
]

如果不存在：

[
\boxed{\text{当前 matching 已经全局 exact}}
]

否则加入 violated edges。

---

## D4. Incremental resolve

不要每轮重新 solve。

保留：

```text
matching
dual
workspace
candidate graph
```

继续增广。

这就是 Wasserstein 对应 Bottleneck `blocked incremental` 的路线。

---

# 九、Pricing oracle 也要做多条实验分支

不要默认 k-d tree。

## E1. Full SIMD pricing

大块扫描：

[
\max_j(s_{ij}-\beta_j).
]

dense/clustered 可能反而最好。

---

## E2. Rotated sweep pricing

W1-L∞ 和 W2-L2 可以利用你已有 sweep。

---

## E3. KD-tree pricing

generic metric 和 large sparse 数据使用。

节点维护类似：

[
\max_{j\in node}(b_j-\beta_j).
]

配合距离 lower bound 剪枝。

---

## E4. Adaptive pricing

最后根据：

[
N,\rho,\text{previous hit rate}
]

选择：

```text
SIMD scan
sweep
KD-tree
```

---

# 十、Duplicate / mass compression 单独实验

Persistence diagram 很可能有重复或高度相同的点。

实验：

```text
no compression
exact duplicate compression
small integer multiplicity solver
```

但必须单独报告：

* duplicate-heavy 数据
* normal 数据

不能因为 synthetic duplicate 数据提速就默认打开。

---

# 十一、PreparedDiagram 要对标 Bottleneck 的水平

不是只缓存 midpoint。

最终建议至少：

```cpp
birth
death

u
v

diag_cost_q

sorted_by_u

optional KD tree

duplicate/mass info
```

并测试：

```text
one-shot
one-to-many
many-to-many
```

因为 prepared diagram 的收益主要不是一次距离计算，而是：

[
X\text{ 与 }Y_1,\dots,Y_K
]

这种场景。

---

# 十二、工程实验组也要全部跑

直接对标 Bottleneck。

分别实验：

```text
workspace reuse
caller output buffer
native one-to-many
blocked candidate generation
AVX2
parallel component solving
parallel candidate generation
custom allocator / arena
LTO
PGO（后期）
```

每一个独立 benchmark。

LTO 如果还是 0.96×：

> 淘汰。

不要因为“理论上应该快”保留。

---

# 十三、最终 benchmark 矩阵必须比现在大得多

规模至少：

[
N=
8,16,32,64,128,256,512,
1024,2048,4096,8192.
]

条件允许加：

[
16384.
]

数据分布：

```text
uniform
near-diagonal
separated
clustered
duplicate-heavy
imbalanced n:m
adversarial dense
adversarial sparse
real persistence diagrams
```

还需要专门测：

[
n\ne m
]

例如：

```text
32 × 512
128 × 2048
512 × 4096
```

因为 persistence diagrams 很少保证等长。

---

# 十四、benchmark 不只报告 total time

每个配置必须拆：

[
T=
T_{\rm prepare}
+
T_{\rm candidate}
+
T_{\rm graph}
+
T_{\rm component}
+
T_{\rm solver}.
]

同时记录：

```text
candidate density
edge count
average degree
max degree
component count
largest component
augmentations
pricing rounds
number of priced edges
peak memory
median
p95
```

否则又会出现：

> “候选只剩 0.6%，为什么才 1.1×？”

但不知道时间究竟花在哪里。

---

# 十五、基线要有四层

不能只对当前 C++ baseline。

统一跑：

```text
1. naive/reference implementation
2. 当前 dense C++ kernel
3. GUDHI/POT exact
4. Hera exact
```

如果之后做 approximate：

```text
5. Hera delta=0.01
```

再加入。

你的优化目标最终不是：

[
\text{比自己第一版快}
]

而是：

[
\boxed{
\text{成为可独立使用的 Wasserstein kernel}
}
]

---

# 十六、正确性标准照 Bottleneck 的力度来

Small-N：

[
n,m\le 7\text{ 或 }8
]

做 exhaustive matching oracle。

然后大量：

```text
random differential
duplicate
empty diagram
single point
all diagonal-near
equal costs
ties
extreme coordinates
n ≠ m
essential points
```

再和 GUDHI/POT 对拍。

最好最后做到百万级甚至千万级随机 case。

不能因为性能分支多了就降低正确性标准。

---

# 十七、最后才产生 adaptive dispatcher

在所有独立实验结束之后再决定。

不是现在提前写死：

```cpp
if N < 128 ...
```

最终 dispatcher 应该看更多结构量：

[
N,\quad M,\quad E,\quad
\rho=\frac{E}{NM},
]

[
d_{\rm avg},\quad
d_{\max},\quad
C_{\max},
]

以及 metric：

[
(q,p).
]

最终可能形成：

```text
tiny
→ SIMD dense + dense JV

small sparse
→ bitmask/lazy weights + tiny/component solver

medium dense
→ SIMD dense + JV

medium sparse
→ sweep/CSR + sparse SAP

large sparse
→ lazy top-k + pricing + incremental solver

large dense
→ blocked SIMD dense solver

repeated queries
→ PreparedDiagram + reusable workspace
```

这才和 Bottleneck 的：

```text
adaptive adjacency
+
adaptive matcher
```

处于同一个层次。

---

# 重新划分你现在的进度

之前那份几分钟的结果，现在只能放到：

```text
Round A — Candidate experiments

A0 dense baseline                 ✅
A1 DSR                            ✅ 基础版
A2 W1-L∞ rotated sweep            ✅ 基础版
A3 W2-L2 rotated sweep            ✅ 基础版

A4 blocked scan                   ⬜
A5 SIMD scan                      ⬜
A6 sweep variants                 ⬜
A7 KD/range search                ⬜
A8 adaptive candidate generator   ⬜
```

也就是说如果对标 Bottleneck 那轮：

[
\boxed{\text{现在可能只有 10%-15%}}
]

比较合理。

---

# 我建议实际执行顺序

第一轮先把最有信息量的 12 个实验做完：

1. `dense Hungarian`
2. `dense JV/SAP`
3. `DSR + dense`
4. `AVX2 DSR + dense`
5. `sweep + dense`
6. `DSR + CSR + sparse SAP`
7. `sweep + CSR + sparse SAP`
8. `component + dense`
9. `component + sparse`
10. `tiny-component + sparse`
11. `greedy warm-start + sparse`
12. `adaptive dense/sparse`

这一轮结束后，你应该已经能回答三个关键问题：

[
\boxed{\text{small-N 最快是谁？}}
]

[
\boxed{\text{sparse-N 最快是谁？}}
]

[
\boxed{\text{dense/clustered 怎么避免负优化？}}
]

然后第二轮再专门攻：

```text
large-N active-set
top-k
pricing
incremental matching
KD-tree
duplicate mass
batch
workspace
```

这才是我认为真正**对标你 Bottleneck 优化流程**的 Wasserstein 方案。

最终完成标准也应该和 Bottleneck 一样：不是“某个 Phase 写完”，而是**每条有希望的独立优化路径都做过实验，负优化被明确淘汰，赢家被组合成一个经过多轮随机顺序 median/p95 验证的默认自适应内核**。这时才能说第一轮 Wasserstein 内核优化完成。
