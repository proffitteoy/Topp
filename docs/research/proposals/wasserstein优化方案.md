可以。并且我建议 Wasserstein 不要简单复制你现在 Bottleneck 的路线，而是把你 Bottleneck 里已经验证有效的几个思想——`prepared diagram / candidate clipping / adaptive dispatch / component decomposition / workspace reuse / one-to-many`——迁移过来，再针对 Wasserstein 的“加权匹配”重新设计核心。

当前 GUDHI 的精确 Wasserstein 路径，本质上先构造 ((n+1)\times(m+1)) 的稠密 cost matrix，再交给 POT 的 `emd/emd2`；Hera 路径则是 (\varepsilon)-scaling auction + weighted k-d tree，默认 `delta=0.01` 是近似算法。([GUDHI library][1]) Hera 的论文已经证明：单纯“把 auction 换成 k-d tree”这条路它做得相当深入，包括 Gauss–Seidel auction、价格加权 k-d tree和 epsilon scaling。

所以我建议你把新的 Wasserstein 内核设计成：

> **Diagonal-Savings Reduction + Adaptive Weighted Matching**
>
> 简称可以叫 **DSR-AWM**。

---

## 1. 第一刀不是换 solver，而是先把 Wasserstein 问题改写

这是整个方案最重要的地方。

设两个 persistence diagrams 为

[
X={x_1,\dots,x_n},\qquad
Y={y_1,\dots,y_m}.
]

定义

[
a_i=d_p(x_i,\Delta)^q,\qquad
b_j=d_p(y_j,\Delta)^q,
]

以及交叉匹配代价

[
c_{ij}=|x_i-y_j|_p^q.
]

那么 Wasserstein 的 (q) 次方实际上可以直接写成

[
W_{q,p}(X,Y)^q
==============

\min_M
\left[
\sum_{(i,j)\in M}c_{ij}
+
\sum_{i\notin M}a_i
+
\sum_{j\notin M}b_j
\right],
]

其中 (M) 只是 (X,Y) 离对角点之间的一个 partial matching；没被匹配的点直接扔给 diagonal。这就是 persistence-diagram Wasserstein 定义里的对角线机制。([GUDHI library][1])

现在定义“全扔给 diagonal”的 baseline：

[
D=\sum_i a_i+\sum_j b_j.
]

如果 (x_i) 和 (y_j) 相互匹配，而不是分别扔给 diagonal，那么节省的 cost 是

[
s_{ij}=a_i+b_j-c_{ij}.
]

于是得到一个非常关键的等价式：

[
\boxed{
W_{q,p}(X,Y)^q
==============

D-
\max_M
\sum_{(i,j)\in M}s_{ij}
}
]

而且

[
\boxed{
s_{ij}\le 0
\quad\Longrightarrow\quad
(i,j)\text{ 永远不需要进入最优匹配}
}
]

因为这种边还不如两边分别匹配 diagonal。

这相当于 Wasserstein 版的 `clipped candidates`。

你 Bottleneck 里面是：

[
d(i,j)\le\tau
]

才建边。

现在 Wasserstein 里变成：

[
\boxed{
c_{ij}<a_i+b_j
}
]

才建边。

这是**精确剪枝**，不是 heuristic。

所以第一阶段甚至可以完全不碰 auction。

---

# 2. 默认 W1-L∞ 可以再砍一刀

这个非常值得专门优化，因为：

```python
order=1
internal_p=inf
```

正是 GUDHI Wasserstein/Hera API 的默认形式。([GUDHI library][1])

对 persistence point

[
x=(b,d)
]

变换坐标：

[
u=\frac{b+d}{2},\qquad
v=\frac{d-b}{2}.
]

其中 (v) 就是 (L^\infty) 到 diagonal 的距离。

而且

[
|(b_1,d_1)-(b_2,d_2)|_\infty
============================

|u_1-u_2|+|v_1-v_2|.
]

于是 W1-L∞ 的 saving 直接变成

[
\begin{aligned}
s_{ij}
&=
v_i+v_j-
\left(
|u_i-u_j|+|v_i-v_j|
\right)\
&=
\boxed{
2\min(v_i,v_j)-|u_i-u_j|
}.
\end{aligned}
]

因此候选边存在的条件直接退化成

[
\boxed{
|u_i-u_j|<2\min(v_i,v_j)
}.
]

这太适合你现在 Bottleneck 已经写过的 `x-sweep` 体系了。

根本不用首先上二维 k-d tree。

把两个 diagrams 按 (u) 排序，然后 sweep：

[
|u_i-u_j|<2v_i
]

先缩出一个很小的窗口，再检查

[
|u_i-u_j|<2v_j.
]

通过后直接计算

[
s_{ij}=2\min(v_i,v_j)-|u_i-u_j|.
]

所以默认 Wasserstein kernel 可以变成：

```text
rotate to (u,v)
        ↓
sort by u
        ↓
x-sweep candidate window
        ↓
positive-saving filter
        ↓
weighted matching
```

这条路径我认为应该成为你的第一个 flagship kernel。

---

# 3. W2-L2 也有一个漂亮的专用内核

如果

[
q=2,\qquad p=2,
]

使用正交变换

[
u=\frac{b+d}{\sqrt2},
\qquad
v=\frac{d-b}{\sqrt2}.
]

那么

[
d_2(x,\Delta)=v,
]

并且

[
c_{ij}
======

(u_i-u_j)^2+(v_i-v_j)^2.
]

因此 saving 为

[
\begin{aligned}
s_{ij}
&=
v_i^2+v_j^2
-----------

\left[
(u_i-u_j)^2+(v_i-v_j)^2
\right]\
&=
\boxed{
2v_iv_j-(u_i-u_j)^2
}.
\end{aligned}
]

候选条件就是

[
\boxed{
(u_i-u_j)^2<2v_iv_j
}.
]

于是 W2-L2 也从一个“二维全部两两距离问题”变成了一个高度一维化的候选问题。

而且这里没有：

* `sqrt`
* `pow`
* Minkowski distance

热循环就是乘法、减法、比较。

非常适合 AVX2/AVX-512。

所以你的 Wasserstein library 第一批专用 fast path 应该就是：

```text
W1 + L∞
W2 + L2
W2 + L∞
generic q,p
```

而不是一开始所有 (q,p) 一视同仁。

---

# 4. Solver 不要只选一个，要继续走 Bottleneck 的 adaptive 路线

Bottleneck 里你现在已经证明“不同规模用不同 matcher”是对的。

Wasserstein 同样如此，但 matcher 要从 unweighted matching 换成 weighted matching。

我建议按照这个顺序做：

1. **Tiny / Small-N：dense primal-dual / JV。** 当 (N\lesssim64) 或 (128) 时，不建 k-d tree，不搞复杂 sparse structure。直接 SIMD 生成 savings matrix，然后用 shortest augmenting path / Jonker–Volgenant 风格的 weighted assignment。未匹配相当于 weight (0)。这一段目标就是干掉 GUDHI 的 Python/SciPy/POT overhead。

2. **Medium-N sparse：positive-saving sparse graph。** 只存

[
s_{ij}>0
]

的边，然后解 maximum-weight bipartite matching。可以用 sparse primal-dual augmenting path / min-cost flow。这里复杂度开始主要取决于

[
E=#{(i,j):s_{ij}>0},
]

而不是 (nm)。

3. **Connected components。** positive-saving graph 如果分成

[
G=C_1\cup C_2\cup\cdots\cup C_k,
]

每个 component 可以完全独立求 weighted matching：

[
S^*=\sum_r S^*(C_r).
]

这和你 Bottleneck 当前 component branch 的思想完全一样，而且 Wasserstein 这里甚至更自然：所有孤立点直接匹配 diagonal，不参与 solver。

4. **Large-N exact：不要完整 materialize 图，做 geometric column generation。** 这是我认为真正可能让你的实现和现有库拉开距离的部分。

---

# 5. Large-N exact：Geometric Column Generation

最大权匹配写成 LP：

[
\max\sum_{ij}s_{ij}z_{ij}
]

满足

[
\sum_jz_{ij}\le1,\qquad
\sum_iz_{ij}\le1,\qquad
z_{ij}\ge0.
]

它的 dual 可以写成

[
\min
\sum_i\alpha_i+\sum_j\beta_j
]

满足

[
\alpha_i+\beta_j\ge s_{ij}.
]

因此我们根本不必一开始生成所有 positive-saving edges。

先给每个点找到少数几个候选，比如 top-4/top-8：

```text
restricted graph
      ↓
weighted matching
      ↓
得到 dual α, β
```

然后检查有没有漏掉违反 dual constraint 的边：

[
\boxed{
s_{ij}-\alpha_i-\beta_j>0
}
]

如果没有，那么当前结果已经是**全图精确最优解**。

如果有，把这些 violated edges 批量加入，再解一轮。

所以：

```text
small candidate graph
       ↓
solve
       ↓
dual pricing
       ↓
find violated edges
       ↓
add edges
       ↓
solve again
```

直到

[
\max_{i,j}
(s_{ij}-\alpha_i-\beta_j)\le0.
]

这就给你一个 exactness certificate。

这比直接存 (O(nm)) cost matrix 更符合你现在 Bottleneck 内核的设计哲学。

---

# 6. Pricing oracle 再用 Hera 那套“geometry”，但用途不一样

Hera 的核心是在 auction 中不断寻找

[
\min_j(c_{ij}+p_j),
]

所以他们在 k-d tree 节点上维护价格下界，从而剪掉整棵子树。

你的 exact column generation 里可以做类似的事情，但搜索目标是

[
\max_j
\left(
s_{ij}-\beta_j
\right).
]

因为

[
s_{ij}
======

a_i+b_j-c_{ij},
]

相当于搜索

[
\max_j
\left(
b_j-\beta_j-c_{ij}
\right).
]

所以每个 k-d tree node 存：

[
\max_{j\in node}(b_j-\beta_j).
]

如果 query (x_i) 到这个 bounding box 的距离下界为 (L)，那么这个 node 的 reduced-saving 上界就是

[
a_i+
\max_{j\in node}(b_j-\beta_j)
-----------------------------

## L^q

\alpha_i.
]

若

[
\boxed{
a_i+
B_{\max}
-L^q-\alpha_i\le0
}
]

整棵子树直接 prune。

这其实是：

> **Hera weighted k-d tree + 你的 exact primal-dual certification**

但用途从 approximate auction 变成 exact pricing oracle。

这个方向我比“自己再写一遍 Hera”更看好。

---

# 7. 继续复用你 Bottleneck 已经完成的工程优化

你现在 Bottleneck 那套基础设施几乎都能留下：

`PreparedDiagram` 继续保留，只增加 `diag_distance`、旋转后的 (u,v)、按 (u) 的排序索引、flat k-d tree、duplicate multiplicity。距离计算继续 SoA；workspace、output buffer、one-to-many、native batch 全部继续用。essential points 在进入主 kernel 前单独拆分；GUDHI 当前也把 essential parts 分开处理，其中单有限坐标的 essential points 可以退化成 1D 排序匹配。([GUDHI library][2])

SIMD fast path 则直接专门写：

[
p=\infty,q=1:
\quad
c=\max(|dx|,|dy|)
]

[
p=\infty,q=2:
\quad
c=\max(|dx|,|dy|)^2
]

[
p=2,q=2:
\quad
c=dx^2+dy^2.
]

generic `pow()` 只留给非标准 (p,q)。

另外，重复 persistence points 可以压成

[
(x,y,\text{mass}),
]

不要复制成十几个相同顶点。Hera 自己也专门实现过 integer-mass auction，并报告在重复点较多时这样处理明显更有价值。

---

# 8. Approximate 分支最后再做

如果用户显式：

```python
delta > 0
```

再进入：

```text
Gauss-Seidel auction
+ epsilon scaling
+ weighted KD-tree
```

Hera 就是这一路，论文中他们使用逐轮降低 (\varepsilon)，并给出了相对误差终止条件；其几何版本在他们的实验中，随着 diagram 增大，相比非几何 auction 的优势显著扩大。

所以 approximate branch 不需要发明新理论。

你的改进点主要应该放在：

```text
adaptive brute-force vs KD
flat tree
SIMD leaves
prepared tree reuse
workspace reuse
integer mass
batch one-to-many
```

而不是重新研究 auction 本身。

更极端的 (10^4\sim10^5) 级 diagram，后期再考虑 quadtree / WSPD / graph sparsification 这种近似算法。已有工作确实在 persistence-diagram (W_1) 上走过 near-linear quadtree approximation 和 graph sparsification/min-cost-flow 路线。([arXiv][3]) 但我不建议现在就进这一层。

---

# 我给你的实际开发顺序

我会把整个 Wasserstein 内核分成六阶段：

```text
Phase 1
DSR 数学降维
D = Σdiag_cost
sij = ai + bj - cij
只保留 sij > 0

        ↓

Phase 2
默认 metric 专用化
W1-L∞ : rotated 1D sweep
W2-L2 : rotated 1D sweep
generic : SIMD / KD

        ↓

Phase 3
Adaptive exact solver
small dense JV
medium sparse primal-dual
component decomposition

        ↓

Phase 4
Large-N exact
lazy candidates
geometric pricing
column generation
dual optimality certificate

        ↓

Phase 5
工程层
PreparedDiagram
workspace reuse
one-to-many
batch
duplicate mass
output buffer
AVX2

        ↓

Phase 6
Approximate
Gauss-Seidel auction
epsilon scaling
weighted KD-tree
```

其中我认为**Phase 1 + Phase 2 是最应该先做的**。

因为这两步不是单纯“把现有代码写快一点”，而是直接改变问题的有效规模。

---

## 最终 dispatcher 我会设计成这样

```text
                    Wasserstein(X,Y)
                           │
                    PreparedDiagram
                           │
                 essential part split
                           │
                 Diagonal-Savings Reduction
                           │
                estimate candidate density
                 /          |           \
                /           |            \
         tiny/small      sparse        large
            │               │             │
        dense JV       sparse PD      lazy GCG
            │               │             │
            └──────── exact result ───────┘

delta > 0
    │
    └──── geometric auction
```

而 `internal_p/order` 再决定 candidate generator：

```text
q=1, p=∞
    → rotated x-sweep

q=2, p=2
    → rotated squared x-sweep

other common cases
    → SIMD block scan / adaptive KD

generic
    → weighted KD
```

这和你 Bottleneck 的最终结构会非常统一：

```text
Bottleneck
candidate generation
→ unweighted feasibility matching

Wasserstein
candidate generation
→ weighted optimal matching
```

底层的 `PreparedDiagram / geometry / SIMD / workspace / batch` 可以共享。

最关键的一点是：**不要把 Wasserstein 项目定义成“我要做一个比 Hera 更快的 auction”。**

更有价值的定位应该是：

> **一个 exact-first、adaptive、sparse persistence-diagram matching kernel；小数据用 dense SIMD，大数据用 diagonal-savings sparsification + geometric primal-dual，近似模式才走 Hera-style auction。**

这样你 Bottleneck 和 Wasserstein 最后是真的可以进同一个通用 C++ matching core，而不是两个互不相关的算法实现。

[1]: https://gudhi.inria.fr/python/latest/wasserstein_distance_user.html "Wasserstein distance user manual — gudhi v3.13.0 documentation"
[2]: https://gudhi.inria.fr/python/latest/_modules/gudhi/wasserstein/wasserstein.html "gudhi.wasserstein.wasserstein — gudhi v3.13.0 documentation"
[3]: https://arxiv.org/abs/2104.07710?utm_source=chatgpt.com "Approximation algorithms for 1-Wasserstein distance between persistence diagrams"
