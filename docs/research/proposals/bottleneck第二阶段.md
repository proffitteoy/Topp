对，这个题应该按“研究项目”来设计，而不是先认定某个 Wasserstein 技巧一定能优化 Bottleneck。

核心问题应当改成：

[
\boxed{
\text{哪些 Wasserstein 信息能降低 exact Bottleneck 的总成本，在哪些数据结构下有效？}
}
]

并且必须把“共享 Wasserstein 的算法结构”和“真的先计算 Wasserstein”分开研究。Hera 本身就证明了两类距离可以共享几何加速思想：Bottleneck 用 geometric Hopcroft–Karp，Wasserstein 用 geometric auction，两者都通过几何邻域查询避免纯组合扫描。([应用数学学会][1]) 另外，Wasserstein 方向已有 condensation、multiplicity、min-cost flow、WSPD sparsification 等成熟思路，这些都值得作为 Bottleneck 的候选优化组件，而不是只研究 (p\to\infty)。([arXiv][2])

我建议整套实验按下面这个框架做。

---

## 一、先明确最终要得到什么

最后不是得到一个 `wasserstein_bottleneck()` 算法，而是得到一个自适应 portfolio：

[
\boxed{
\text{diagram features}
\longrightarrow
\text{router}
\longrightarrow
\text{best exact Bottleneck kernel}
}
]

例如最终可能变成：

| 输入结构              | 最终路径                                       |
| ----------------- | ------------------------------------------ |
| small-(N)         | 当前 Bottleneck，什么都不加                        |
| dense uniform     | dense/bitset Bottleneck                    |
| clustered dense   | dense matching，禁用 Wasserstein              |
| near-diagonal     | mandatory-flow + W-bound                   |
| duplicate-heavy   | multiplicity compression + W-bound         |
| separated         | 几何 certificate，禁用完整 Wasserstein            |
| 极度不等长             | rectangular/mandatory matching             |
| sparse large-(N)  | geometric search + components              |
| batch one-to-many | prepared diagram + amortized W information |

也就是说，论文/项目最终最好回答的是：

> Wasserstein-assisted Bottleneck 不是统一算法，而是一个根据 persistence diagram 几何结构自适应选择优化策略的 exact solver。

---

# 二、实验 0：冻结一个绝对基线

先别加任何新东西。

至少固定：

[
B_0=\text{当前 Topp exact Bottleneck}
]

以及外部参照：

[
B_G=\text{GUDHI exact},\qquad
B_H=\text{Hera}.
]

GUDHI 明确区分 `e=0` 的 expensive exact 路径和默认的有限精度快速路径。([GUDHI library][3])

但你现在已经发现 GUDHI `e=0` 在退化输入上存在方向相关错误，所以后续不能再把它作为唯一 oracle。你的 1280 对测试里唯一异常恰好就是这个问题，Topp 与剩余可靠 oracle 的 1279/1279 图对一致。

因此正确性 oracle 应该是：

[
\text{tiny brute force}
+
\text{Topp/Hera/GUDHI cross-check}
+
\text{metamorphic properties}.
]

特别测：

[
d(X,Y)=d(Y,X),\qquad
d(X,X)=0,
]

排列不变性、正比例缩放、同时沿 diagonal 平移、decision monotonicity。

后面每个优化实验都必须：

[
\text{optimization ON}
\quad\text{vs}\quad
B_0
]

单独比较。

否则最后根本不知道到底是哪项优化产生收益。

---

# 三、输入空间不能只测五种分布

你现在五种合成数据已经非常有价值：

* uniform
* near-diagonal
* clustered
* duplicate-heavy
* separated

因为它们已经表现出完全不同的算法性质。

但下一阶段建议扩成一个真正的“diagram phase space”。

### 规模轴

测：

[
N=8,16,32,64,128,256,512,1024,2048,4096,8192.
]

不要只测 (m=n)。

加入：

[
m:n=
1:1,\ 1:2,\ 1:4,\ 1:8,\ 1:16.
]

因为 rectangular matching 本身就是 Wasserstein 优化的重要来源，不等长图可能是新算法最容易拉开差距的地方。

### 几何结构轴

除了现有五类，再加：

| 分布                     | 要攻击的问题             |
| ---------------------- | ------------------ |
| small perturbation     | 两图非常相似             |
| diagonal-band          | 大部分点贴 diagonal     |
| off-diagonal           | 大量高 persistence 点  |
| lattice/grid           | 大量相同候选距离           |
| concentric clusters    | 高局部密度              |
| threshold-shell        | 大量距离集中在真正 (d_B) 附近 |
| adversarial dense      | 几乎所有点互为候选          |
| adversarial sparse     | 每点只有常数个可能邻居        |
| high dynamic range     | 距离跨很多数量级           |
| exact duplicates       | multiplicity       |
| near duplicates        | compression 不能直接使用 |
| asymmetric cardinality | partial transport  |
| essential points       | infinity 处理        |

其中 `threshold-shell` 很重要。

人为构造大量

[
c_{ij}\approx d_B.
]

它专门攻击 candidate threshold search。

如果你的算法在 uniform 上快，但这种输入突然崩掉，就说明优化依赖候选距离分散，而不是真正 general-purpose。

---

# 四、第一大方向：先研究 Wasserstein 与 Bottleneck 的“信息相关性”

这是整套项目最基础的实验。

先完全不优化 Bottleneck。

对每一对 persistence diagrams 同时求：

[
d_B,\quad W_1,\quad W_2,\quad W_4,\quad W_8,\quad W_{16},
]

这里 ground metric 必须统一用 (L_\infty)，否则是在混入另一个变量。

设非零匹配边数安全上界为 (N)，有

[
\frac{W_q}{N^{1/q}}
\le
d_B
\le
W_q.
]

如果 Wasserstein solver 返回匹配 (M_q)，更重要的是

[
d_B
\le
U_q
:=
\max_{e\in M_q}c_e.
]

因此真正研究：

[
L_q=\frac{W_q}{N^{1/q}},
\qquad
U_q=\max_{e\in M_q}c_e.
]

不要只看数值差。

记录四个指标：

[
\frac{L_q}{d_B},
\qquad
\frac{U_q}{d_B},
\qquad
T_{W_q},
\qquad
\frac{#C[L_q,U_q]}{#C},
]

这里 (C) 是全部 Bottleneck candidate radii。

最后那个指标尤其关键。

例如数值区间从

[
[0.1,10]
]

缩成

[
[0.9,1.1]
]

不一定重要。

真正影响 exact search 的是：

> 原来 100000 个候选半径，现在还剩几个？

所以应该定义

[
R_{\mathrm{cand}}
=================

1-
\frac{#C[L_q,U_q]}{#C}.
]

这是 Wasserstein 对 Bottleneck 的第一核心评价指标。

---

# 五、实验方向 A：完整 Wasserstein prepass

最简单：

[
X,Y
\rightarrow W_q
\rightarrow[L_q,U_q]
\rightarrow\text{exact Bottleneck}.
]

分别测

[
q=1,2,4,8,16.
]

实验不是为了证明 (q) 越大越好，而是找：

[
\boxed{
q^*=
\arg\min_q
(T_{W_q}+T_{B\mid[L_q,U_q]})
}
]

很可能 (q=8,16) 虽然界更紧，但 Wasserstein 自己更慢，最后反而不划算。

甚至很可能结论是：

[
q=1\text{ 或 }2
]

最好。

这本身就是有意义的实验结果。

---

# 六、什么时候完整 Wasserstein 明确不应该跑

你现有 benchmark 已经能先划掉很多情况。

把当前 Bottleneck 时间除以 Wasserstein 时间：

[
H=\frac{T_B}{T_W}.
]

这是完整 Wasserstein prepass 的“理论最大空间”。

因为即使后面的 Bottleneck refinement 变成 **0 ms**，总时间仍至少是 (T_W)。

根据你当前数据：

| (N) | 分布              |      (T_B/T_W) | 判断   |
| --: | --------------- | -------------: | ---- |
|   8 | 全部分布            | (0.45\sim0.98) | 完全禁用 |
|  32 | 全部分布            | (1.11\sim1.93) | 基本禁用 |
| 128 | uniform         |           1.20 | 禁用   |
| 128 | near-diagonal   |       **5.69** | 值得实验 |
| 128 | clustered       |           0.92 | 必须禁用 |
| 128 | duplicate-heavy |       **6.38** | 值得实验 |
| 128 | separated       |           1.13 | 禁用   |
| 512 | uniform         |           0.59 | 必须禁用 |
| 512 | near-diagonal   |      **22.75** | 强候选  |
| 512 | clustered       |           0.46 | 必须禁用 |
| 512 | duplicate-heavy |      **43.44** | 强候选  |
| 512 | separated       |           1.02 | 必须禁用 |

原始时间可以直接看到这一结构：512 点 Wasserstein 的 uniform/clustered 分别比当前 GUDHI 基线表现很差，而 near-diagonal、duplicate-heavy、separated 出现极大结构性收益。

但关键是 separated。

你的 Bottleneck separated 已经只有：

[
9.786\text{ ms},
]

Wasserstein 是：

[
9.592\text{ ms}.
]

 

即使 Wasserstein 算完以后 Bottleneck 一纳秒都不用，

[
\frac{9.786}{9.592}\approx1.02.
]

所以：

[
\boxed{\text{separated 虽然 Wasserstein 极快，却绝对不应该先算 Wasserstein。}}
]

这是后续 gating 实验必须体现出来的。

“Wasserstein 快”不等于“Wasserstein-assisted Bottleneck 值得”。

---

# 七、实验方向 B：只取 Wasserstein matching，不追求 Wasserstein value

这个方向我觉得比完整 (W_q) 更有潜力。

假设 Wasserstein 算法很快构造出一个 feasible matching (M)。

立刻有

[
d_B\le \max_{e\in M}c_e.
]

那么实验：

[
\text{greedy Wasserstein matching},
]

[
\text{partial SAP matching},
]

[
\text{full optimal Wasserstein matching}.
]

比较它们产生的 Bottleneck upper bound：

[
U_{\mathrm{greedy}},
\quad
U_{\mathrm{partial}},
\quad
U_{W}.
]

看是否存在一种现象：

> Wasserstein 只跑 5%–10% 的计算，就已经得到了和最终最优 Wasserstein matching 差不多的 Bottleneck upper bound。

如果成立，那你根本没必要完整算 (W_q)。

---

# 八、实验方向 C：Anytime Wasserstein

这个是上一项进一步发展。

Wasserstein solver 运行过程中不断记录：

[
t\mapsto L(t),
\qquad
t\mapsto U(t).
]

如果有合法 dual lower bound (D_q(t))，则

[
L(t)
====

\left(
\frac{D_q(t)}{N}
\right)^{1/q}.
]

一旦有完整 primal feasible matching (M(t))，

[
U(t)=\max_{e\in M(t)}c_e.
]

然后画：

[
t
\longmapsto
#C[L(t),U(t)].
]

这比“Wasserstein 最终多快”重要得多。

你真正要找的是：

[
\boxed{
\text{time-to-useful-bound}
}
]

而不是 time-to-optimal-Wasserstein。

可以测试停止条件：

[
R_{\mathrm{cand}}>50%,70%,90%,95%
]

时立即中止 Wasserstein，切回 Bottleneck。

---

# 九、实验方向 D：Wasserstein matching 做 Bottleneck warm start

这和 bounds 是两个独立实验，必须拆开。

给定 (M_W)，在 Bottleneck threshold (r) 下保留：

[
M_W(r)
======

{e\in M_W:c_e\le r}.
]

它天然是一个合法 partial matching。

然后从这里开始 augment。

比较：

| 方法 | 初始 matching                          |
| -- | ------------------------------------ |
| B0 | empty                                |
| B1 | greedy                               |
| B2 | previous threshold                   |
| B3 | Wasserstein matching surviving edges |
| B4 | Wasserstein + previous threshold     |

记录：

[
#\text{augmentations},
]

[
#\text{BFS/DFS edge visits},
]

[
\text{initial matched fraction},
]

[
T_{\mathrm{matching}}.
]

特别研究

[
S(r)=
\frac{|M_W(r)|}{|M_W|}.
]

如果在

[
r\approx d_B
]

时仍然有

[
S(r)>0.8,
]

Wasserstein matching 就是非常好的 warm start。

如果只有 20%，就没有必要。

---

# 十、实验方向 E：(q)-continuation / homotopy

这才是比较接近你最开始说的“从 Wasserstein 推 Bottleneck”。

做：

[
q=1\rightarrow2\rightarrow4\rightarrow8\rightarrow16.
]

但重点不是只测 (W_q)。

同时测 matching 如何变化：

[
M_1,M_2,M_4,M_8,M_{16}.
]

记录：

[
U_q=\max_{e\in M_q}c_e,
]

以及

[
\frac{U_q-d_B}{d_B}.
]

如果发现：

[
U_4\approx U_8\approx d_B,
]

那么没有必要继续到 16。

进一步研究：

[
M_q\rightarrow M_{2q}
]

的 warm start 是否有效。

最终可能得到：

[
W_1
\rightarrow W_4
\rightarrow d_B
]

而不是无限往上算。

这部分才是真正可以形成一个“continuation algorithm”的东西。

---

# 十一、实验方向 F：直接迁移 Wasserstein 的 multiplicity / flow

这个不需要计算 Wasserstein。

这是“算法结构迁移”。

Wasserstein 文献已经证明 condensation、multiplicity 与 transshipment/min-cost-flow 能显著减少图规模。([arXiv][2])

对于 Bottleneck，专门测：

[
n_{\mathrm{raw}}
\rightarrow
n_{\mathrm{unique}}.
]

定义

[
\rho_{\mathrm{dup}}
===================

1-
\frac{n_{\mathrm{unique}}}{n_{\mathrm{raw}}}.
]

实验：

* no compression
* hash exact duplicate compression
* sorted compression
* multiplicity-capacity flow

扫描：

[
\rho_{\mathrm{dup}}
===================

0,0.1,0.25,0.5,0.75,0.9,0.99.
]

找 crossover：

[
\rho_{\mathrm{dup}}^*.
]

例如最后可能得到：

[
\rho_{\mathrm{dup}}>0.2
\Rightarrow
\text{capacity kernel}.
]

这就成为 router 的规则。

---

# 十二、实验方向 G：partial transport / mandatory-flow

这个也属于 Wasserstein transport formulation 的迁移，而不是 W-value。

对于 threshold (r)，大量点如果能够送到 diagonal，就不需要参与完整 cross matching。

定义：

[
\rho_\Delta(r)
==============

\frac{
#{x:d(x,\Delta)\le r}
+
#{y:d(y,\Delta)\le r}
}{
m+n
}.
]

扫描不同 (\rho_\Delta)。

比较：

[
\text{full augmented graph}
]

和

[
\text{mandatory-only graph}.
]

预计：

* near-diagonal：强收益；
* dense off-diagonal：收益有限；
* 高 cardinality imbalance：可能强收益。

于是 runtime gating 可以使用：

[
\rho_\Delta(\hat r).
]

---

# 十三、实验方向 H：Wasserstein 几何稀疏化移植

这是另一条独立路线。

Hera 的主要经验就是用几何查询替代组合图中的线性扫描。([应用数学学会][1])

PDoptFlow 则进一步用了 WSPD 做 arc sparsification，不过它解决的是 approximate Wasserstein。([arXiv][2])

对 exact Bottleneck 不能直接随便删 edge。

但可以做安全的 speculative scheme：

[
G_1\subset G_2\subset\cdots\subset G.
]

先在稀疏图 (G_1) 中找匹配。

如果找到：

[
\text{feasible in }G_1
\Rightarrow
\text{feasible in }G.
]

这是安全 certificate。

如果失败：

[
G_1\text{ infeasible}
]

不能说原图 infeasible，只继续扩大：

[
G_1\rightarrow G_2\rightarrow\cdots.
]

比较：

* full adjacency
* x-sweep
* kd-tree
* grid
* WSPD-inspired neighborhood
* adaptive radius expansion

重点扫 candidate density：

[
\rho_E=\frac{|E|}{mn}.
]

最后应该找到：

[
\rho_E<\rho^*
\Rightarrow
\text{geometry},
]

[
\rho_E>\rho^*
\Rightarrow
\text{dense bitset}.
]

而不是永远使用 kd-tree。

---

# 十四、实验方向 I：component decomposition

这也应该单独做。

在当前 threshold graph 上计算：

[
G(r)=G_1\cup G_2\cup\cdots\cup G_k.
]

分别匹配。

研究两个变量：

[
k,
\qquad
s_{\max}
========

\frac{\max_i|G_i|}{|G|}.
]

如果

[
k\gg1,\qquad s_{\max}\ll1,
]

component 很有价值。

如果 clustered 数据最终形成一个 giant component：

[
s_{\max}\approx1,
]

做 component detection 本身可能就是额外开销。

因此也需要 gating。

---

# 十五、实验方向 J：Wasserstein bound 改变 Bottleneck 的搜索算法

这个方向容易被忽略，但很重要。

传统 exact Bottleneck 本质是在候选半径中寻找最小 feasible threshold。GUDHI 的定义本身就是寻找存在合法 perfect matching 的最小半径。([GUDHI library][3])

比较四种：

[
\text{binary search},
]

[
\text{ascending incremental sweep},
]

[
\text{galloping + binary},
]

[
\text{W-bound + ascending sweep}.
]

如果 Wasserstein 已经把候选范围缩成：

[
100000\rightarrow20,
]

那么 binary search 可能反而不是最好的。

直接从 lower bound 向上：

[
r_1<r_2<\cdots
]

增量加入 edges，并持续复用 matching。

这还能规避你在 GUDHI 里发现的状态复用问题：你定位到的错误正是 exact 搜索回退后继续复用了不合适的 matching state。

所以实验里必须比较：

[
\boxed{\text{monotone-state reuse}}
]

和

[
\boxed{\text{bidirectional search + rebuild}}
]

而不能再做不受约束的状态复用。

---

# 十六、batch 是另一维度，不能和 single-pair 混在一起

你现在的 benchmark 是：

[
1\text{ query}
\rightarrow64\text{ targets},
]

而且 Topp 时间包括 `prepare_diagram(query)` 和原生 batch API。

下一阶段必须同时测：

[
B=1,4,16,64,256,1024.
]

因为 Wasserstein preparation 的成本可能在 single-pair 下不划算，却能在 one-to-many 下被摊薄。

比如：

[
T=
T_{\mathrm{prepare}}
+
B,T_{\mathrm{pair}}.
]

如果 (B=1)，prepared spatial index 可能是负优化。

如果 (B=256)，它可能极其划算。

所以最终 router 甚至要包含：

[
\text{batch size}.
]

---

# 十七、最终 runtime router 应该看什么

不要一开始上机器学习。

先收集廉价 feature：

[
\boxed{
F=
(
n_{\mathrm{raw}},
n_{\mathrm{unique}},
m/n,
\rho_{\mathrm{dup}},
\rho_\Delta,
\hat\rho_E,
k,
s_{\max},
\text{bbox separation},
B
)
}
]

其中：

* (n_{\mathrm{unique}})：压缩后规模；
* (\rho_{\mathrm{dup}})：重复率；
* (\rho_\Delta)：near-diagonal 程度；
* (\hat\rho_E)：candidate graph density 的抽样估计；
* (k,s_{\max})：component 情况；
* bbox separation：两图是否明显分离；
* (B)：batch size。

第一版做手写 rule router。

数据够多以后，再做第二个实验：

[
F\rightarrow
\arg\min_j T_j.
]

用一个很小的 decision tree 去预测 kernel。

然后比较：

[
\text{oracle router}
]

[
\text{hand router}
]

[
\text{learned router}.
]

这里的 oracle router 是事后知道所有算法时间以后选择最快，只用来给出理论上限。

如果：

[
T_{\mathrm{hand}}
\approx T_{\mathrm{oracle}},
]

根本没必要上 ML。

---

# 十八、最终实验矩阵

最后我建议正式整理成下面十个方向：

| 编号  | 实验                       | 回答的问题                                |
| --- | ------------------------ | ------------------------------------ |
| E1  | (W_q) bound quality      | Wasserstein 界到底紧不紧                   |
| E2  | full W prepass           | 完整算 Wasserstein 值不值得                 |
| E3  | anytime W                | 能否只跑一小部分                             |
| E4  | W matching warm-start    | matching 能否迁移                        |
| E5  | (q)-continuation         | (W_q\to W_\infty) 是否工程上成立            |
| E6  | multiplicity             | duplicate-heavy 如何处理                 |
| E7  | mandatory partial flow   | near-diagonal / unequal size         |
| E8  | geometric sparsification | sparse large-(N)                     |
| E9  | components               | separated / clustered                |
| E10 | search strategy          | W bounds 是否改变 exact threshold search |
| E11 | batch reuse              | one-to-many 是否值得预处理                  |
| E12 | adaptive router          | 最终哪些优化应该开                            |

其中 E1–E5 才是真正的“利用 Wasserstein 数值信息”。

E6–E11 是“把 Wasserstein solver 里已经验证过的算法结构迁移到 Bottleneck”。

这两类必须分开汇报。

---

# 十九、我现在对哪些场景值得做的预判

根据你已有结果，我会先这样排优先级。

**第一优先：512+ duplicate-heavy。**

这是最漂亮的 Wasserstein-assisted 候选：

[
T_B=709.542\text{ ms},
\qquad
T_W=16.332\text{ ms}.
]

有超过 (43\times) 的理论 prepass 空间。 

这里一起测：

[
\text{multiplicity}
+
W_1\text{ bounds}
+
W\text{ warm-start}.
]

**第二优先：512+ near-diagonal。**

[
262.888\text{ ms}
\quad\text{vs}\quad
11.557\text{ ms}.
]

 

重点测：

[
\text{mandatory flow}
+
\text{anytime Wasserstein}.
]

**第三优先：dense uniform / clustered。**

但这里不是研究“怎么使用 Wasserstein prepass”，而是研究：

[
\boxed{\text{如何自动识别并禁止 Wasserstein。}}
]

因为 512 uniform：

[
T_W>T_B,
]

512 clustered 也是一样。

这些 case 是 router 的 negative control。

**第四优先：separated。**

也是 negative control。

不是因为 Wasserstein 慢，而是 Bottleneck 本身已经太快。

这里应该研究：

[
\text{geometric separation certificate},
]

而不是跑一次完整 Wasserstein。

**small-(N)** 基本也作为 negative control。

8 点、32 点你现在主要赢在低常数和 Python/API 工程，额外运行 Wasserstein 极难回本。当前 Bottleneck 在 8/32 点已经约比 GUDHI 快 (9.70\times/9.98\times)。

---

最后，我会把整个项目的核心判据定成一句话：

[
\boxed{
\text{优化是否启用}
\iff
\text{预计节省的 Bottleneck 工作量}

>

\text{该优化自身的成本}.
}
]

而不是“这个数学技巧理论上可以用于 Bottleneck，所以就启用”。

如果按实际开发顺序，我建议就是：

1. 先写 E1：扫描 (W_q) 的 bound tightness，不改 Bottleneck；
2. 做 E2/E3：full vs anytime Wasserstein；
3. 做 E4：matching warm start；
4. 同时做 E6/E7：multiplicity + mandatory-flow；
5. 做 E8/E9：density/component adaptive kernel；
6. 再研究 E5 的 (q)-continuation；
7. 最后统一跑 (N\times)distribution(\times)cardinality-ratio(\times)batch-size；
8. 用全部结果反推 router。

这样即使最后证明“(W_p\to W_\infty) continuation 本身没什么用”，项目也不会失败，因为你仍然能得到一整套关于 **Wasserstein 与 Bottleneck 共享结构、适用区域和自适应算法选择** 的完整结果。

[1]: https://epubs.siam.org/doi/10.1137/1.9781611974317.9?utm_source=chatgpt.com "2016 Proceedings of the Meeting on Algorithm Engineering and Experiments (ALENEX) | Geometry Helps to Compare Persistence Diagrams"
[2]: https://arxiv.org/abs/2110.14734 "Approximating 1-Wasserstein Distance between Persistence Diagrams by Graph Sparsification"
[3]: https://gudhi.inria.fr/python/latest/bottleneck_distance_user.html "Bottleneck distance user manual — gudhi v3.13.0 documentation"
