我倾向于：**Large-N 的算法基础应该转向 Hera，但不要把整个内核“换成 Hera 源码”。更好的做法是把 Hera 的 geometric matching 思路做成你们 dispatcher 的一个 Large-N backend。**

这其实和你最早的路线一致：原计划就已经预判到 (n\ge 200) 后应该考虑 “KD-tree / Hera style”，而不是继续无限扩展 Small-N bitset。

原因很直接。你们现在最成功的路线本质上是在吃“小中规模问题可以塞进 cache / machine words”的红利：

[
\text{fixed bitset}
+\text{quickselect}
+\text{greedy Kuhn}
+\text{x-sweep}.
]

到了真正的 Large-N，例如

[
n=512,1024,2048,4096,\ldots,
]

继续显式或半显式维护越来越大的 adjacency，其根本问题开始变成 (n^2) 级候选、内存流量和图扫描，而不是 `popcount`、Kuhn 还是 HK 再快几个百分点。

Hera 恰好解决的是这一层问题。它的 bottleneck matching 采用几何 Hopcroft–Karp，通过 k-d tree 的近邻查询避免显式线性扫描整个候选邻接；论文的核心思想就是用 persistence diagram 的二维几何结构替代纯组合图扫描。 论文中几何方法在大规模实例上的 scaling 明显优于 combinatorial 方法，并测试到了每张图约 25,000 个点。

所以我现在会把架构改成这种，而不是简单的：

[
n<512:\text{ours},\qquad n\ge512:\text{Hera}.
]

更合理的是二维 dispatcher：

[
\boxed{
\operatorname{backend}
======================

f(n,m,\rho,\text{geometry})
}
]

例如最终可能变成：

[
\begin{cases}
\text{fixed-bitset + greedy Kuhn},&
N\text{ small},[2mm]
\text{x-sweep + bitset},&
N\text{ medium},[2mm]
\text{mandatory-flow},&
N\text{ large and }\rho\text{ very low},[2mm]
\text{Hera-style geometric HK},&
N\text{ large/general}.
\end{cases}
]

这里特别不能把 mandatory-flow 扔掉。你们已经观察到它在 256 点极稀疏分布上反超 expanded matching，因此在更大的极稀疏图上它反而可能继续胜过 Hera-style KD-tree。Large-N 也不是一个单一相区。

### 我更推荐“Hera-style”，而不是直接依赖 Hera

这是很重要的区别。

现代 Hera 仓库是 BSD 许可，可以作为 C++ header-only library 使用，但当前仓库还带 Boost、TBB 等依赖。([GitHub][1]) 如果你的目标最终是一个轻量级 Python/C++ bottleneck 库，我不太建议为了 Large-N 直接把整个 Hera dependency tree 拉进来。

更合理的是：

[
\boxed{
\text{读 Hera}
\rightarrow
\text{复现核心 geometric algorithm}
\rightarrow
\text{接入你们现有统一接口}
}
]

保留你们已经有的：

* `PreparedDiagram`
* exact threshold contract
* essential-point handling
* caller workspace
* one-to-many
* instrumentation
* dispatcher

只把 Large-N 的：

[
\text{adjacency/query + matching}
]

换成 geometric backend。

这样你才能继续做独立 ablation，而不是突然出现一个黑盒：

> `if n > 512: call hera()`

那会破坏你整个项目目前最有价值的“所有优化均可独立开关和归因”的实验结构。

### Hera 最值得借的其实有两层

第一层是 geometric matching。

普通 threshold matching 在半径 (r) 下需要不断找：

[
q\in V,\qquad
|p-q|_\infty\le r.
]

Small-N 可以直接 bitset：

[
N(p)=\text{bitmask}.
]

Large-N 则改成：

[
N_r(p)
\longrightarrow
\text{geometric range / near-neighbor query}.
]

Hera 使用 k-d tree，并通过 lazy deletion 和 subtree remaining-count 等技巧支持 Hopcroft–Karp 中的动态近邻查询。

第二层我认为更值得你借：**approximate first → exact refinement**。

Hera 的 exact 思路不是一开始就把所有 (O(n^2)) candidate distances 全部算出来。它先得到一个相对误差区间

[
a<d_B\le b,
]

然后 exact bottleneck 必须落在这个窄区间内，再只生成满足

[
a<|p-q|_\infty\le b
]

的 candidate pairs，最后对这些候选做精确搜索。

这对你们现在的架构非常契合。

你现在已经有：

[
\text{clipping}
+\text{quickselect}
+\text{x-sweep}.
]

Large-N 可以进一步变成：

[
\boxed{
\text{cheap geometric approximation}
\rightarrow
[a,b]
\rightarrow
\text{x/y annulus candidate extraction}
\rightarrow
\text{exact refinement}
}
]

这可能比“直接复制 Hera KD-tree matcher”还重要。

有一个语义问题必须注意：GUDHI 当前暴露的 `gudhi.hera.bottleneck_distance` 默认 `delta=0.01`，这是乘法近似，不是你们现在要求的 exact 语义。([GUDHI library][2]) 因此你们不能因为换 Large-N backend 就把默认语义偷偷变成 approximate。应当保持：

[
\boxed{\text{最终返回 exact}}
]

approximation 只作为内部缩小搜索区间的手段。

---

所以如果让我现在重画整个项目，我会定成三个阶段：

[
\boxed{\text{Small-N kernel}}
]

大约 (4\sim64/128)：

[
\text{fixed bitset + quickselect + greedy Kuhn}.
]

这一块你们现在已经非常成熟，而且 64 点已经有约 **19.93× reference**。

然后：

[
\boxed{\text{Medium-N kernel}}
]

大约 (128\sim512/1024)：

[
\text{x-sweep}
+\text{adaptive bitset}
+\text{mandatory-flow}.
]

你现在实际上正在完成这一层。

最后增加：

[
\boxed{\text{Large-N geometric kernel}}
]

大约 (512/1024+)：

[
\text{Hera-style geometric HK}
+
\text{approximation-assisted exact refinement}.
]

但这些数字只是 dispatcher calibration 的初始实验点，**不要硬编码 512**。真正的 crossover 应该通过：

[
(n,m,\rho)
]

benchmark 找出来。

我的判断因此是：

> **值得。甚至我认为这是下一阶段最自然的方向。**

你们已经把 Small/Medium-N combinatorial kernel 做得很深了。继续在 1024 点还死磕更大的 `uint64_t[]` bitset，边际价值会越来越低；Large-N 改研究 Hera 那套 geometric formulation，才是真正从“Small-N optimizer”迈向“通用 bottleneck solver”的关键一步。

[1]: https://github.com/anigmetov/hera "GitHub - anigmetov/hera · GitHub"
[2]: https://gudhi.inria.fr/python/latest/bottleneck_distance_user.html "Bottleneck distance user manual — gudhi v3.13.0 documentation"
