# 数学约定

[English](MATHEMATICS.en.md) · [仓库 README](https://github.com/proffitteoy/Topp/blob/main/README.md) · [API](API.md)

本文定义 Topp `v0.1.0` 公开 Python API 的距离语义。这里的 **exact** 表示算法不引入近似参数或近似容差；输入和返回值仍使用 IEEE 754 `float64`。

## Persistence diagram

一个 diagram 是带重数的点集，每个点写作 $x=(b,d)$，且 $b\leq d$。有限且严格位于对角线上方的点构成 $X_f$。重复行表示不同的点实例，其 multiplicity 会参与匹配。

有限对角点 $(a,a)$ 的持久度为零，会被距离计算忽略。对角线为

$$
\Delta=\{(t,t):t\in\mathbb{R}\}.
$$

在匹配中视为具有无限重数，因此两个 diagram 的有限点数可以不同。

## Bottleneck Distance

对有限部分，Topp 使用点间 $L^\infty$ 距离和对角线无限重数：

$$
d_B(X,Y)
=
\inf_{\gamma:X_f\cup\Delta\to Y_f\cup\Delta}
\sup_{x\in X_f\cup\Delta}
\lVert x-\gamma(x)\rVert_\infty.
$$

其中 $\gamma$ 是在 $X_f\cup\Delta$ 与 $Y_f\cup\Delta$ 之间保持重数的双射。有限点 $x=(b,d)$ 到对角线的代价为

$$
d_\infty(x,\Delta)
=
\inf_{z\in\Delta}\lVert x-z\rVert_\infty
=
\frac{d-b}{2}.
$$

`bottleneck_within(X, Y, t)` 精确判断上述距离是否满足 `d_B(X, Y) <= t`，不会先计算近似距离再比较。

## Wasserstein Distance

对于 $1\leq p,q<\infty$，Topp 使用以下约定：

$$
W_{q,p}(X,Y)
=
\left(
\inf_{\gamma:X_f\cup\Delta\to Y_f\cup\Delta}
\sum_{x\in X_f\cup\Delta}
\lVert x-\gamma(x)\rVert_p^q
\right)^{1/q}.
$$

其中 `order` 对应 $q$，`internal_p` 对应 $p$，匹配同样包含无限重数的对角线。`v0.1.0` 仅实现：

- `order=1, internal_p=np.inf`，即 $W_{1,\infty}$；
- `order=2, internal_p=2`，即 $W_{2,2}$。

对于有限点 $x=(b,d)$，两个已支持度量的对角线代价为：

$$
d_\infty(x,\Delta)=\frac{d-b}{2},
\qquad
d_2(x,\Delta)=\frac{d-b}{\sqrt{2}}.
$$

因此 $W_{2,2}$ 中该对角线匹配对平方和的贡献为 $(d-b)^2/2$。其他 `(order, internal_p)` 组合会抛出 `NotImplementedError`，不会被映射到最接近的已支持度量。

## Essential points

支持三种 essential point：

| 类型 | 形式 | 同类型点间代价 |
|---|---|---|
| 正向 essential | $(b,+\infty)$ | $|b_1-b_2|$ |
| 负向 essential | $(-\infty,d)$ | $|d_1-d_2|$ |
| 完全 essential | $(-\infty,+\infty)$ | $0$ |

Essential point 只能与相同类型匹配，不能匹配到对角线。任一类型的 multiplicity 在两个 diagram 中不相等时，Bottleneck 和 Wasserstein 均返回 `inf`。

当 multiplicity 相等时，同类型的有限坐标按排序后的顺序配对。若得到 essential 代价 $e_1,\ldots,e_k$，则 Bottleneck 使用 $\max_i e_i$，$W_{1,\infty}$ 使用 $\sum_i e_i$，$W_{2,2}$ 将 $\sum_i e_i^2$ 计入最终开平方前的总和。完全 essential points 只要求两侧数量相等。

## 边界行为

- 空 diagram 与空 diagram 的距离为 `0`；
- 重复的非对角点保留 multiplicity，不会按集合去重；
- 对角点可出现在输入和 `PreparedDiagram.n_points` 中，但不计入 `n_finite_points`；
- NaN、`birth > death`、`birth=+inf`、`death=-inf` 及非法无穷组合会在 Python 层抛出 `ValueError`；
- 距离可能为 `inf`，但合法输入不会被静默交换、删除或修正。

具体调用形式和异常类型见 [API 文档](API.md)。
