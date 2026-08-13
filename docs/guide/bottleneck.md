# Bottleneck 距离

[English](https://proffitteoy.github.io/Topp/en/guide/bottleneck.html)

`bottleneck_distance(X, Y)` 使用 $L^\infty$ 点代价计算精确的 Bottleneck 距离。

## 成对调用

```python
import numpy as np
import topp

x = np.array([[0.0, 1.0], [0.2, 0.7]])
y = np.array([[0.0, 1.1], [0.3, 0.8]])

distance = topp.bottleneck_distance(x, y)
print(distance)
```

结果是 Python `float`。空 diagram、对角点、重复点和受支持的 essential points 都是合法输入。

## 匹配到对角线

有限点 $x=(b,d)$ 匹配到对角线的代价为

$$
d_\infty(x,\Delta)=\frac{d-b}{2}.
$$

例如：

```pycon
>>> topp.bottleneck_distance([[0.0, 2.0]], [])
1.0
```

## 对称性与预处理输入

```python
prepared = topp.prepare_diagram(x)

forward = topp.bottleneck_distance(prepared, y)
reverse = topp.bottleneck_distance(y, prepared)
assert forward == reverse
```

预处理不会改变数学结果；重复使用同一个 diagram 时，它可以避免重复执行输入预处理。

完整定义见[数学约定](../MATHEMATICS.md#bottleneck-distance)。
