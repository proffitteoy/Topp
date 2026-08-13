# 线程与 GIL

[English](https://proffitteoy.github.io/Topp/en/guide/threads.html)

原生距离计算会释放 Python GIL。`PreparedDiagram` 不可变，可以被多个线程并发读取。

```python
from concurrent.futures import ThreadPoolExecutor

import topp

query = topp.prepare_diagram([[0.0, 1.0], [0.25, 0.75]])
targets = [
    [[0.0, 1.1], [0.30, 0.80]],
    [[0.0, 2.0]],
    [],
]

with ThreadPoolExecutor(max_workers=3) as pool:
    distances = list(
        pool.map(lambda target: topp.bottleneck_distance(query, target), targets)
    )
```

不要在没有应用层同步的情况下并发写入同一个 `out` 数组。不同的输出数组相互独立。

只有周边工作负载和 diagram 规模足以抵消调度开销时，线程才有帮助。请对真实应用进行 benchmark，不要假设更多 worker 一定更快。
