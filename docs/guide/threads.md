# Threads and the GIL

Native distance calculations release the Python GIL. A `PreparedDiagram` is immutable and may be read concurrently.

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

Do not concurrently write to the same `out` array without application-level synchronization. Separate output arrays are independent.

Threading is useful only when the surrounding workload and diagram sizes justify its scheduling overhead. Benchmark the real application rather than assuming that more workers are faster.
