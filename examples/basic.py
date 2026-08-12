import numpy as np

import topp


query = np.array([[0.0, 1.0], [0.3, 0.8]])
target = np.array([[0.0, 1.1], [0.4, 0.9]])

prepared = topp.prepare_diagram(query)

print("bottleneck:", topp.bottleneck_distance(prepared, target))
print(
    "wasserstein W2-L2:",
    topp.wasserstein_distance(prepared, target, order=2, internal_p=2),
)
print(
    "batch:",
    topp.bottleneck_distances(prepared, [target, np.empty((0, 2))]),
)
print("within 0.1:", topp.bottleneck_within(prepared, target, 0.1))
