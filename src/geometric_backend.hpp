#pragma once

#include <bottleneck/core.hpp>

namespace bottleneck::detail {

[[nodiscard]] bool geometric_bottleneck_within(
    const PreparedDiagram& first,
    const PreparedDiagram& second,
    double threshold,
    SolverStats* stats);

[[nodiscard]] double geometric_refined_bottleneck_distance(
    const PreparedDiagram& first,
    const PreparedDiagram& second,
    SolverStats* stats);

}  // namespace bottleneck::detail
