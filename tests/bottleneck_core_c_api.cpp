#include <bottleneck/core.hpp>

#include <cstddef>
#include <limits>

#if defined(_WIN32)
#define BOTTLENECK_EXPORT __declspec(dllexport)
#else
#define BOTTLENECK_EXPORT __attribute__((visibility("default")))
#endif

namespace {

bottleneck::Diagram make_diagram(const double* values, std::size_t size) {
  bottleneck::Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    result.push_back({values[index * 2], values[index * 2 + 1]});
  }
  return result;
}

bottleneck::SolverConfig make_config(int candidates, int threshold, int distance, int adjacency,
                                     int matcher, int order) {
  return {
      candidates == 0 ? bottleneck::CandidateStrategy::sort_all
                      : candidates == 1 ? bottleneck::CandidateStrategy::sort_unique
                      : candidates == 2 ? bottleneck::CandidateStrategy::sort_unique_clipped
                                        : bottleneck::CandidateStrategy::sort_unique_greedy_clipped,
      threshold == 0 ? bottleneck::ThresholdStrategy::binary
                     : threshold == 1 ? bottleneck::ThresholdStrategy::gudhi_alpha
                     : threshold == 2 ? bottleneck::ThresholdStrategy::exponential
                     : threshold == 3 ? bottleneck::ThresholdStrategy::incremental
                     : threshold == 4 ? bottleneck::ThresholdStrategy::incremental_blocked
                     : threshold == 6 ? bottleneck::ThresholdStrategy::geometric_refinement
                     : threshold == 7 ? bottleneck::ThresholdStrategy::adaptive
                                      : bottleneck::ThresholdStrategy::quickselect,
      distance == 0 ? bottleneck::DistanceStrategy::dense_aos
                    : distance == 1 ? bottleneck::DistanceStrategy::dense_soa
                    : distance == 2 ? bottleneck::DistanceStrategy::recompute_aos
                    : distance == 3 ? bottleneck::DistanceStrategy::recompute_soa
                                    : bottleneck::DistanceStrategy::dense_soa_avx2,
      adjacency == 0 ? bottleneck::AdjacencyStrategy::on_demand
                     : adjacency == 1 ? bottleneck::AdjacencyStrategy::dense_byte
                     : adjacency == 2 ? bottleneck::AdjacencyStrategy::sparse_csr
                     : adjacency == 3 ? bottleneck::AdjacencyStrategy::bitset64
                     : adjacency == 4 ? bottleneck::AdjacencyStrategy::dynamic_bitset
                     : adjacency == 5 ? bottleneck::AdjacencyStrategy::bitset128
                     : adjacency == 6 ? bottleneck::AdjacencyStrategy::bitset256
                     : adjacency == 7 ? bottleneck::AdjacencyStrategy::bitset_auto
                     : adjacency == 8 ? bottleneck::AdjacencyStrategy::x_sweep_csr
                                      : bottleneck::AdjacencyStrategy::adaptive,
      matcher == 0 ? bottleneck::MatcherStrategy::kuhn
                   : matcher == 1 ? bottleneck::MatcherStrategy::greedy_kuhn
                   : matcher == 2 ? bottleneck::MatcherStrategy::constraint_kuhn
                   : matcher == 3 ? bottleneck::MatcherStrategy::hopcroft_karp
                   : matcher == 4 ? bottleneck::MatcherStrategy::greedy_hopcroft_karp
                   : matcher == 5 ? bottleneck::MatcherStrategy::component_kuhn
                   : matcher == 6 ? bottleneck::MatcherStrategy::reusable_greedy_kuhn
                   : matcher == 7 ? bottleneck::MatcherStrategy::mandatory_flow
                   : matcher == 8 ? bottleneck::MatcherStrategy::fixed_greedy_kuhn
                   : matcher == 10 ? bottleneck::MatcherStrategy::geometric_hopcroft_karp
                                  : bottleneck::MatcherStrategy::adaptive,
      order == 0 ? bottleneck::VertexOrder::natural : bottleneck::VertexOrder::degree_ascending,
  };
}

}  // namespace

extern "C" BOTTLENECK_EXPORT double bottleneck_core_distance(
    const double* first, std::size_t first_size, const double* second, std::size_t second_size,
    int candidates, int threshold, int distance, int adjacency, int matcher, int order) noexcept {
  try {
    return bottleneck::bottleneck_distance(make_diagram(first, first_size),
                                           make_diagram(second, second_size),
                                           make_config(candidates, threshold, distance, adjacency, matcher, order));
  } catch (...) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}
