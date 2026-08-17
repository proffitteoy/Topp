#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace bottleneck {

struct Point {
  double birth;
  double death;
};

using Diagram = std::vector<Point>;

enum class CandidateStrategy {
  sort_all,
  sort_unique,
  sort_unique_clipped,
  sort_unique_greedy_clipped,
  x_sweep_clipped,
};

enum class ThresholdStrategy {
  binary,
  gudhi_alpha,
  exponential,
  quickselect,
  incremental,
  incremental_blocked,
  geometric_refinement,
  adaptive,
};

enum class DistanceStrategy {
  dense_aos,
  dense_soa,
  dense_soa_avx2,
  recompute_aos,
  recompute_soa,
};

enum class AdjacencyStrategy {
  on_demand,
  dense_byte,
  sparse_csr,
  x_sweep_csr,
  bitset64,
  bitset128,
  bitset256,
  bitset_auto,
  adaptive,
  dynamic_bitset,
};

enum class MatcherStrategy {
  kuhn,
  greedy_kuhn,
  reusable_greedy_kuhn,
  fixed_greedy_kuhn,
  constraint_kuhn,
  component_kuhn,
  mandatory_flow,
  multiplicity_flow,
  adaptive,
  hopcroft_karp,
  greedy_hopcroft_karp,
  geometric_hopcroft_karp,
  mandatory_sparse_flow,
};

enum class VertexOrder {
  natural,
  degree_ascending,
};

struct SolverConfig {
  CandidateStrategy candidates = CandidateStrategy::sort_unique_clipped;
  ThresholdStrategy threshold = ThresholdStrategy::adaptive;
  DistanceStrategy distance = DistanceStrategy::dense_aos;
  AdjacencyStrategy adjacency = AdjacencyStrategy::adaptive;
  MatcherStrategy matcher = MatcherStrategy::adaptive;
  VertexOrder vertex_order = VertexOrder::degree_ascending;
  double lower_bound_hint = 0.0;
};

struct SolverStats {
  std::uint64_t raw_candidates = 0;
  std::uint64_t retained_candidates = 0;
  std::uint64_t clipped_candidates = 0;
  std::uint64_t threshold_decisions = 0;
  std::uint64_t adjacency_checks = 0;
  std::uint64_t emitted_edges = 0;
  std::uint64_t augment_searches = 0;
  std::uint64_t bfs_phases = 0;
  std::uint64_t constraint_rounds = 0;
  std::uint64_t forced_matches = 0;
  std::uint64_t component_count = 0;
  std::uint64_t component_rejects = 0;
  std::uint64_t largest_component_vertices = 0;
  std::uint64_t lower_bound_rejects = 0;
  std::uint64_t geometric_queries = 0;
  std::uint64_t kd_nodes_visited = 0;
  std::uint64_t refinement_rounds = 0;
  std::uint64_t incremental_groups_added = 0;
  std::uint64_t incremental_rollbacks = 0;
  std::uint64_t multiplicity_groups = 0;
  std::uint64_t multiplicity_points_removed = 0;
  std::uint64_t capacity_edges = 0;
  std::uint64_t x_window_candidates = 0;
  std::uint64_t router_x_window_pairs = 0;
  std::uint64_t mandatory_vertices = 0;
  std::uint64_t optional_pairs_pruned = 0;
  std::uint64_t router_identity_shortcuts = 0;
  std::uint64_t router_no_cross_shortcuts = 0;
  std::uint64_t router_multiplicity_routes = 0;
  std::uint64_t router_mandatory_routes = 0;
  std::uint64_t router_refinement_routes = 0;
  std::uint64_t router_quickselect_routes = 0;
  std::uint64_t lower_bound_hints_tested = 0;
  std::uint64_t lower_bound_hints_accepted = 0;
};

class PreparedDiagram {
 public:
  explicit PreparedDiagram(const Diagram& diagram);

  [[nodiscard]] const Diagram& finite_points() const noexcept;
  [[nodiscard]] const std::vector<double>& finite_births() const noexcept;
  [[nodiscard]] const std::vector<double>& finite_deaths() const noexcept;
  [[nodiscard]] const std::vector<double>& finite_diagonal_distances() const noexcept;
  [[nodiscard]] const std::vector<double>& finite_midpoints() const noexcept;
  [[nodiscard]] const std::vector<double>& finite_half_persistences() const noexcept;
  [[nodiscard]] const std::vector<double>& sorted_finite_midpoints() const noexcept;
  [[nodiscard]] const std::vector<std::size_t>& finite_midpoint_order() const noexcept;
  [[nodiscard]] double max_finite_half_persistence() const noexcept;
  [[nodiscard]] const std::vector<double>& sorted_finite_births() const noexcept;
  [[nodiscard]] const std::vector<std::size_t>& finite_birth_order() const noexcept;
  [[nodiscard]] const std::vector<std::size_t>&
  finite_duplicate_representatives() const noexcept;
  [[nodiscard]] const std::vector<std::size_t>&
  finite_duplicate_multiplicities() const noexcept;
  [[nodiscard]] double max_finite_diagonal_distance() const noexcept;
  [[nodiscard]] const std::vector<double>& positive_infinite_births() const noexcept;
  [[nodiscard]] const std::vector<double>& negative_infinite_deaths() const noexcept;
  [[nodiscard]] std::size_t fully_infinite_count() const noexcept;

 private:
  Diagram finite_points_;
  std::vector<double> finite_births_;
  std::vector<double> finite_deaths_;
  std::vector<double> finite_diagonal_distances_;
  std::vector<double> finite_midpoints_;
  std::vector<double> finite_half_persistences_;
  std::vector<double> sorted_finite_midpoints_;
  std::vector<std::size_t> finite_midpoint_order_;
  double max_finite_half_persistence_ = 0.0;
  std::vector<double> sorted_finite_births_;
  std::vector<std::size_t> finite_birth_order_;
  std::vector<std::size_t> finite_duplicate_representatives_;
  std::vector<std::size_t> finite_duplicate_multiplicities_;
  double max_finite_diagonal_distance_ = 0.0;
  std::vector<double> positive_infinite_births_;
  std::vector<double> negative_infinite_deaths_;
  std::size_t fully_infinite_count_ = 0;
};

[[nodiscard]] double bottleneck_distance(
    const PreparedDiagram& first,
    const PreparedDiagram& second,
    const SolverConfig& config = {},
    SolverStats* stats = nullptr);

[[nodiscard]] double bottleneck_distance(
    const Diagram& first,
    const Diagram& second,
    const SolverConfig& config = {},
    SolverStats* stats = nullptr);

[[nodiscard]] bool bottleneck_within(
    const PreparedDiagram& first,
    const PreparedDiagram& second,
    double threshold,
    const SolverConfig& config = {},
    SolverStats* stats = nullptr);

void bottleneck_distances(
    const PreparedDiagram& query,
    std::span<const PreparedDiagram> diagrams,
    std::span<double> output,
    const SolverConfig& config = {},
    SolverStats* stats = nullptr);

[[nodiscard]] std::vector<double> bottleneck_distances(
    const PreparedDiagram& query,
    std::span<const PreparedDiagram> diagrams,
    const SolverConfig& config = {},
    SolverStats* stats = nullptr);

[[nodiscard]] const char* to_string(CandidateStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(ThresholdStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(DistanceStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(AdjacencyStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(MatcherStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(VertexOrder order) noexcept;

}  // namespace bottleneck
