#pragma once

#include <bottleneck/core.hpp>

#include <cstdint>
#include <memory>

namespace bottleneck {

// The first optimized kernels intentionally cover the two exact metrics from
// Phase 2 of docs/research/proposals/wasserstein优化方案.md. More metrics can be added without
// changing the distance entry points.
enum class WassersteinMetric {
  w1_linf,
  w2_l2,
};

enum class WassersteinCandidateStrategy {
  dense_scalar,
  dense_blocked,
  dense_avx2,
  dense_parallel,
  sweep_binary,
  sweep_two_pointer,
  topk_pricing_full_scan,
  topk_pricing_simd,
  topk_pricing_sweep,
  topk_pricing_kdtree,
  topk_pricing_kdtree_persistent,
  topk_pricing_adaptive,
  topk_pricing_sweep_incremental,
  topk_pricing_sweep_persistent,
  topk_pricing_sweep_dynamic,
  topk_pricing_sweep_dynamic_batched,
  adaptive,
};

enum class WassersteinGraphStrategy {
  dense_matrix,
  csr,
  row_vectors,
  fixed_degree_4,
  fixed_degree_8,
  fixed_degree_16,
  fixed_degree_32,
  bitmask_lazy,
  block_sparse_16,
  block_sparse_32,
  adaptive,
};

enum class WassersteinMatcherStrategy {
  dense_hungarian,
  dense_sap,
  dense_sap_row_reduction,
  dense_sap_jv_reduction,
  sparse_sap,
  sparse_sap_arena,
  adaptive,
};

enum class WassersteinComponentStrategy {
  none,
  dense,
  sparse,
  tiny_sparse,
  parallel_dense,
  parallel_sparse,
  adaptive,
};

enum class WassersteinWarmStart {
  none,
  row_max,
  global_descending,
};

enum class WassersteinDuplicateStrategy {
  none,
  exact,
  adaptive,
};

struct WassersteinConfig {
  WassersteinMetric metric = WassersteinMetric::w1_linf;
  WassersteinCandidateStrategy candidates = WassersteinCandidateStrategy::adaptive;
  WassersteinGraphStrategy graph = WassersteinGraphStrategy::adaptive;
  WassersteinMatcherStrategy matcher = WassersteinMatcherStrategy::adaptive;
  WassersteinComponentStrategy components = WassersteinComponentStrategy::adaptive;
  WassersteinWarmStart warm_start = WassersteinWarmStart::global_descending;
  WassersteinDuplicateStrategy duplicates = WassersteinDuplicateStrategy::adaptive;
  std::size_t top_k = 8;
};

struct WassersteinStats {
  std::uint64_t possible_pairs = 0;
  std::uint64_t candidate_pairs = 0;
  std::uint64_t positive_edges = 0;
  std::uint64_t pruned_pairs = 0;
  std::uint64_t active_rows = 0;
  std::uint64_t active_columns = 0;
  std::uint64_t augmentations = 0;
  std::uint64_t pricing_rounds = 0;
  std::uint64_t pricing_full_scan_rounds = 0;
  std::uint64_t pricing_simd_rounds = 0;
  std::uint64_t pricing_sweep_rounds = 0;
  std::uint64_t pricing_kdtree_rounds = 0;
  std::uint64_t pricing_kdtree_builds = 0;
  std::uint64_t pricing_kdtree_updates = 0;
  std::uint64_t priced_edges = 0;
  std::uint64_t pricing_violations = 0;
  std::uint64_t peak_materialized_edges = 0;
  std::uint64_t max_degree = 0;
  std::uint64_t sparse_fallbacks = 0;
  std::uint64_t sparse_arena_builds = 0;
  std::uint64_t sparse_scratch_reuses = 0;
  std::uint64_t sparse_heap_growths = 0;
  std::uint64_t peak_sparse_arena_bytes = 0;
  std::uint64_t dynamic_inserted_edges = 0;
  std::uint64_t dynamic_dijkstra_runs = 0;
  std::uint64_t dynamic_batch_groups = 0;
  std::uint64_t jv_fallbacks = 0;
  std::uint64_t component_count = 0;
  std::uint64_t largest_component = 0;
  std::uint64_t tiny_components = 0;
  std::uint64_t greedy_matches = 0;
  std::uint64_t warm_start_certificates = 0;
  std::uint64_t duplicate_groups = 0;
  std::uint64_t duplicate_points_removed = 0;
  std::uint64_t graph_bytes = 0;
  std::uint64_t peak_graph_bytes = 0;
  std::uint64_t prepare_time_ns = 0;
  std::uint64_t candidate_time_ns = 0;
  std::uint64_t graph_time_ns = 0;
  std::uint64_t component_time_ns = 0;
  std::uint64_t solver_time_ns = 0;
  std::uint64_t pricing_time_ns = 0;
};

class WassersteinWorkspace {
 public:
  struct Impl;

  WassersteinWorkspace();
  ~WassersteinWorkspace();
  WassersteinWorkspace(WassersteinWorkspace&&) noexcept;
  WassersteinWorkspace& operator=(WassersteinWorkspace&&) noexcept;
  WassersteinWorkspace(const WassersteinWorkspace&) = delete;
  WassersteinWorkspace& operator=(const WassersteinWorkspace&) = delete;

 private:
  std::unique_ptr<Impl> impl_;

  friend double wasserstein_distance(
      const PreparedDiagram&, const PreparedDiagram&, WassersteinWorkspace&,
      const WassersteinConfig&, WassersteinStats*);
  friend void wasserstein_distances(
      const PreparedDiagram&, std::span<const PreparedDiagram>, std::span<double>,
      WassersteinWorkspace&, const WassersteinConfig&, WassersteinStats*);
};

[[nodiscard]] double wasserstein_distance(
    const PreparedDiagram& first,
    const PreparedDiagram& second,
    const WassersteinConfig& config = {},
    WassersteinStats* stats = nullptr);

[[nodiscard]] double wasserstein_distance(
    const PreparedDiagram& first,
    const PreparedDiagram& second,
    WassersteinWorkspace& workspace,
    const WassersteinConfig& config = {},
    WassersteinStats* stats = nullptr);

[[nodiscard]] double wasserstein_distance(
    const Diagram& first,
    const Diagram& second,
    const WassersteinConfig& config = {},
    WassersteinStats* stats = nullptr);

void wasserstein_distances(
    const PreparedDiagram& query,
    std::span<const PreparedDiagram> diagrams,
    std::span<double> output,
    const WassersteinConfig& config = {},
    WassersteinStats* stats = nullptr);

void wasserstein_distances(
    const PreparedDiagram& query,
    std::span<const PreparedDiagram> diagrams,
    std::span<double> output,
    WassersteinWorkspace& workspace,
    const WassersteinConfig& config = {},
    WassersteinStats* stats = nullptr);

[[nodiscard]] std::vector<double> wasserstein_distances(
    const PreparedDiagram& query,
    std::span<const PreparedDiagram> diagrams,
    const WassersteinConfig& config = {},
    WassersteinStats* stats = nullptr);

[[nodiscard]] const char* to_string(WassersteinMetric metric) noexcept;
[[nodiscard]] const char* to_string(WassersteinCandidateStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinGraphStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinMatcherStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinComponentStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinWarmStart strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinDuplicateStrategy strategy) noexcept;

}  // namespace bottleneck
