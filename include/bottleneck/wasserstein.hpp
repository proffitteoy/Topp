#pragma once

#include <bottleneck/core.hpp>

#include <cstdint>

namespace bottleneck {

// The first optimized kernels intentionally cover the two exact metrics from
// Phase 2 of docs/wasserstein优化方案.md. More metrics can be added without
// changing the distance entry points.
enum class WassersteinMetric {
  w1_linf,
  w2_l2,
};

enum class WassersteinCandidateStrategy {
  dense_scalar,
  dense_blocked,
  dense_avx2,
  sweep_binary,
  sweep_two_pointer,
  adaptive,
};

enum class WassersteinGraphStrategy {
  dense_matrix,
  csr,
  row_vectors,
  adaptive,
};

enum class WassersteinMatcherStrategy {
  dense_hungarian,
  dense_sap,
  sparse_sap,
  adaptive,
};

enum class WassersteinComponentStrategy {
  none,
  dense,
  sparse,
  tiny_sparse,
  adaptive,
};

enum class WassersteinWarmStart {
  none,
  row_max,
  global_descending,
};

struct WassersteinConfig {
  WassersteinMetric metric = WassersteinMetric::w1_linf;
  WassersteinCandidateStrategy candidates = WassersteinCandidateStrategy::dense_scalar;
  WassersteinGraphStrategy graph = WassersteinGraphStrategy::dense_matrix;
  WassersteinMatcherStrategy matcher = WassersteinMatcherStrategy::dense_hungarian;
  WassersteinComponentStrategy components = WassersteinComponentStrategy::none;
  WassersteinWarmStart warm_start = WassersteinWarmStart::none;
};

struct WassersteinStats {
  std::uint64_t possible_pairs = 0;
  std::uint64_t candidate_pairs = 0;
  std::uint64_t positive_edges = 0;
  std::uint64_t pruned_pairs = 0;
  std::uint64_t active_rows = 0;
  std::uint64_t active_columns = 0;
  std::uint64_t augmentations = 0;
  std::uint64_t component_count = 0;
  std::uint64_t largest_component = 0;
  std::uint64_t tiny_components = 0;
  std::uint64_t greedy_matches = 0;
  std::uint64_t warm_start_certificates = 0;
  std::uint64_t graph_bytes = 0;
  std::uint64_t prepare_time_ns = 0;
  std::uint64_t candidate_time_ns = 0;
  std::uint64_t graph_time_ns = 0;
  std::uint64_t component_time_ns = 0;
  std::uint64_t solver_time_ns = 0;
};

[[nodiscard]] double wasserstein_distance(
    const PreparedDiagram& first,
    const PreparedDiagram& second,
    const WassersteinConfig& config = {},
    WassersteinStats* stats = nullptr);

[[nodiscard]] double wasserstein_distance(
    const Diagram& first,
    const Diagram& second,
    const WassersteinConfig& config = {},
    WassersteinStats* stats = nullptr);

[[nodiscard]] const char* to_string(WassersteinMetric metric) noexcept;
[[nodiscard]] const char* to_string(WassersteinCandidateStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinGraphStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinMatcherStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinComponentStrategy strategy) noexcept;
[[nodiscard]] const char* to_string(WassersteinWarmStart strategy) noexcept;

}  // namespace bottleneck
