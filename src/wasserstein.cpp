#include <bottleneck/wasserstein.hpp>

#include <algorithm>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#if defined(BOTTLENECK_HAVE_WASSERSTEIN_AVX2) && defined(_MSC_VER)
#include <intrin.h>
#endif

#if defined(BOTTLENECK_HAVE_WASSERSTEIN_AVX2)
namespace bottleneck::detail {
void fill_w1_savings_avx2(double first_u, double first_v, const double* second_u,
                          const double* second_v, std::size_t size,
                          double* output) noexcept;
void fill_w2_savings_avx2(double first_u, double first_v, const double* second_u,
                          const double* second_v, std::size_t size,
                          double* output) noexcept;
}  // namespace bottleneck::detail
#endif

namespace bottleneck {

namespace {

using Clock = std::chrono::steady_clock;
using Weight = long double;

constexpr std::size_t dense_block_size = 32;
constexpr std::size_t tiny_component_limit = 8;

std::uint64_t elapsed_ns(Clock::time_point start, Clock::time_point stop) {
  return static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count());
}

[[maybe_unused]] bool cpu_has_avx2() noexcept {
#if defined(BOTTLENECK_HAVE_WASSERSTEIN_AVX2) && defined(_MSC_VER) && defined(_M_X64)
  static const bool available = [] {
    int registers[4]{};
    __cpuid(registers, 1);
    constexpr int osxsave = 1 << 27;
    constexpr int avx = 1 << 28;
    if ((registers[2] & (osxsave | avx)) != (osxsave | avx) ||
        (_xgetbv(0) & 0x6) != 0x6) {
      return false;
    }
    __cpuidex(registers, 7, 0);
    return (registers[1] & (1 << 5)) != 0;
  }();
  return available;
#else
  return false;
#endif
}

struct Edge {
  std::uint32_t column = 0;
  double saving = 0.0;
};

Weight raw_saving(const Point& first, const Point& second,
                  double first_half_persistence,
                  double second_half_persistence,
                  WassersteinMetric metric);

struct CandidateRows {
  CandidateRows(std::size_t row_count, std::size_t column_count)
      : rows(row_count), columns(column_count), edges(row_count) {}

  void reset(std::size_t row_count, std::size_t column_count) {
    rows = row_count;
    columns = column_count;
    edges.resize(row_count);
    for (auto& row : edges) {
      row.clear();
    }
    edge_count = 0;
  }

  void add(std::size_t row, std::size_t column, Weight saving,
           WassersteinStats* stats) {
    if (saving <= Weight{0}) {
      return;
    }
    edges[row].push_back(
        {static_cast<std::uint32_t>(column), static_cast<double>(saving)});
    ++edge_count;
    if (stats != nullptr) {
      ++stats->positive_edges;
    }
  }

  std::size_t rows;
  std::size_t columns;
  std::vector<std::vector<Edge>> edges;
  std::size_t edge_count = 0;
};

enum class StoredGraphKind {
  dense,
  csr,
  rows,
  fixed_degree,
  bitmask_lazy,
  block_sparse,
};

struct WeightedGraph {
  void reset() {
    kind = StoredGraphKind::rows;
    rows = 0;
    columns = 0;
    edge_count = 0;
    dense.clear();
    offsets.clear();
    column_indices.clear();
    weights.clear();
    row_edges.clear();
    fixed_capacity = 0;
    fixed_degrees.clear();
    fixed_columns.clear();
    fixed_weights.clear();
    lazy_first = nullptr;
    lazy_second = nullptr;
    mask_words = 0;
    positive_masks.clear();
    block_size = 0;
    block_row_offsets.clear();
    block_columns.clear();
    block_weights.clear();
  }

  template <class Function>
  void for_each_edge(std::size_t row, Function&& function) const {
    if (kind == StoredGraphKind::dense) {
      const std::size_t offset = row * columns;
      for (std::size_t column = 0; column < columns; ++column) {
        const double saving = dense[offset + column];
        if (saving > 0.0) {
          function(column, saving);
        }
      }
      return;
    }
    if (kind == StoredGraphKind::csr) {
      for (std::size_t index = offsets[row]; index < offsets[row + 1]; ++index) {
        function(static_cast<std::size_t>(column_indices[index]), weights[index]);
      }
      return;
    }
    if (kind == StoredGraphKind::fixed_degree) {
      const std::uint8_t degree = fixed_degrees[row];
      if (degree != fixed_degree_overflow) {
        const std::size_t base = row * fixed_capacity;
        for (std::size_t index = 0; index < degree; ++index) {
          function(static_cast<std::size_t>(fixed_columns[base + index]),
                   fixed_weights[base + index]);
        }
      } else {
        for (std::size_t index = offsets[row]; index < offsets[row + 1]; ++index) {
          function(static_cast<std::size_t>(column_indices[index]), weights[index]);
        }
      }
      return;
    }
    if (kind == StoredGraphKind::bitmask_lazy) {
      const auto& first_points = lazy_first->finite_points();
      const auto& second_points = lazy_second->finite_points();
      const auto& first_v = lazy_first->finite_half_persistences();
      const auto& second_v = lazy_second->finite_half_persistences();
      const std::size_t mask_base = row * mask_words;
      for (std::size_t word = 0; word < mask_words; ++word) {
        std::uint64_t remaining = positive_masks[mask_base + word];
        while (remaining != 0) {
          const unsigned bit = std::countr_zero(remaining);
          const std::size_t column = word * 64 + bit;
          const Weight saving = raw_saving(first_points[row], second_points[column],
                                           first_v[row], second_v[column], lazy_metric);
          if (saving > Weight{0}) {
            function(column, static_cast<double>(saving));
          }
          remaining &= remaining - 1;
        }
      }
      return;
    }
    if (kind == StoredGraphKind::block_sparse) {
      const std::size_t block_row = row / block_size;
      const std::size_t local_row = row % block_size;
      for (std::size_t tile = block_row_offsets[block_row];
           tile < block_row_offsets[block_row + 1]; ++tile) {
        const std::size_t column_base =
            static_cast<std::size_t>(block_columns[tile]) * block_size;
        const std::size_t column_count =
            (std::min)(block_size, columns - column_base);
        const std::size_t weight_base =
            tile * block_size * block_size + local_row * block_size;
        for (std::size_t local_column = 0; local_column < column_count;
             ++local_column) {
          const double saving = block_weights[weight_base + local_column];
          if (saving > 0.0) {
            function(column_base + local_column, saving);
          }
        }
      }
      return;
    }
    for (const Edge& edge : row_edges[row]) {
      function(static_cast<std::size_t>(edge.column), edge.saving);
    }
  }

  StoredGraphKind kind = StoredGraphKind::rows;
  std::size_t rows = 0;
  std::size_t columns = 0;
  std::size_t edge_count = 0;
  std::vector<double> dense;
  std::vector<std::size_t> offsets;
  std::vector<std::uint32_t> column_indices;
  std::vector<double> weights;
  std::vector<std::vector<Edge>> row_edges;
  static constexpr std::uint8_t fixed_degree_overflow = 0xff;
  std::size_t fixed_capacity = 0;
  std::vector<std::uint8_t> fixed_degrees;
  std::vector<std::uint32_t> fixed_columns;
  std::vector<double> fixed_weights;
  const PreparedDiagram* lazy_first = nullptr;
  const PreparedDiagram* lazy_second = nullptr;
  WassersteinMetric lazy_metric = WassersteinMetric::w1_linf;
  std::size_t mask_words = 0;
  std::vector<std::uint64_t> positive_masks;
  std::size_t block_size = 0;
  std::vector<std::size_t> block_row_offsets;
  std::vector<std::uint32_t> block_columns;
  std::vector<double> block_weights;
};

}  // namespace

struct WassersteinWorkspace::Impl {
  CandidateRows candidates{0, 0};
  WeightedGraph graph;
  std::vector<long double> short_potential;
  std::vector<long double> long_potential;
  std::vector<long double> minimum;
  std::vector<std::size_t> matched_short;
  std::vector<std::size_t> predecessor;
  std::vector<std::uint8_t> used;
  std::vector<std::uint8_t> initially_matched;
};

namespace {

Weight diagonal_power(double half_persistence, WassersteinMetric metric) {
  const Weight persistence = static_cast<Weight>(half_persistence);
  return metric == WassersteinMetric::w1_linf
             ? persistence
             : Weight{2} * persistence * persistence;
}

Weight cross_power(const Point& first, const Point& second,
                   WassersteinMetric metric) {
  const Weight birth_delta =
      std::fabs(static_cast<Weight>(first.birth) - second.birth);
  const Weight death_delta =
      std::fabs(static_cast<Weight>(first.death) - second.death);
  return metric == WassersteinMetric::w1_linf
             ? (std::max)(birth_delta, death_delta)
             : birth_delta * birth_delta + death_delta * death_delta;
}

Weight raw_saving(const Point& first, const Point& second,
                  double first_half_persistence, double second_half_persistence,
                  WassersteinMetric metric) {
  const Weight first_diagonal = diagonal_power(first_half_persistence, metric);
  const Weight second_diagonal = diagonal_power(second_half_persistence, metric);
  return first_diagonal + second_diagonal - cross_power(first, second, metric);
}

Weight rotated_saving(double first_u, double first_v, double second_u,
                      double second_v, WassersteinMetric metric) {
  const Weight u_delta = static_cast<Weight>(first_u) - second_u;
  if (metric == WassersteinMetric::w1_linf) {
    return Weight{2} *
               (std::min)(static_cast<Weight>(first_v),
                          static_cast<Weight>(second_v)) -
           std::fabs(u_delta);
  }
  return Weight{4} * static_cast<Weight>(first_v) * second_v -
         Weight{2} * u_delta * u_delta;
}

void dense_scalar_candidates(const PreparedDiagram& first,
                             const PreparedDiagram& second,
                             WassersteinMetric metric, WassersteinStats* stats,
                             CandidateRows& result) {
  const auto& first_points = first.finite_points();
  const auto& second_points = second.finite_points();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_v = second.finite_half_persistences();
  result.reset(first_points.size(), second_points.size());
  for (std::size_t row = 0; row < first_points.size(); ++row) {
    for (std::size_t column = 0; column < second_points.size(); ++column) {
      if (stats != nullptr) {
        ++stats->candidate_pairs;
      }
      result.add(row, column,
                 raw_saving(first_points[row], second_points[column], first_v[row],
                            second_v[column], metric),
                 stats);
    }
  }
}

void dense_parallel_candidates(const PreparedDiagram& first,
                               const PreparedDiagram& second,
                               WassersteinMetric metric,
                               WassersteinStats* stats,
                               CandidateRows& result) {
  const auto& first_points = first.finite_points();
  const auto& second_points = second.finite_points();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_v = second.finite_half_persistences();
  if (first_points.size() < 64 || second_points.empty()) {
    dense_scalar_candidates(first, second, metric, stats, result);
    return;
  }
  result.reset(first_points.size(), second_points.size());
  const std::size_t hardware =
      (std::max)(std::size_t{2},
                 static_cast<std::size_t>(std::thread::hardware_concurrency()));
  const std::size_t thread_count =
      (std::min)({std::size_t{8}, hardware,
                  (std::max)(std::size_t{2}, first_points.size() / 32)});
  std::vector<std::size_t> candidate_counts(thread_count, 0);
  std::vector<std::size_t> positive_counts(thread_count, 0);
  std::vector<std::thread> workers;
  workers.reserve(thread_count);
  for (std::size_t worker = 0; worker < thread_count; ++worker) {
    const std::size_t begin = worker * first_points.size() / thread_count;
    const std::size_t end = (worker + 1) * first_points.size() / thread_count;
    workers.emplace_back([&, worker, begin, end] {
      std::size_t candidates = 0;
      std::size_t positives = 0;
      for (std::size_t row = begin; row < end; ++row) {
        auto& edges = result.edges[row];
        for (std::size_t column = 0; column < second_points.size(); ++column) {
          ++candidates;
          const Weight saving =
              raw_saving(first_points[row], second_points[column], first_v[row],
                         second_v[column], metric);
          if (saving > Weight{0}) {
            edges.push_back({static_cast<std::uint32_t>(column),
                             static_cast<double>(saving)});
            ++positives;
          }
        }
      }
      candidate_counts[worker] = candidates;
      positive_counts[worker] = positives;
    });
  }
  for (std::thread& worker : workers) {
    worker.join();
  }
  result.edge_count =
      std::accumulate(positive_counts.begin(), positive_counts.end(), std::size_t{0});
  if (stats != nullptr) {
    stats->candidate_pairs +=
        std::accumulate(candidate_counts.begin(), candidate_counts.end(),
                        std::size_t{0});
    stats->positive_edges += result.edge_count;
  }
}

void dense_blocked_candidates(const PreparedDiagram& first,
                              const PreparedDiagram& second,
                              WassersteinMetric metric, WassersteinStats* stats,
                              CandidateRows& result) {
  const auto& first_u = first.finite_midpoints();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_u = second.finite_midpoints();
  const auto& second_v = second.finite_half_persistences();
  result.reset(first_u.size(), second_u.size());
  for (std::size_t row_block = 0; row_block < first_u.size();
       row_block += dense_block_size) {
    const std::size_t row_end =
        (std::min)(first_u.size(), row_block + dense_block_size);
    for (std::size_t column_block = 0; column_block < second_u.size();
         column_block += dense_block_size) {
      const std::size_t column_end =
          (std::min)(second_u.size(), column_block + dense_block_size);
      for (std::size_t row = row_block; row < row_end; ++row) {
        for (std::size_t column = column_block; column < column_end; ++column) {
          if (stats != nullptr) {
            ++stats->candidate_pairs;
          }
          result.add(row, column,
                     rotated_saving(first_u[row], first_v[row], second_u[column],
                                    second_v[column], metric),
                     stats);
        }
      }
    }
  }
}

void dense_avx2_candidates(const PreparedDiagram& first,
                           const PreparedDiagram& second,
                           WassersteinMetric metric, WassersteinStats* stats,
                           CandidateRows& result) {
#if defined(BOTTLENECK_HAVE_WASSERSTEIN_AVX2)
  if (cpu_has_avx2()) {
    const auto& first_u = first.finite_midpoints();
    const auto& first_v = first.finite_half_persistences();
    const auto& second_u = second.finite_midpoints();
    const auto& second_v = second.finite_half_persistences();
    result.reset(first_u.size(), second_u.size());
    std::vector<double> savings(second_u.size());
    for (std::size_t row = 0; row < first_u.size(); ++row) {
      if (metric == WassersteinMetric::w1_linf) {
        detail::fill_w1_savings_avx2(first_u[row], first_v[row], second_u.data(),
                                     second_v.data(), second_u.size(), savings.data());
      } else {
        detail::fill_w2_savings_avx2(first_u[row], first_v[row], second_u.data(),
                                     second_v.data(), second_u.size(), savings.data());
      }
      for (std::size_t column = 0; column < second_u.size(); ++column) {
        if (stats != nullptr) {
          ++stats->candidate_pairs;
        }
        result.add(row, column, static_cast<Weight>(savings[column]), stats);
      }
    }
    return;
  }
#endif
  dense_blocked_candidates(first, second, metric, stats, result);
}

void sweep_binary_candidates(const PreparedDiagram& first,
                             const PreparedDiagram& second,
                             WassersteinMetric metric, WassersteinStats* stats,
                             CandidateRows& result) {
  const auto& first_u = first.finite_midpoints();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_u = second.finite_midpoints();
  const auto& second_v = second.finite_half_persistences();
  const auto& sorted_u = second.sorted_finite_midpoints();
  const auto& order = second.finite_midpoint_order();
  result.reset(first_u.size(), second_u.size());
  const double maximum_second_v = second.max_finite_half_persistence();

  for (std::size_t row = 0; row < first_u.size(); ++row) {
    const double radius =
        metric == WassersteinMetric::w1_linf
            ? 2.0 * first_v[row]
            : std::sqrt(2.0 * first_v[row] * maximum_second_v);
    const auto begin =
        std::lower_bound(sorted_u.begin(), sorted_u.end(), first_u[row] - radius);
    const auto end =
        std::upper_bound(begin, sorted_u.end(), first_u[row] + radius);
    for (auto iterator = begin; iterator != end; ++iterator) {
      const std::size_t sorted_index =
          static_cast<std::size_t>(iterator - sorted_u.begin());
      const std::size_t column = order[sorted_index];
      if (stats != nullptr) {
        ++stats->candidate_pairs;
      }
      result.add(row, column,
                 rotated_saving(first_u[row], first_v[row], second_u[column],
                                second_v[column], metric),
                 stats);
    }
  }
}

void sweep_two_pointer_candidates(const PreparedDiagram& first,
                                  const PreparedDiagram& second,
                                  WassersteinMetric metric,
                                  WassersteinStats* stats,
                                  CandidateRows& result) {
  if (metric == WassersteinMetric::w2_l2) {
    sweep_binary_candidates(first, second, metric, stats, result);
    return;
  }
  struct Event {
    double u = 0.0;
    double v = 0.0;
    std::size_t index = 0;
    bool first = false;
  };
  const auto& first_u = first.finite_midpoints();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_u = second.finite_midpoints();
  const auto& second_v = second.finite_half_persistences();
  std::vector<Event> events;
  events.reserve(first_u.size() + second_u.size());
  for (std::size_t index = 0; index < first_u.size(); ++index) {
    events.push_back({first_u[index], first_v[index], index, true});
  }
  for (std::size_t index = 0; index < second_u.size(); ++index) {
    events.push_back({second_u[index], second_v[index], index, false});
  }
  std::stable_sort(events.begin(), events.end(), [](const Event& left, const Event& right) {
    return left.u < right.u;
  });

  result.reset(first_u.size(), second_u.size());
  std::vector<Event> active_first;
  std::vector<Event> active_second;
  for (const Event& current : events) {
    auto& opposite = current.first ? active_second : active_first;
    opposite.erase(
        std::remove_if(opposite.begin(), opposite.end(), [&](const Event& candidate) {
          return current.u - candidate.u >= 2.0 * candidate.v;
        }),
        opposite.end());
    for (const Event& candidate : opposite) {
      if (stats != nullptr) {
        ++stats->candidate_pairs;
      }
      const std::size_t row = current.first ? current.index : candidate.index;
      const std::size_t column = current.first ? candidate.index : current.index;
      result.add(row, column,
                 rotated_saving(first_u[row], first_v[row], second_u[column],
                                second_v[column], metric),
                 stats);
    }
    (current.first ? active_first : active_second).push_back(current);
  }
}

double sampled_window_density(const PreparedDiagram& first,
                              const PreparedDiagram& second,
                              WassersteinMetric metric) {
  const std::size_t sample_rows = (std::min)(std::size_t{4}, first.finite_points().size());
  if (sample_rows == 0 || second.finite_points().empty()) {
    return 0.0;
  }
  std::size_t candidates = 0;
  const auto& first_u = first.finite_midpoints();
  const auto& first_v = first.finite_half_persistences();
  const auto& sorted_second_u = second.sorted_finite_midpoints();
  for (std::size_t sample = 0; sample < sample_rows; ++sample) {
    const std::size_t row = sample * first_u.size() / sample_rows;
    const double radius =
        metric == WassersteinMetric::w1_linf
            ? 2.0 * first_v[row]
            : std::sqrt(2.0 * first_v[row] *
                        second.max_finite_half_persistence());
    const auto begin = std::lower_bound(sorted_second_u.begin(), sorted_second_u.end(),
                                        first_u[row] - radius);
    const auto end = std::upper_bound(begin, sorted_second_u.end(),
                                      first_u[row] + radius);
    candidates += static_cast<std::size_t>(end - begin);
  }
  return static_cast<double>(candidates) /
         static_cast<double>(sample_rows * sorted_second_u.size());
}

void generate_candidates(const PreparedDiagram& first,
                         const PreparedDiagram& second,
                         WassersteinConfig config, WassersteinStats* stats,
                         CandidateRows& result) {
  WassersteinCandidateStrategy strategy = config.candidates;
  if (strategy == WassersteinCandidateStrategy::adaptive) {
    const std::size_t pairs =
        first.finite_points().size() * second.finite_points().size();
    if (pairs <= 1024) {
      strategy = WassersteinCandidateStrategy::dense_scalar;
    } else {
      const double window_density =
          sampled_window_density(first, second, config.metric);
      if (window_density >= 0.75) {
        if (pairs >= 262144) {
          strategy = WassersteinCandidateStrategy::dense_parallel;
        } else {
          strategy = cpu_has_avx2() ? WassersteinCandidateStrategy::dense_avx2
                                    : WassersteinCandidateStrategy::dense_blocked;
        }
      } else {
        strategy = WassersteinCandidateStrategy::sweep_binary;
      }
    }
  }
  switch (strategy) {
    case WassersteinCandidateStrategy::dense_scalar:
      dense_scalar_candidates(first, second, config.metric, stats, result);
      return;
    case WassersteinCandidateStrategy::dense_blocked:
      dense_blocked_candidates(first, second, config.metric, stats, result);
      return;
    case WassersteinCandidateStrategy::dense_avx2:
      dense_avx2_candidates(first, second, config.metric, stats, result);
      return;
    case WassersteinCandidateStrategy::dense_parallel:
      dense_parallel_candidates(first, second, config.metric, stats, result);
      return;
    case WassersteinCandidateStrategy::sweep_binary:
      sweep_binary_candidates(first, second, config.metric, stats, result);
      return;
    case WassersteinCandidateStrategy::sweep_two_pointer:
      sweep_two_pointer_candidates(first, second, config.metric, stats, result);
      return;
    case WassersteinCandidateStrategy::topk_pricing_full_scan:
    case WassersteinCandidateStrategy::topk_pricing_simd:
    case WassersteinCandidateStrategy::topk_pricing_sweep:
    case WassersteinCandidateStrategy::topk_pricing_kdtree:
    case WassersteinCandidateStrategy::topk_pricing_kdtree_persistent:
    case WassersteinCandidateStrategy::topk_pricing_adaptive:
    case WassersteinCandidateStrategy::topk_pricing_sweep_incremental:
    case WassersteinCandidateStrategy::topk_pricing_sweep_persistent:
    case WassersteinCandidateStrategy::topk_pricing_sweep_dynamic:
    case WassersteinCandidateStrategy::topk_pricing_sweep_dynamic_batched:
      break;
    case WassersteinCandidateStrategy::adaptive:
      break;
  }
  dense_scalar_candidates(first, second, config.metric, stats, result);
}

void build_graph_into(const CandidateRows& candidates,
                      WassersteinGraphStrategy requested,
                      WassersteinStats* stats, const PreparedDiagram* first,
                      const PreparedDiagram* second, WassersteinMetric metric,
                      WeightedGraph& graph) {
  const std::uint64_t graph_bytes_before =
      stats == nullptr ? 0 : stats->graph_bytes;
  if (stats != nullptr) {
    for (const auto& row : candidates.edges) {
      stats->max_degree =
          (std::max)(stats->max_degree, static_cast<std::uint64_t>(row.size()));
    }
  }
  WassersteinGraphStrategy strategy = requested;
  const double density = candidates.rows == 0 || candidates.columns == 0
                             ? 0.0
                             : static_cast<double>(candidates.edge_count) /
                                   static_cast<double>(candidates.rows * candidates.columns);
  if (strategy == WassersteinGraphStrategy::adaptive) {
    strategy = candidates.rows * candidates.columns <= 1024 ||
                       density >= 0.15
                   ? WassersteinGraphStrategy::dense_matrix
                   : WassersteinGraphStrategy::csr;
  }

  graph.reset();
  graph.rows = candidates.rows;
  graph.columns = candidates.columns;
  graph.edge_count = candidates.edge_count;
  if (strategy == WassersteinGraphStrategy::dense_matrix) {
    graph.kind = StoredGraphKind::dense;
    graph.dense.assign(graph.rows * graph.columns, 0.0);
    for (std::size_t row = 0; row < graph.rows; ++row) {
      for (const Edge& edge : candidates.edges[row]) {
        graph.dense[row * graph.columns + edge.column] = edge.saving;
      }
    }
    if (stats != nullptr) {
      stats->graph_bytes += graph.dense.size() * sizeof(double);
    }
  } else if (strategy == WassersteinGraphStrategy::csr) {
    graph.kind = StoredGraphKind::csr;
    graph.offsets.reserve(graph.rows + 1);
    graph.offsets.push_back(0);
    graph.column_indices.reserve(graph.edge_count);
    graph.weights.reserve(graph.edge_count);
    for (const auto& row : candidates.edges) {
      for (const Edge& edge : row) {
        graph.column_indices.push_back(edge.column);
        graph.weights.push_back(edge.saving);
      }
      graph.offsets.push_back(graph.column_indices.size());
    }
    if (stats != nullptr) {
      stats->graph_bytes += graph.offsets.size() * sizeof(std::size_t) +
                            graph.column_indices.size() * sizeof(std::uint32_t) +
                            graph.weights.size() * sizeof(double);
    }
  } else if (strategy == WassersteinGraphStrategy::bitmask_lazy) {
    if (first == nullptr || second == nullptr) {
      throw std::invalid_argument(
          "bitmask_lazy graph requires prepared diagrams");
    }
    graph.kind = StoredGraphKind::bitmask_lazy;
    graph.lazy_first = first;
    graph.lazy_second = second;
    graph.lazy_metric = metric;
    graph.mask_words = (graph.columns + 63) / 64;
    graph.positive_masks.assign(graph.rows * graph.mask_words, 0);
    for (std::size_t row = 0; row < graph.rows; ++row) {
      const std::size_t base = row * graph.mask_words;
      for (const Edge& edge : candidates.edges[row]) {
        graph.positive_masks[base + edge.column / 64] |=
            std::uint64_t{1} << (edge.column % 64);
      }
    }
    if (stats != nullptr) {
      stats->graph_bytes +=
          graph.positive_masks.size() * sizeof(std::uint64_t);
    }
  } else if (strategy == WassersteinGraphStrategy::block_sparse_16 ||
             strategy == WassersteinGraphStrategy::block_sparse_32) {
    graph.kind = StoredGraphKind::block_sparse;
    graph.block_size =
        strategy == WassersteinGraphStrategy::block_sparse_16 ? 16 : 32;
    const std::size_t block_rows =
        (graph.rows + graph.block_size - 1) / graph.block_size;
    const std::size_t block_columns =
        (graph.columns + graph.block_size - 1) / graph.block_size;
    std::vector<std::uint8_t> occupied(block_rows * block_columns, 0);
    for (std::size_t row = 0; row < graph.rows; ++row) {
      const std::size_t block_row = row / graph.block_size;
      for (const Edge& edge : candidates.edges[row]) {
        occupied[block_row * block_columns + edge.column / graph.block_size] = 1;
      }
    }
    const std::size_t missing = (std::numeric_limits<std::size_t>::max)();
    std::vector<std::size_t> tile_lookup(block_rows * block_columns, missing);
    graph.block_row_offsets.reserve(block_rows + 1);
    graph.block_row_offsets.push_back(0);
    for (std::size_t block_row = 0; block_row < block_rows; ++block_row) {
      for (std::size_t block_column = 0; block_column < block_columns;
           ++block_column) {
        const std::size_t lookup = block_row * block_columns + block_column;
        if (occupied[lookup] == 0) {
          continue;
        }
        tile_lookup[lookup] = graph.block_columns.size();
        graph.block_columns.push_back(
            static_cast<std::uint32_t>(block_column));
      }
      graph.block_row_offsets.push_back(graph.block_columns.size());
    }
    graph.block_weights.assign(
        graph.block_columns.size() * graph.block_size * graph.block_size, 0.0);
    for (std::size_t row = 0; row < graph.rows; ++row) {
      const std::size_t block_row = row / graph.block_size;
      const std::size_t local_row = row % graph.block_size;
      for (const Edge& edge : candidates.edges[row]) {
        const std::size_t block_column = edge.column / graph.block_size;
        const std::size_t tile =
            tile_lookup[block_row * block_columns + block_column];
        const std::size_t local_column = edge.column % graph.block_size;
        graph.block_weights[tile * graph.block_size * graph.block_size +
                            local_row * graph.block_size + local_column] =
            edge.saving;
      }
    }
    if (stats != nullptr) {
      stats->graph_bytes +=
          graph.block_row_offsets.size() * sizeof(std::size_t) +
          graph.block_columns.size() * sizeof(std::uint32_t) +
          graph.block_weights.size() * sizeof(double);
    }
  } else if (strategy == WassersteinGraphStrategy::fixed_degree_4 ||
             strategy == WassersteinGraphStrategy::fixed_degree_8 ||
             strategy == WassersteinGraphStrategy::fixed_degree_16 ||
             strategy == WassersteinGraphStrategy::fixed_degree_32) {
    graph.kind = StoredGraphKind::fixed_degree;
    switch (strategy) {
      case WassersteinGraphStrategy::fixed_degree_4:
        graph.fixed_capacity = 4;
        break;
      case WassersteinGraphStrategy::fixed_degree_8:
        graph.fixed_capacity = 8;
        break;
      case WassersteinGraphStrategy::fixed_degree_16:
        graph.fixed_capacity = 16;
        break;
      case WassersteinGraphStrategy::fixed_degree_32:
        graph.fixed_capacity = 32;
        break;
      default:
        break;
    }
    graph.fixed_degrees.resize(graph.rows);
    graph.fixed_columns.resize(graph.rows * graph.fixed_capacity);
    graph.fixed_weights.resize(graph.rows * graph.fixed_capacity);
    graph.offsets.reserve(graph.rows + 1);
    graph.offsets.push_back(0);
    for (std::size_t row = 0; row < graph.rows; ++row) {
      const auto& edges = candidates.edges[row];
      if (edges.size() <= graph.fixed_capacity) {
        graph.fixed_degrees[row] = static_cast<std::uint8_t>(edges.size());
        const std::size_t base = row * graph.fixed_capacity;
        for (std::size_t index = 0; index < edges.size(); ++index) {
          graph.fixed_columns[base + index] = edges[index].column;
          graph.fixed_weights[base + index] = edges[index].saving;
        }
      } else {
        graph.fixed_degrees[row] = WeightedGraph::fixed_degree_overflow;
        for (const Edge& edge : edges) {
          graph.column_indices.push_back(edge.column);
          graph.weights.push_back(edge.saving);
        }
      }
      graph.offsets.push_back(graph.column_indices.size());
    }
    if (stats != nullptr) {
      stats->graph_bytes += graph.fixed_degrees.size() * sizeof(std::uint8_t) +
                            graph.fixed_columns.size() * sizeof(std::uint32_t) +
                            graph.fixed_weights.size() * sizeof(double) +
                            graph.offsets.size() * sizeof(std::size_t) +
                            graph.column_indices.size() * sizeof(std::uint32_t) +
                            graph.weights.size() * sizeof(double);
    }
  } else {
    graph.kind = StoredGraphKind::rows;
    graph.row_edges = candidates.edges;
    if (stats != nullptr) {
      stats->graph_bytes += graph.row_edges.size() * sizeof(std::vector<Edge>) +
                            graph.edge_count * sizeof(Edge);
    }
  }
  if (stats != nullptr) {
    stats->peak_graph_bytes =
        (std::max)(stats->peak_graph_bytes,
                   stats->graph_bytes - graph_bytes_before);
  }
}

WeightedGraph build_graph(const CandidateRows& candidates,
                          WassersteinGraphStrategy requested,
                          WassersteinStats* stats,
                          const PreparedDiagram* first = nullptr,
                          const PreparedDiagram* second = nullptr,
                          WassersteinMetric metric = WassersteinMetric::w1_linf) {
  WeightedGraph graph;
  build_graph_into(candidates, requested, stats, first, second, metric, graph);
  return graph;
}

struct ActiveVertices {
  std::vector<std::size_t> rows;
  std::vector<std::size_t> columns;
};

struct MatchingResult {
  explicit MatchingResult(std::size_t rows = 0)
      : row_to_column(rows, -1) {}

  Weight saving = 0;
  std::vector<int> row_to_column;
};

struct SparseDualState {
  std::vector<Weight> row_potential;
  std::vector<Weight> column_potential;
  bool valid = false;
};

bool build_matching_dual(const WeightedGraph& graph,
                         const MatchingResult& matching,
                         SparseDualState& dual);

ActiveVertices active_vertices(const WeightedGraph& graph) {
  std::vector<bool> row_active(graph.rows, false);
  std::vector<bool> column_active(graph.columns, false);
  for (std::size_t row = 0; row < graph.rows; ++row) {
    graph.for_each_edge(row, [&](std::size_t column, double) {
      row_active[row] = true;
      column_active[column] = true;
    });
  }
  ActiveVertices result;
  for (std::size_t row = 0; row < graph.rows; ++row) {
    if (row_active[row]) {
      result.rows.push_back(row);
    }
  }
  for (std::size_t column = 0; column < graph.columns; ++column) {
    if (column_active[column]) {
      result.columns.push_back(column);
    }
  }
  return result;
}

MatchingResult dense_assignment(const WeightedGraph& graph, bool compact,
                                WassersteinStats* stats) {
  std::vector<double> dense_lookup;
  if (graph.kind != StoredGraphKind::dense) {
    dense_lookup.assign(graph.rows * graph.columns, 0.0);
    for (std::size_t row = 0; row < graph.rows; ++row) {
      graph.for_each_edge(row, [&](std::size_t column, double saving) {
        dense_lookup[row * graph.columns + column] = saving;
      });
    }
  }
  const auto saving_at = [&](std::size_t row, std::size_t column) {
    return graph.kind == StoredGraphKind::dense
               ? graph.dense[row * graph.columns + column]
               : dense_lookup[row * graph.columns + column];
  };
  ActiveVertices active;
  if (compact) {
    active = active_vertices(graph);
  } else {
    active.rows.resize(graph.rows);
    active.columns.resize(graph.columns);
    std::iota(active.rows.begin(), active.rows.end(), std::size_t{0});
    std::iota(active.columns.begin(), active.columns.end(), std::size_t{0});
  }
  if (stats != nullptr) {
    stats->active_rows += active.rows.size();
    stats->active_columns += active.columns.size();
  }
  const std::size_t dimension =
      (std::max)(active.rows.size(), active.columns.size());
  if (dimension == 0) {
    return MatchingResult(graph.rows);
  }

  std::vector<Weight> row_potential(dimension + 1, Weight{0});
  std::vector<Weight> column_potential(dimension + 1, Weight{0});
  std::vector<std::size_t> matched_row(dimension + 1, 0);
  std::vector<std::size_t> predecessor(dimension + 1, 0);
  for (std::size_t row = 1; row <= dimension; ++row) {
    matched_row[0] = row;
    std::size_t column = 0;
    std::vector<Weight> minimum(dimension + 1,
                                std::numeric_limits<Weight>::infinity());
    std::vector<bool> used(dimension + 1, false);
    do {
      used[column] = true;
      const std::size_t current_row = matched_row[column];
      Weight delta = std::numeric_limits<Weight>::infinity();
      std::size_t next_column = 0;
      for (std::size_t candidate = 1; candidate <= dimension; ++candidate) {
        if (used[candidate]) {
          continue;
        }
        Weight saving = Weight{0};
        if (current_row <= active.rows.size() &&
            candidate <= active.columns.size()) {
          saving = saving_at(active.rows[current_row - 1],
                             active.columns[candidate - 1]);
        }
        const Weight reduced = -saving - row_potential[current_row] -
                               column_potential[candidate];
        if (reduced < minimum[candidate]) {
          minimum[candidate] = reduced;
          predecessor[candidate] = column;
        }
        if (minimum[candidate] < delta) {
          delta = minimum[candidate];
          next_column = candidate;
        }
      }
      for (std::size_t candidate = 0; candidate <= dimension; ++candidate) {
        if (used[candidate]) {
          row_potential[matched_row[candidate]] += delta;
          column_potential[candidate] -= delta;
        } else {
          minimum[candidate] -= delta;
        }
      }
      column = next_column;
    } while (matched_row[column] != 0);
    do {
      const std::size_t previous = predecessor[column];
      matched_row[column] = matched_row[previous];
      column = previous;
    } while (column != 0);
    if (stats != nullptr) {
      ++stats->augmentations;
    }
  }

  MatchingResult result(graph.rows);
  for (std::size_t column = 1; column <= active.columns.size(); ++column) {
    const std::size_t row = matched_row[column];
    if (row != 0 && row <= active.rows.size()) {
      const std::size_t original_row = active.rows[row - 1];
      const std::size_t original_column = active.columns[column - 1];
      const Weight saving = saving_at(original_row, original_column);
      if (saving > Weight{0}) {
        result.saving += saving;
        result.row_to_column[original_row] =
            static_cast<int>(original_column);
      }
    }
  }
  return result;
}

// Rectangular shortest-augmenting-path assignment. Unlike dense_assignment,
// this kernel never pads the smaller side to a max(rows, columns) square and
// reuses its scan buffers across augmentations. Zero-saving implicit edges make
// assigning every vertex on the smaller side equivalent to allowing unmatched
// vertices in the original positive-profit matching.
MatchingResult dense_rectangular_sap(const WeightedGraph& graph, bool compact,
                                     WassersteinStats* stats,
                                     WassersteinWorkspace::Impl* workspace,
                                     bool row_reduction,
                                     bool jv_reduction = false) {
  std::vector<double> dense_lookup;
  if (graph.kind != StoredGraphKind::dense) {
    dense_lookup.assign(graph.rows * graph.columns, 0.0);
    for (std::size_t row = 0; row < graph.rows; ++row) {
      graph.for_each_edge(row, [&](std::size_t column, double saving) {
        dense_lookup[row * graph.columns + column] = saving;
      });
    }
  }
  const auto saving_at = [&](std::size_t row, std::size_t column) {
    return graph.kind == StoredGraphKind::dense
               ? graph.dense[row * graph.columns + column]
               : dense_lookup[row * graph.columns + column];
  };

  ActiveVertices active;
  if (compact) {
    active = active_vertices(graph);
  } else {
    active.rows.resize(graph.rows);
    active.columns.resize(graph.columns);
    std::iota(active.rows.begin(), active.rows.end(), std::size_t{0});
    std::iota(active.columns.begin(), active.columns.end(), std::size_t{0});
  }
  if (stats != nullptr) {
    stats->active_rows += active.rows.size();
    stats->active_columns += active.columns.size();
  }
  if (active.rows.empty() || active.columns.empty()) {
    return MatchingResult(graph.rows);
  }

  const bool transposed = active.rows.size() > active.columns.size();
  const std::size_t short_count =
      transposed ? active.columns.size() : active.rows.size();
  const std::size_t long_count =
      transposed ? active.rows.size() : active.columns.size();
  const auto local_saving = [&](std::size_t short_index,
                                std::size_t long_index) {
    return transposed
               ? saving_at(active.rows[long_index], active.columns[short_index])
               : saving_at(active.rows[short_index], active.columns[long_index]);
  };

  std::vector<Weight> local_short_potential;
  std::vector<Weight> local_long_potential;
  std::vector<Weight> local_minimum;
  std::vector<std::size_t> local_matched_short;
  std::vector<std::size_t> local_predecessor;
  std::vector<std::uint8_t> local_used;
  std::vector<std::uint8_t> local_initially_matched;
  auto& short_potential =
      workspace == nullptr ? local_short_potential : workspace->short_potential;
  auto& long_potential =
      workspace == nullptr ? local_long_potential : workspace->long_potential;
  auto& minimum = workspace == nullptr ? local_minimum : workspace->minimum;
  auto& matched_short =
      workspace == nullptr ? local_matched_short : workspace->matched_short;
  auto& predecessor =
      workspace == nullptr ? local_predecessor : workspace->predecessor;
  auto& used = workspace == nullptr ? local_used : workspace->used;
  auto& initially_matched = workspace == nullptr
                                ? local_initially_matched
                                : workspace->initially_matched;
  short_potential.assign(short_count + 1, Weight{0});
  long_potential.assign(long_count + 1, Weight{0});
  minimum.resize(long_count + 1);
  matched_short.assign(long_count + 1, 0);
  predecessor.resize(long_count + 1);
  used.resize(long_count + 1);
  initially_matched.assign(short_count + 1, std::uint8_t{0});

  if (jv_reduction) {
    std::vector<std::size_t> row_column(short_count + 1, 0);
    for (std::size_t column = long_count; column > 0; --column) {
      Weight best_cost = std::numeric_limits<Weight>::infinity();
      std::size_t best_row = 0;
      for (std::size_t short_row = 1; short_row <= short_count; ++short_row) {
        const Weight cost = -static_cast<Weight>(
            local_saving(short_row - 1, column - 1));
        if (cost < best_cost) {
          best_cost = cost;
          best_row = short_row;
        }
      }
      long_potential[column] = best_cost;
      if (row_column[best_row] != 0) {
        matched_short[row_column[best_row]] = 0;
      }
      row_column[best_row] = column;
      matched_short[column] = best_row;
    }

    // Reduction transfer: move as much dual mass as possible from the
    // assigned column to its row while preserving feasibility and tightness.
    for (std::size_t short_row = 1; short_row <= short_count; ++short_row) {
      const std::size_t assigned_column = row_column[short_row];
      if (assigned_column == 0) {
        continue;
      }
      Weight second_reduced = std::numeric_limits<Weight>::infinity();
      for (std::size_t column = 1; column <= long_count; ++column) {
        if (column == assigned_column) {
          continue;
        }
        const Weight cost = -static_cast<Weight>(
            local_saving(short_row - 1, column - 1));
        second_reduced =
            (std::min)(second_reduced, cost - long_potential[column]);
      }
      if (!std::isfinite(second_reduced)) {
        second_reduced = Weight{0};
      }
      const Weight assigned_cost = -static_cast<Weight>(
          local_saving(short_row - 1, assigned_column - 1));
      const Weight assigned_reduced =
          assigned_cost - long_potential[assigned_column];
      const Weight transfer = second_reduced - assigned_reduced;
      long_potential[assigned_column] -= transfer;
    }

    std::vector<std::size_t> free_rows;
    free_rows.reserve(short_count);
    for (std::size_t short_row = 1; short_row <= short_count; ++short_row) {
      if (row_column[short_row] == 0) {
        free_rows.push_back(short_row);
      }
    }
    // Two augmenting-row-reduction passes. Displaced rows are handled by the
    // next pass; any rows still free afterwards enter the exact SAP loop.
    for (int pass = 0; pass < 2 && !free_rows.empty(); ++pass) {
      std::vector<std::size_t> next_free;
      next_free.reserve(free_rows.size());
      for (std::size_t short_row : free_rows) {
        Weight minimum_value = std::numeric_limits<Weight>::infinity();
        Weight second_value = std::numeric_limits<Weight>::infinity();
        std::size_t best_column = 0;
        std::size_t second_column = 0;
        for (std::size_t column = 1; column <= long_count; ++column) {
          const Weight cost = -static_cast<Weight>(
              local_saving(short_row - 1, column - 1));
          const Weight reduced = cost - long_potential[column];
          if (reduced < minimum_value) {
            second_value = minimum_value;
            second_column = best_column;
            minimum_value = reduced;
            best_column = column;
          } else if (reduced < second_value) {
            second_value = reduced;
            second_column = column;
          }
        }
        if (!std::isfinite(second_value)) {
          second_value = minimum_value;
        }
        std::size_t displaced = matched_short[best_column];
        if (minimum_value < second_value) {
          long_potential[best_column] -= second_value - minimum_value;
        } else if (displaced != 0 && second_column != 0) {
          best_column = second_column;
          displaced = matched_short[best_column];
        }
        matched_short[best_column] = short_row;
        row_column[short_row] = best_column;
        if (displaced != 0) {
          row_column[displaced] = 0;
          next_free.push_back(displaced);
        }
      }
      free_rows = std::move(next_free);
    }
    std::fill(row_column.begin(), row_column.end(), std::size_t{0});
    for (std::size_t column = 1; column <= long_count; ++column) {
      if (matched_short[column] != 0) {
        row_column[matched_short[column]] = column;
      }
    }
    for (std::size_t short_row = 1; short_row <= short_count; ++short_row) {
      const std::size_t assigned_column = row_column[short_row];
      if (assigned_column != 0) {
        const Weight cost = -static_cast<Weight>(
            local_saving(short_row - 1, assigned_column - 1));
        short_potential[short_row] =
            cost - long_potential[assigned_column];
      } else {
        Weight minimum_reduced = std::numeric_limits<Weight>::infinity();
        for (std::size_t column = 1; column <= long_count; ++column) {
          const Weight cost = -static_cast<Weight>(
              local_saving(short_row - 1, column - 1));
          minimum_reduced =
              (std::min)(minimum_reduced, cost - long_potential[column]);
        }
        short_potential[short_row] = minimum_reduced;
      }
    }
    bool feasible = true;
    std::vector<std::uint8_t> seen_rows(short_count + 1,
                                        std::uint8_t{0});
    for (std::size_t column = 1; column <= long_count && feasible; ++column) {
      const std::size_t short_row = matched_short[column];
      if (short_row == 0) {
        continue;
      }
      if (seen_rows[short_row] != 0) {
        feasible = false;
        break;
      }
      seen_rows[short_row] = 1;
      const Weight cost = -static_cast<Weight>(
          local_saving(short_row - 1, column - 1));
      const Weight reduced =
          cost - short_potential[short_row] - long_potential[column];
      const Weight scale =
          (std::max)({Weight{1}, std::fabs(cost),
                      std::fabs(short_potential[short_row]),
                      std::fabs(long_potential[column])});
      if (std::fabs(reduced) >
          Weight{256} * std::numeric_limits<Weight>::epsilon() * scale) {
        feasible = false;
      }
    }
    for (std::size_t short_row = 1;
         short_row <= short_count && feasible; ++short_row) {
      for (std::size_t column = 1; column <= long_count; ++column) {
        const Weight cost = -static_cast<Weight>(
            local_saving(short_row - 1, column - 1));
        const Weight reduced =
            cost - short_potential[short_row] - long_potential[column];
        const Weight scale =
            (std::max)({Weight{1}, std::fabs(cost),
                        std::fabs(short_potential[short_row]),
                        std::fabs(long_potential[column])});
        if (reduced <
            -Weight{256} * std::numeric_limits<Weight>::epsilon() * scale) {
          feasible = false;
          break;
        }
      }
    }
    if (!feasible) {
      short_potential.assign(short_count + 1, Weight{0});
      long_potential.assign(long_count + 1, Weight{0});
      matched_short.assign(long_count + 1, 0);
      if (stats != nullptr) {
        ++stats->jv_fallbacks;
      }
    }
    for (std::size_t column = 1; column <= long_count; ++column) {
      const std::size_t short_row = matched_short[column];
      if (short_row != 0) {
        initially_matched[short_row] = 1;
      }
    }
  } else if (row_reduction) {
    for (std::size_t short_row = 1; short_row <= short_count; ++short_row) {
      Weight best_cost = std::numeric_limits<Weight>::infinity();
      std::size_t best_column = 0;
      for (std::size_t column = 1; column <= long_count; ++column) {
        const Weight cost = -static_cast<Weight>(
            local_saving(short_row - 1, column - 1));
        if (cost < best_cost) {
          best_cost = cost;
          best_column = column;
        }
      }
      short_potential[short_row] = best_cost;
      if (matched_short[best_column] == 0) {
        matched_short[best_column] = short_row;
        initially_matched[short_row] = 1;
      }
    }
  }

  for (std::size_t short_row = 1; short_row <= short_count; ++short_row) {
    if (initially_matched[short_row] != 0) {
      continue;
    }
    matched_short[0] = short_row;
    std::size_t column = 0;
    std::fill(minimum.begin(), minimum.end(),
              std::numeric_limits<Weight>::infinity());
    std::fill(used.begin(), used.end(), std::uint8_t{0});
    do {
      used[column] = 1;
      const std::size_t current_short = matched_short[column];
      Weight delta = std::numeric_limits<Weight>::infinity();
      std::size_t next_column = 0;
      for (std::size_t candidate = 1; candidate <= long_count; ++candidate) {
        if (used[candidate] != 0) {
          continue;
        }
        const Weight saving = static_cast<Weight>(
            local_saving(current_short - 1, candidate - 1));
        const Weight reduced = -saving - short_potential[current_short] -
                               long_potential[candidate];
        if (reduced < minimum[candidate]) {
          minimum[candidate] = reduced;
          predecessor[candidate] = column;
        }
        if (minimum[candidate] < delta) {
          delta = minimum[candidate];
          next_column = candidate;
        }
      }
      for (std::size_t candidate = 0; candidate <= long_count; ++candidate) {
        if (used[candidate] != 0) {
          short_potential[matched_short[candidate]] += delta;
          long_potential[candidate] -= delta;
        } else {
          minimum[candidate] -= delta;
        }
      }
      column = next_column;
    } while (matched_short[column] != 0);
    do {
      const std::size_t previous = predecessor[column];
      matched_short[column] = matched_short[previous];
      column = previous;
    } while (column != 0);
    if (stats != nullptr) {
      ++stats->augmentations;
    }
  }

  MatchingResult result(graph.rows);
  for (std::size_t column = 1; column <= long_count; ++column) {
    const std::size_t short_row = matched_short[column];
    if (short_row != 0) {
      const std::size_t original_row =
          transposed ? active.rows[column - 1] : active.rows[short_row - 1];
      const std::size_t original_column =
          transposed ? active.columns[short_row - 1]
                     : active.columns[column - 1];
      const Weight saving = saving_at(original_row, original_column);
      if (saving > Weight{0}) {
        result.saving += saving;
        result.row_to_column[original_row] =
            static_cast<int>(original_column);
      }
    }
  }
  if (jv_reduction) {
    SparseDualState certificate;
    if (!build_matching_dual(graph, result, certificate)) {
      if (stats != nullptr) {
        ++stats->jv_fallbacks;
      }
      return dense_rectangular_sap(graph, compact, stats, workspace, false,
                                   false);
    }
  }
  return result;
}

struct ResidualEdge {
  int destination = 0;
  int reverse = 0;
  int capacity = 0;
  Weight cost = 0;
};

void add_residual_edge(std::vector<std::vector<ResidualEdge>>& network, int source,
                       int destination, Weight cost, int capacity = 1) {
  const int source_reverse = static_cast<int>(network[destination].size());
  const int destination_reverse = static_cast<int>(network[source].size());
  network[source].push_back({destination, source_reverse, capacity, cost});
  network[destination].push_back({source, destination_reverse, 0, -cost});
}

MatchingResult sparse_shortest_augmenting_path(const WeightedGraph& graph,
                                               WassersteinStats* stats,
                                               bool include_all_vertices = false,
                                               SparseDualState* dual = nullptr) {
  ActiveVertices active;
  if (include_all_vertices) {
    active.rows.resize(graph.rows);
    active.columns.resize(graph.columns);
    std::iota(active.rows.begin(), active.rows.end(), std::size_t{0});
    std::iota(active.columns.begin(), active.columns.end(), std::size_t{0});
  } else {
    active = active_vertices(graph);
  }
  if (dual != nullptr) {
    dual->valid = false;
  }
  if (stats != nullptr) {
    stats->active_rows += active.rows.size();
    stats->active_columns += active.columns.size();
  }
  if (active.rows.empty() || active.columns.empty()) {
    return MatchingResult(graph.rows);
  }
  std::vector<int> row_map(graph.rows, -1);
  std::vector<int> column_map(graph.columns, -1);
  for (std::size_t index = 0; index < active.rows.size(); ++index) {
    row_map[active.rows[index]] = static_cast<int>(index);
  }
  for (std::size_t index = 0; index < active.columns.size(); ++index) {
    column_map[active.columns[index]] = static_cast<int>(index);
  }
  const int source = 0;
  const int row_base = 1;
  const int column_base = row_base + static_cast<int>(active.rows.size());
  const int sink = column_base + static_cast<int>(active.columns.size());
  std::vector<std::vector<ResidualEdge>> network(static_cast<std::size_t>(sink + 1));
  for (std::size_t index = 0; index < active.rows.size(); ++index) {
    add_residual_edge(network, source, row_base + static_cast<int>(index), 0);
  }
  for (std::size_t index = 0; index < active.columns.size(); ++index) {
    add_residual_edge(network, column_base + static_cast<int>(index), sink, 0);
  }
  for (std::size_t row : active.rows) {
    graph.for_each_edge(row, [&](std::size_t column, double saving) {
      add_residual_edge(network, row_base + row_map[row],
                        column_base + column_map[column], -static_cast<Weight>(saving));
    });
  }

  Weight total_saving = 0;
  const std::size_t node_count = network.size();
  std::vector<Weight> potential(node_count, Weight{0});
  Weight minimum_column_potential = Weight{0};
  for (std::size_t column = 0; column < active.columns.size(); ++column) {
    Weight maximum_saving = Weight{0};
    const int column_node = column_base + static_cast<int>(column);
    for (const ResidualEdge& reverse_edge : network[column_node]) {
      if (reverse_edge.destination >= row_base &&
          reverse_edge.destination < column_base) {
        maximum_saving = (std::max)(maximum_saving, reverse_edge.cost);
      }
    }
    potential[static_cast<std::size_t>(column_node)] = -maximum_saving;
    minimum_column_potential =
        (std::min)(minimum_column_potential, -maximum_saving);
  }
  potential[static_cast<std::size_t>(sink)] = minimum_column_potential;

  while (true) {
    std::vector<Weight> distance(node_count,
                                 std::numeric_limits<Weight>::infinity());
    std::vector<int> previous_node(node_count, -1);
    std::vector<int> previous_edge(node_count, -1);
    using QueueItem = std::pair<Weight, int>;
    std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<>> queue;
    distance[source] = 0;
    queue.push({0, source});
    while (!queue.empty()) {
      const auto [queued_distance, node] = queue.top();
      queue.pop();
      if (queued_distance != distance[static_cast<std::size_t>(node)]) {
        continue;
      }
      for (std::size_t edge_index = 0; edge_index < network[node].size(); ++edge_index) {
        const ResidualEdge& edge = network[node][edge_index];
        if (edge.capacity == 0) {
          continue;
        }
        // Feasible primal-dual potentials make every reduced cost non-negative.
        // Clamp a negative zero caused by double rounding on MSVC; a materially
        // invalid potential falls back to the dense exact solver below.
        const Weight raw_reduced =
            edge.cost + potential[static_cast<std::size_t>(node)] -
            potential[static_cast<std::size_t>(edge.destination)];
        if (!std::isfinite(raw_reduced)) {
          if (stats != nullptr) {
            ++stats->sparse_fallbacks;
          }
          return dense_assignment(graph, true, stats);
        }
        const Weight reduced = (std::max)(Weight{0}, raw_reduced);
        const Weight candidate =
            distance[static_cast<std::size_t>(node)] + reduced;
        if (distance[static_cast<std::size_t>(edge.destination)] <= candidate) {
          continue;
        }
        distance[static_cast<std::size_t>(edge.destination)] = candidate;
        previous_node[static_cast<std::size_t>(edge.destination)] = node;
        previous_edge[static_cast<std::size_t>(edge.destination)] =
            static_cast<int>(edge_index);
        queue.push({candidate, edge.destination});
      }
    }
    if (previous_node[static_cast<std::size_t>(sink)] < 0) {
      break;
    }
    const Weight path_cost = distance[static_cast<std::size_t>(sink)] -
                             potential[static_cast<std::size_t>(source)] +
                             potential[static_cast<std::size_t>(sink)];
    if (path_cost >= Weight{0}) {
      for (std::size_t node = 0; node < node_count; ++node) {
        if (std::isfinite(distance[node])) {
          potential[node] += distance[node];
        }
      }
      break;
    }
    for (std::size_t node = 0; node < node_count; ++node) {
      if (std::isfinite(distance[node])) {
        potential[node] += distance[node];
      }
    }
    for (int node = sink; node != source;
         node = previous_node[static_cast<std::size_t>(node)]) {
      const int parent = previous_node[static_cast<std::size_t>(node)];
      const int edge_index = previous_edge[static_cast<std::size_t>(node)];
      ResidualEdge& edge = network[parent][static_cast<std::size_t>(edge_index)];
      --edge.capacity;
      ++network[node][static_cast<std::size_t>(edge.reverse)].capacity;
    }
    total_saving -= path_cost;
    if (stats != nullptr) {
      ++stats->augmentations;
    }
  }
  MatchingResult result(graph.rows);
  result.saving = total_saving;
  for (std::size_t local_row = 0; local_row < active.rows.size(); ++local_row) {
    const int row_node = row_base + static_cast<int>(local_row);
    for (const ResidualEdge& edge : network[static_cast<std::size_t>(row_node)]) {
      if (edge.destination >= column_base && edge.destination < sink &&
          edge.capacity == 0) {
        const std::size_t local_column =
            static_cast<std::size_t>(edge.destination - column_base);
        result.row_to_column[active.rows[local_row]] =
            static_cast<int>(active.columns[local_column]);
        break;
      }
    }
  }
  if (dual != nullptr) {
    dual->valid = build_matching_dual(graph, result, *dual);
  }
  return result;
}

struct ArenaResidualEdge {
  int destination = 0;
  std::size_t reverse = 0;
  int capacity = 0;
  Weight cost = 0;
};

// Exact sparse SAP experiment with the same primal-dual path search as the
// vector-of-vectors baseline. The topology lives in one exactly sized edge
// arena, while Dijkstra buffers and heap capacity are reused across
// augmentations. This remains an explicit matcher until end-to-end evidence
// shows a stable win.
MatchingResult sparse_shortest_augmenting_path_arena(
    const WeightedGraph& graph, WassersteinStats* stats,
    bool include_all_vertices = false) {
  ActiveVertices active;
  if (include_all_vertices) {
    active.rows.resize(graph.rows);
    active.columns.resize(graph.columns);
    std::iota(active.rows.begin(), active.rows.end(), std::size_t{0});
    std::iota(active.columns.begin(), active.columns.end(), std::size_t{0});
  } else {
    active = active_vertices(graph);
  }
  if (stats != nullptr) {
    stats->active_rows += active.rows.size();
    stats->active_columns += active.columns.size();
  }
  if (active.rows.empty() || active.columns.empty()) {
    return MatchingResult(graph.rows);
  }

  std::vector<int> row_map(graph.rows, -1);
  std::vector<int> column_map(graph.columns, -1);
  for (std::size_t index = 0; index < active.rows.size(); ++index) {
    row_map[active.rows[index]] = static_cast<int>(index);
  }
  for (std::size_t index = 0; index < active.columns.size(); ++index) {
    column_map[active.columns[index]] = static_cast<int>(index);
  }

  const int source = 0;
  const int row_base = 1;
  const int column_base = row_base + static_cast<int>(active.rows.size());
  const int sink = column_base + static_cast<int>(active.columns.size());
  const std::size_t node_count = static_cast<std::size_t>(sink + 1);
  std::vector<std::size_t> degree(node_count, 0);
  for (std::size_t index = 0; index < active.rows.size(); ++index) {
    ++degree[static_cast<std::size_t>(source)];
    ++degree[static_cast<std::size_t>(row_base + static_cast<int>(index))];
  }
  for (std::size_t index = 0; index < active.columns.size(); ++index) {
    ++degree[static_cast<std::size_t>(column_base + static_cast<int>(index))];
    ++degree[static_cast<std::size_t>(sink)];
  }
  for (std::size_t row : active.rows) {
    graph.for_each_edge(row, [&](std::size_t column, double) {
      ++degree[static_cast<std::size_t>(row_base + row_map[row])];
      ++degree[static_cast<std::size_t>(column_base + column_map[column])];
    });
  }

  std::vector<std::size_t> offsets(node_count + 1, 0);
  for (std::size_t node = 0; node < node_count; ++node) {
    offsets[node + 1] = offsets[node] + degree[node];
  }
  std::vector<std::size_t> cursor(offsets.begin(), offsets.end() - 1);
  std::vector<ArenaResidualEdge> edges(offsets.back());
  const auto add_arena_edge = [&](int from, int destination, Weight cost) {
    const std::size_t forward = cursor[static_cast<std::size_t>(from)]++;
    const std::size_t reverse = cursor[static_cast<std::size_t>(destination)]++;
    edges[forward] = {destination, reverse, 1, cost};
    edges[reverse] = {from, forward, 0, -cost};
  };
  for (std::size_t index = 0; index < active.rows.size(); ++index) {
    add_arena_edge(source, row_base + static_cast<int>(index), Weight{0});
  }
  for (std::size_t index = 0; index < active.columns.size(); ++index) {
    add_arena_edge(column_base + static_cast<int>(index), sink, Weight{0});
  }
  for (std::size_t row : active.rows) {
    graph.for_each_edge(row, [&](std::size_t column, double saving) {
      add_arena_edge(row_base + row_map[row], column_base + column_map[column],
                     -static_cast<Weight>(saving));
    });
  }

  std::vector<Weight> potential(node_count, Weight{0});
  std::vector<Weight> distance(node_count);
  std::vector<int> previous_node(node_count);
  std::vector<std::size_t> previous_edge(node_count);
  using QueueItem = std::pair<Weight, int>;
  std::vector<QueueItem> heap;
  heap.reserve((std::min)(edges.size(), node_count * 8));
  const auto arena_bytes = [&]() {
    return degree.capacity() * sizeof(degree.front()) +
           offsets.capacity() * sizeof(offsets.front()) +
           cursor.capacity() * sizeof(cursor.front()) +
           edges.capacity() * sizeof(edges.front()) +
           potential.capacity() * sizeof(potential.front()) +
           distance.capacity() * sizeof(distance.front()) +
           previous_node.capacity() * sizeof(previous_node.front()) +
           previous_edge.capacity() * sizeof(previous_edge.front()) +
           heap.capacity() * sizeof(QueueItem);
  };
  const auto record_arena_bytes = [&]() {
    if (stats != nullptr) {
      stats->peak_sparse_arena_bytes =
          (std::max)(stats->peak_sparse_arena_bytes,
                     static_cast<std::uint64_t>(arena_bytes()));
    }
  };
  if (stats != nullptr) {
    ++stats->sparse_arena_builds;
  }
  record_arena_bytes();

  Weight minimum_column_potential = Weight{0};
  for (std::size_t column = 0; column < active.columns.size(); ++column) {
    Weight maximum_saving = Weight{0};
    const std::size_t column_node = static_cast<std::size_t>(
        column_base + static_cast<int>(column));
    for (std::size_t edge_index = offsets[column_node];
         edge_index < offsets[column_node + 1]; ++edge_index) {
      const ArenaResidualEdge& reverse_edge = edges[edge_index];
      if (reverse_edge.destination >= row_base &&
          reverse_edge.destination < column_base) {
        maximum_saving = (std::max)(maximum_saving, reverse_edge.cost);
      }
    }
    potential[column_node] = -maximum_saving;
    minimum_column_potential =
        (std::min)(minimum_column_potential, -maximum_saving);
  }
  potential[static_cast<std::size_t>(sink)] = minimum_column_potential;

  Weight total_saving = 0;
  bool scratch_initialized = false;
  const std::greater<QueueItem> compare;
  while (true) {
    if (scratch_initialized && stats != nullptr) {
      ++stats->sparse_scratch_reuses;
    }
    scratch_initialized = true;
    std::fill(distance.begin(), distance.end(),
              std::numeric_limits<Weight>::infinity());
    std::fill(previous_node.begin(), previous_node.end(), -1);
    std::fill(previous_edge.begin(), previous_edge.end(), std::size_t{0});
    heap.clear();
    const auto push = [&](QueueItem item) {
      const std::size_t old_capacity = heap.capacity();
      heap.push_back(item);
      if (heap.capacity() != old_capacity && stats != nullptr) {
        ++stats->sparse_heap_growths;
        record_arena_bytes();
      }
      std::push_heap(heap.begin(), heap.end(), compare);
    };
    distance[static_cast<std::size_t>(source)] = Weight{0};
    push({Weight{0}, source});
    while (!heap.empty()) {
      std::pop_heap(heap.begin(), heap.end(), compare);
      const auto [queued_distance, node] = heap.back();
      heap.pop_back();
      if (queued_distance != distance[static_cast<std::size_t>(node)]) {
        continue;
      }
      const std::size_t node_index = static_cast<std::size_t>(node);
      for (std::size_t edge_index = offsets[node_index];
           edge_index < offsets[node_index + 1]; ++edge_index) {
        const ArenaResidualEdge& edge = edges[edge_index];
        if (edge.capacity == 0) {
          continue;
        }
        const Weight raw_reduced =
            edge.cost + potential[node_index] -
            potential[static_cast<std::size_t>(edge.destination)];
        if (!std::isfinite(raw_reduced)) {
          if (stats != nullptr) {
            ++stats->sparse_fallbacks;
          }
          return dense_assignment(graph, true, stats);
        }
        const Weight reduced = (std::max)(Weight{0}, raw_reduced);
        const Weight candidate = distance[node_index] + reduced;
        const std::size_t destination =
            static_cast<std::size_t>(edge.destination);
        if (distance[destination] <= candidate) {
          continue;
        }
        distance[destination] = candidate;
        previous_node[destination] = node;
        previous_edge[destination] = edge_index;
        push({candidate, edge.destination});
      }
    }
    if (previous_node[static_cast<std::size_t>(sink)] < 0) {
      break;
    }
    const Weight path_cost = distance[static_cast<std::size_t>(sink)] -
                             potential[static_cast<std::size_t>(source)] +
                             potential[static_cast<std::size_t>(sink)];
    for (std::size_t node = 0; node < node_count; ++node) {
      if (std::isfinite(distance[node])) {
        potential[node] += distance[node];
      }
    }
    if (path_cost >= Weight{0}) {
      break;
    }
    for (int node = sink; node != source;
         node = previous_node[static_cast<std::size_t>(node)]) {
      ArenaResidualEdge& edge =
          edges[previous_edge[static_cast<std::size_t>(node)]];
      --edge.capacity;
      ++edges[edge.reverse].capacity;
    }
    total_saving -= path_cost;
    if (stats != nullptr) {
      ++stats->augmentations;
    }
  }

  MatchingResult result(graph.rows);
  result.saving = total_saving;
  for (std::size_t local_row = 0; local_row < active.rows.size(); ++local_row) {
    const std::size_t row_node = static_cast<std::size_t>(
        row_base + static_cast<int>(local_row));
    for (std::size_t edge_index = offsets[row_node];
         edge_index < offsets[row_node + 1]; ++edge_index) {
      const ArenaResidualEdge& edge = edges[edge_index];
      if (edge.destination >= column_base && edge.destination < sink &&
          edge.capacity == 0) {
        const std::size_t local_column =
            static_cast<std::size_t>(edge.destination - column_base);
        result.row_to_column[active.rows[local_row]] =
            static_cast<int>(active.columns[local_column]);
        break;
      }
    }
  }
  record_arena_bytes();
  return result;
}

bool build_matching_dual(const WeightedGraph& graph,
                         const MatchingResult& matching,
                         SparseDualState& dual) {
  struct Constraint {
    std::size_t destination = 0;
    Weight delta = 0;
  };
  dual.row_potential.assign(graph.rows, Weight{0});
  dual.column_potential.assign(graph.columns, Weight{0});
  std::vector<int> column_match(graph.columns, -1);
  std::vector<Weight> matched_saving(graph.columns, Weight{0});
  std::vector<std::uint8_t> matched_row(graph.rows, std::uint8_t{0});
  for (std::size_t row = 0; row < graph.rows; ++row) {
    const int column = matching.row_to_column[row];
    if (column < 0) {
      continue;
    }
    matched_row[row] = 1;
    column_match[static_cast<std::size_t>(column)] = static_cast<int>(row);
    graph.for_each_edge(row, [&](std::size_t candidate, double saving) {
      if (candidate == static_cast<std::size_t>(column)) {
        matched_saving[candidate] = static_cast<Weight>(saving);
      }
    });
  }

  std::vector<std::vector<Constraint>> constraints(graph.rows);
  for (std::size_t row = 0; row < graph.rows; ++row) {
    graph.for_each_edge(row, [&](std::size_t column, double saving_value) {
      const Weight saving = static_cast<Weight>(saving_value);
      const int matched = column_match[column];
      if (matched < 0) {
        dual.row_potential[row] =
            (std::max)(dual.row_potential[row], saving);
      } else {
        const std::size_t source = static_cast<std::size_t>(matched);
        constraints[source].push_back(
            {row, saving - matched_saving[column]});
      }
    });
  }

  std::queue<std::size_t> queue;
  std::vector<std::uint8_t> queued(graph.rows, std::uint8_t{1});
  std::vector<std::size_t> relaxations(graph.rows, 0);
  for (std::size_t row = 0; row < graph.rows; ++row) {
    if (matched_row[row] == 0 && dual.row_potential[row] > Weight{0}) {
      return false;
    }
    queue.push(row);
  }
  while (!queue.empty()) {
    const std::size_t source = queue.front();
    queue.pop();
    queued[source] = 0;
    for (const Constraint& constraint : constraints[source]) {
      const Weight candidate =
          dual.row_potential[source] + constraint.delta;
      const Weight scale =
          (std::max)({Weight{1}, std::fabs(candidate),
                      std::fabs(dual.row_potential[constraint.destination])});
      if (candidate <= dual.row_potential[constraint.destination] +
                           Weight{64} *
                               std::numeric_limits<Weight>::epsilon() * scale) {
        continue;
      }
      if (matched_row[constraint.destination] == 0) {
        return false;
      }
      dual.row_potential[constraint.destination] = candidate;
      if (++relaxations[constraint.destination] > graph.rows) {
        return false;
      }
      if (queued[constraint.destination] == 0) {
        queued[constraint.destination] = 1;
        queue.push(constraint.destination);
      }
    }
  }

  for (std::size_t column = 0; column < graph.columns; ++column) {
    const int row = column_match[column];
    if (row < 0) {
      continue;
    }
    const Weight value =
        matched_saving[column] -
        dual.row_potential[static_cast<std::size_t>(row)];
    const Weight scale =
        (std::max)({Weight{1}, std::fabs(value),
                    std::fabs(matched_saving[column])});
    if (value < -Weight{64} * std::numeric_limits<Weight>::epsilon() * scale) {
      return false;
    }
    dual.column_potential[column] = (std::max)(Weight{0}, value);
  }

  bool feasible = true;
  for (std::size_t row = 0; row < graph.rows && feasible; ++row) {
    graph.for_each_edge(row, [&](std::size_t column, double saving_value) {
      const Weight saving = static_cast<Weight>(saving_value);
      const Weight bound =
          dual.row_potential[row] + dual.column_potential[column];
      const Weight scale =
          (std::max)({Weight{1}, std::fabs(saving), std::fabs(bound)});
      if (saving - bound >
          Weight{128} * std::numeric_limits<Weight>::epsilon() * scale) {
        feasible = false;
      }
    });
  }
  const Weight dual_objective =
      std::accumulate(dual.row_potential.begin(), dual.row_potential.end(),
                      Weight{0}) +
      std::accumulate(dual.column_potential.begin(), dual.column_potential.end(),
                      Weight{0});
  const Weight objective_scale =
      (std::max)({Weight{1}, std::fabs(dual_objective),
                  std::fabs(matching.saving)});
  return feasible &&
         std::fabs(dual_objective - matching.saving) <=
             Weight{512} * std::numeric_limits<Weight>::epsilon() *
                 objective_scale;
}

MatchingResult sparse_cycle_reoptimize(const WeightedGraph& graph,
                                       const MatchingResult& initial,
                                       WassersteinStats* stats,
                                       SparseDualState& dual) {
  const int source = 0;
  const int row_base = 1;
  const int column_base = row_base + static_cast<int>(graph.rows);
  const int sink = column_base + static_cast<int>(graph.columns);
  std::vector<std::vector<ResidualEdge>> network(static_cast<std::size_t>(sink + 1));
  std::vector<std::size_t> source_edges(graph.rows);
  std::vector<std::size_t> sink_edges(graph.columns);
  std::vector<std::size_t> matched_edges(
      graph.rows, (std::numeric_limits<std::size_t>::max)());
  for (std::size_t row = 0; row < graph.rows; ++row) {
    source_edges[row] = network[static_cast<std::size_t>(source)].size();
    add_residual_edge(network, source, row_base + static_cast<int>(row), 0);
  }
  for (std::size_t column = 0; column < graph.columns; ++column) {
    const int column_node = column_base + static_cast<int>(column);
    sink_edges[column] = network[static_cast<std::size_t>(column_node)].size();
    add_residual_edge(network, column_node, sink, 0);
  }
  for (std::size_t row = 0; row < graph.rows; ++row) {
    const int row_node = row_base + static_cast<int>(row);
    const int matched_column = initial.row_to_column[row];
    graph.for_each_edge(row, [&](std::size_t column, double saving) {
      const std::size_t edge_index =
          network[static_cast<std::size_t>(row_node)].size();
      add_residual_edge(network, row_node,
                        column_base + static_cast<int>(column),
                        -static_cast<Weight>(saving));
      if (matched_column >= 0 &&
          column == static_cast<std::size_t>(matched_column)) {
        matched_edges[row] = edge_index;
      }
    });
  }
  const std::size_t virtual_edge = network[static_cast<std::size_t>(sink)].size();
  add_residual_edge(network, sink, source, 0,
                    static_cast<int>((std::min)(graph.rows, graph.columns)));
  const auto send_flow = [&](int from, std::size_t edge_index) {
    ResidualEdge& edge = network[static_cast<std::size_t>(from)][edge_index];
    --edge.capacity;
    ++network[static_cast<std::size_t>(edge.destination)]
             [static_cast<std::size_t>(edge.reverse)]
                 .capacity;
  };
  for (std::size_t row = 0; row < graph.rows; ++row) {
    const int column = initial.row_to_column[row];
    if (column < 0 ||
        matched_edges[row] == (std::numeric_limits<std::size_t>::max)()) {
      continue;
    }
    send_flow(source, source_edges[row]);
    send_flow(row_base + static_cast<int>(row), matched_edges[row]);
    send_flow(column_base + column,
              sink_edges[static_cast<std::size_t>(column)]);
    send_flow(sink, virtual_edge);
  }

  const std::size_t node_count = network.size();
  while (true) {
    std::vector<Weight> distance(node_count, Weight{0});
    std::vector<int> previous_node(node_count, -1);
    std::vector<int> previous_edge(node_count, -1);
    int updated = -1;
    for (std::size_t pass = 0; pass < node_count; ++pass) {
      updated = -1;
      for (std::size_t node = 0; node < node_count; ++node) {
        for (std::size_t edge_index = 0; edge_index < network[node].size();
             ++edge_index) {
          const ResidualEdge& edge = network[node][edge_index];
          if (edge.capacity == 0) {
            continue;
          }
          const Weight candidate = distance[node] + edge.cost;
          const Weight scale =
              (std::max)({Weight{1}, std::fabs(candidate),
                          std::fabs(distance[static_cast<std::size_t>(
                              edge.destination)])});
          if (candidate + Weight{64} *
                                  std::numeric_limits<Weight>::epsilon() * scale >=
              distance[static_cast<std::size_t>(edge.destination)]) {
            continue;
          }
          distance[static_cast<std::size_t>(edge.destination)] = candidate;
          previous_node[static_cast<std::size_t>(edge.destination)] =
              static_cast<int>(node);
          previous_edge[static_cast<std::size_t>(edge.destination)] =
              static_cast<int>(edge_index);
          updated = edge.destination;
        }
      }
      if (updated < 0) {
        break;
      }
    }
    if (updated < 0) {
      break;
    }
    int cycle = updated;
    for (std::size_t step = 0; step < node_count; ++step) {
      cycle = previous_node[static_cast<std::size_t>(cycle)];
      if (cycle < 0) {
        break;
      }
    }
    if (cycle < 0) {
      break;
    }
    int flow = (std::numeric_limits<int>::max)();
    int node = cycle;
    do {
      const int parent = previous_node[static_cast<std::size_t>(node)];
      const int edge_index = previous_edge[static_cast<std::size_t>(node)];
      flow = (std::min)(flow,
                        network[static_cast<std::size_t>(parent)]
                               [static_cast<std::size_t>(edge_index)]
                                   .capacity);
      node = parent;
    } while (node != cycle);
    node = cycle;
    do {
      const int parent = previous_node[static_cast<std::size_t>(node)];
      const int edge_index = previous_edge[static_cast<std::size_t>(node)];
      ResidualEdge& edge = network[static_cast<std::size_t>(parent)]
                                  [static_cast<std::size_t>(edge_index)];
      edge.capacity -= flow;
      network[static_cast<std::size_t>(node)]
             [static_cast<std::size_t>(edge.reverse)]
                 .capacity += flow;
      node = parent;
    } while (node != cycle);
    if (stats != nullptr) {
      ++stats->augmentations;
    }
  }

  MatchingResult result(graph.rows);
  for (std::size_t row = 0; row < graph.rows; ++row) {
    const int row_node = row_base + static_cast<int>(row);
    for (const ResidualEdge& edge : network[static_cast<std::size_t>(row_node)]) {
      if (edge.destination >= column_base && edge.destination < sink &&
          edge.capacity == 0) {
        result.row_to_column[row] = edge.destination - column_base;
        result.saving -= edge.cost;
        break;
      }
    }
  }
  dual.valid = build_matching_dual(graph, result, dual);
  return result;
}

class PersistentCycleReoptimizer {
 public:
  PersistentCycleReoptimizer(std::size_t rows, std::size_t columns)
      : rows_(rows), columns_(columns), row_base_(1),
        column_base_(row_base_ + static_cast<int>(rows)),
        sink_(column_base_ + static_cast<int>(columns)),
        network_(static_cast<std::size_t>(sink_ + 1)),
        source_edges_(rows), sink_edges_(columns), known_(rows) {
    for (std::size_t row = 0; row < rows_; ++row) {
      source_edges_[row] = network_[0].size();
      add_residual_edge(network_, 0, row_base_ + static_cast<int>(row), 0);
    }
    for (std::size_t column = 0; column < columns_; ++column) {
      const int node = column_base_ + static_cast<int>(column);
      sink_edges_[column] = network_[static_cast<std::size_t>(node)].size();
      add_residual_edge(network_, node, sink_, 0);
    }
    virtual_edge_ = network_[static_cast<std::size_t>(sink_)].size();
    add_residual_edge(network_, sink_, 0, 0,
                      static_cast<int>((std::min)(rows_, columns_)));
    const std::size_t nodes = network_.size();
    distance_.resize(nodes);
    previous_node_.resize(nodes);
    previous_edge_.resize(nodes);
  }

  void initialize(const CandidateRows& candidates,
                  const MatchingResult& matching,
                  const SparseDualState* dual = nullptr) {
    append(candidates);
    for (std::size_t row = 0; row < rows_; ++row) {
      const int column = matching.row_to_column[row];
      if (column < 0) {
        continue;
      }
      const auto iterator = std::lower_bound(
          known_[row].begin(), known_[row].end(),
          static_cast<std::uint32_t>(column),
          [](const KnownEdge& edge, std::uint32_t target) {
            return edge.column < target;
          });
      if (iterator == known_[row].end() ||
          iterator->column != static_cast<std::uint32_t>(column)) {
        continue;
      }
      send_flow(0, source_edges_[row]);
      send_flow(row_base_ + static_cast<int>(row), iterator->edge_index);
      send_flow(column_base_ + column,
                sink_edges_[static_cast<std::size_t>(column)]);
      send_flow(sink_, virtual_edge_);
    }
    potential_.assign(network_.size(), Weight{0});
    if (dual != nullptr && dual->valid &&
        dual->row_potential.size() == rows_ &&
        dual->column_potential.size() == columns_) {
      for (std::size_t row = 0; row < rows_; ++row) {
        potential_[static_cast<std::size_t>(row_base_) + row] =
            dual->row_potential[row];
      }
      for (std::size_t column = 0; column < columns_; ++column) {
        potential_[static_cast<std::size_t>(column_base_) + column] =
            -dual->column_potential[column];
      }
    }
  }

  void append(const CandidateRows& candidates) {
    for (std::size_t row = 0; row < rows_; ++row) {
      auto& known = known_[row];
      for (const Edge& edge : candidates.edges[row]) {
        const auto position = std::lower_bound(
            known.begin(), known.end(), edge.column,
            [](const KnownEdge& item, std::uint32_t target) {
              return item.column < target;
            });
        if (position != known.end() && position->column == edge.column) {
          continue;
        }
        const int node = row_base_ + static_cast<int>(row);
        const std::size_t edge_index =
            network_[static_cast<std::size_t>(node)].size();
        add_residual_edge(network_, node,
                          column_base_ + static_cast<int>(edge.column),
                          -static_cast<Weight>(edge.saving));
        known.insert(position, {edge.column, edge_index});
      }
    }
  }

  MatchingResult reoptimize(WassersteinStats* stats) {
    const std::size_t node_count = network_.size();
    while (true) {
      std::fill(distance_.begin(), distance_.end(), Weight{0});
      std::fill(previous_node_.begin(), previous_node_.end(), -1);
      std::fill(previous_edge_.begin(), previous_edge_.end(), -1);
      int updated = -1;
      for (std::size_t pass = 0; pass < node_count; ++pass) {
        updated = -1;
        for (std::size_t node = 0; node < node_count; ++node) {
          for (std::size_t edge_index = 0; edge_index < network_[node].size();
               ++edge_index) {
            const ResidualEdge& edge = network_[node][edge_index];
            if (edge.capacity == 0) {
              continue;
            }
            const Weight candidate = distance_[node] + edge.cost;
            const Weight scale =
                (std::max)({Weight{1}, std::fabs(candidate),
                            std::fabs(distance_[static_cast<std::size_t>(
                                edge.destination)])});
            if (candidate + Weight{64} *
                                    std::numeric_limits<Weight>::epsilon() *
                                    scale >=
                distance_[static_cast<std::size_t>(edge.destination)]) {
              continue;
            }
            distance_[static_cast<std::size_t>(edge.destination)] = candidate;
            previous_node_[static_cast<std::size_t>(edge.destination)] =
                static_cast<int>(node);
            previous_edge_[static_cast<std::size_t>(edge.destination)] =
                static_cast<int>(edge_index);
            updated = edge.destination;
          }
        }
        if (updated < 0) {
          break;
        }
      }
      if (updated < 0) {
        break;
      }
      int cycle = updated;
      for (std::size_t step = 0; step < node_count; ++step) {
        cycle = previous_node_[static_cast<std::size_t>(cycle)];
        if (cycle < 0) {
          break;
        }
      }
      if (cycle < 0) {
        break;
      }
      int flow = (std::numeric_limits<int>::max)();
      int node = cycle;
      do {
        const int parent = previous_node_[static_cast<std::size_t>(node)];
        const int edge_index = previous_edge_[static_cast<std::size_t>(node)];
        flow = (std::min)(
            flow, network_[static_cast<std::size_t>(parent)]
                          [static_cast<std::size_t>(edge_index)]
                              .capacity);
        node = parent;
      } while (node != cycle);
      node = cycle;
      do {
        const int parent = previous_node_[static_cast<std::size_t>(node)];
        const int edge_index = previous_edge_[static_cast<std::size_t>(node)];
        ResidualEdge& edge =
            network_[static_cast<std::size_t>(parent)]
                    [static_cast<std::size_t>(edge_index)];
        edge.capacity -= flow;
        network_[static_cast<std::size_t>(node)]
                [static_cast<std::size_t>(edge.reverse)]
                    .capacity += flow;
        node = parent;
      } while (node != cycle);
      if (stats != nullptr) {
        ++stats->augmentations;
      }
    }
    return matching();
  }

  bool append_and_reoptimize(const CandidateRows& candidates,
                             WassersteinStats* stats) {
    using QueueItem = std::pair<Weight, int>;
    for (std::size_t row = 0; row < rows_; ++row) {
      auto& known = known_[row];
      for (const Edge& candidate_edge : candidates.edges[row]) {
        const auto position = std::lower_bound(
            known.begin(), known.end(), candidate_edge.column,
            [](const KnownEdge& item, std::uint32_t target) {
              return item.column < target;
            });
        if (position != known.end() &&
            position->column == candidate_edge.column) {
          continue;
        }
        const int from = row_base_ + static_cast<int>(row);
        const int to = column_base_ + static_cast<int>(candidate_edge.column);
        const std::size_t edge_index =
            network_[static_cast<std::size_t>(from)].size();
        add_residual_edge(network_, from, to,
                          -static_cast<Weight>(candidate_edge.saving));
        known.insert(position, {candidate_edge.column, edge_index});
        if (stats != nullptr) {
          ++stats->dynamic_inserted_edges;
        }

        const ResidualEdge& inserted =
            network_[static_cast<std::size_t>(from)][edge_index];
        const Weight inserted_reduced =
            inserted.cost + potential_[static_cast<std::size_t>(from)] -
            potential_[static_cast<std::size_t>(to)];
        const Weight inserted_scale =
            (std::max)({Weight{1}, std::fabs(inserted.cost),
                        std::fabs(potential_[static_cast<std::size_t>(from)]),
                        std::fabs(potential_[static_cast<std::size_t>(to)])});
        const Weight tolerance =
            Weight{128} * std::numeric_limits<Weight>::epsilon() *
            inserted_scale;
        if (inserted_reduced >= -tolerance) {
          continue;
        }

        if (stats != nullptr) {
          ++stats->dynamic_dijkstra_runs;
        }

        std::fill(distance_.begin(), distance_.end(),
                  std::numeric_limits<Weight>::infinity());
        std::fill(previous_node_.begin(), previous_node_.end(), -1);
        std::fill(previous_edge_.begin(), previous_edge_.end(), -1);
        std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<>>
            queue;
        distance_[static_cast<std::size_t>(to)] = Weight{0};
        queue.push({Weight{0}, to});
        while (!queue.empty()) {
          const auto [queued_distance, node] = queue.top();
          queue.pop();
          if (queued_distance != distance_[static_cast<std::size_t>(node)]) {
            continue;
          }
          for (std::size_t residual_index = 0;
               residual_index <
               network_[static_cast<std::size_t>(node)].size();
               ++residual_index) {
            if (node == from && residual_index == edge_index) {
              continue;
            }
            const ResidualEdge& edge =
                network_[static_cast<std::size_t>(node)][residual_index];
            if (edge.capacity == 0) {
              continue;
            }
            const Weight raw_reduced =
                edge.cost + potential_[static_cast<std::size_t>(node)] -
                potential_[static_cast<std::size_t>(edge.destination)];
            const Weight scale =
                (std::max)({Weight{1}, std::fabs(edge.cost),
                            std::fabs(potential_[static_cast<std::size_t>(node)]),
                            std::fabs(potential_[static_cast<std::size_t>(
                                edge.destination)])});
            const Weight edge_tolerance =
                Weight{256} * std::numeric_limits<Weight>::epsilon() * scale;
            if (raw_reduced < -edge_tolerance) {
              return false;
            }
            const Weight reduced = (std::max)(Weight{0}, raw_reduced);
            const Weight next = queued_distance + reduced;
            if (next >=
                distance_[static_cast<std::size_t>(edge.destination)]) {
              continue;
            }
            distance_[static_cast<std::size_t>(edge.destination)] = next;
            previous_node_[static_cast<std::size_t>(edge.destination)] = node;
            previous_edge_[static_cast<std::size_t>(edge.destination)] =
                static_cast<int>(residual_index);
            queue.push({next, edge.destination});
          }
        }
        const Weight path_reduced = distance_[static_cast<std::size_t>(from)];
        if (!std::isfinite(path_reduced)) {
          return false;
        }
        // Cap the distance shift at the inserted edge's endpoint distance.
        // Updating only reachable nodes (or using their full distances) can
        // make an edge entering the reached set acquire negative reduced cost.
        // The capped shift preserves dual feasibility for every old residual
        // edge while still making the chosen shortest path tight.
        for (std::size_t node = 0; node < network_.size(); ++node) {
          potential_[node] +=
              std::isfinite(distance_[node])
                  ? (std::min)(distance_[node], path_reduced)
                  : path_reduced;
        }
        if (inserted_reduced + path_reduced >= -tolerance) {
          continue;
        }
        send_flow(from, edge_index);
        for (int node = from; node != to;
             node = previous_node_[static_cast<std::size_t>(node)]) {
          const int parent = previous_node_[static_cast<std::size_t>(node)];
          const int residual_index =
              previous_edge_[static_cast<std::size_t>(node)];
          if (parent < 0 || residual_index < 0) {
            return false;
          }
          send_flow(parent, static_cast<std::size_t>(residual_index));
        }
        if (stats != nullptr) {
          ++stats->augmentations;
        }
      }
    }
    return true;
  }

  bool append_and_reoptimize_batched(const CandidateRows& candidates,
                                     WassersteinStats* stats) {
    using QueueItem = std::pair<Weight, int>;
    struct PendingEdge {
      std::size_t row = 0;
      std::size_t edge_index = 0;
      Weight tolerance = 0;
    };

    std::vector<std::vector<PendingEdge>> groups(columns_);
    for (std::size_t row = 0; row < rows_; ++row) {
      auto& known = known_[row];
      for (const Edge& candidate_edge : candidates.edges[row]) {
        const auto position = std::lower_bound(
            known.begin(), known.end(), candidate_edge.column,
            [](const KnownEdge& item, std::uint32_t target) {
              return item.column < target;
            });
        if (position != known.end() &&
            position->column == candidate_edge.column) {
          continue;
        }
        const int from = row_base_ + static_cast<int>(row);
        const int to = column_base_ + static_cast<int>(candidate_edge.column);
        const std::size_t edge_index =
            network_[static_cast<std::size_t>(from)].size();
        add_residual_edge(network_, from, to,
                          -static_cast<Weight>(candidate_edge.saving));
        known.insert(position, {candidate_edge.column, edge_index});
        if (stats != nullptr) {
          ++stats->dynamic_inserted_edges;
        }

        const ResidualEdge& inserted =
            network_[static_cast<std::size_t>(from)][edge_index];
        const Weight inserted_reduced =
            inserted.cost + potential_[static_cast<std::size_t>(from)] -
            potential_[static_cast<std::size_t>(to)];
        const Weight inserted_scale =
            (std::max)({Weight{1}, std::fabs(inserted.cost),
                        std::fabs(potential_[static_cast<std::size_t>(from)]),
                        std::fabs(potential_[static_cast<std::size_t>(to)])});
        const Weight tolerance =
            Weight{128} * std::numeric_limits<Weight>::epsilon() *
            inserted_scale;
        if (inserted_reduced < -tolerance) {
          groups[static_cast<std::size_t>(candidate_edge.column)].push_back(
              {row, edge_index, tolerance});
        }
      }
    }

    std::vector<std::vector<std::uint8_t>> blocked_edges(rows_);
    for (std::size_t row = 0; row < rows_; ++row) {
      blocked_edges[row].resize(
          network_[static_cast<std::size_t>(row_base_) + row].size(),
          std::uint8_t{0});
    }
    for (const auto& group : groups) {
      for (const PendingEdge& edge : group) {
        blocked_edges[edge.row][edge.edge_index] = std::uint8_t{1};
      }
    }
    for (std::size_t column = 0; column < columns_; ++column) {
      auto& pending = groups[column];
      if (pending.empty()) {
        continue;
      }
      if (stats != nullptr) {
        ++stats->dynamic_batch_groups;
      }
      const int to = column_base_ + static_cast<int>(column);
      while (!pending.empty()) {
        std::fill(distance_.begin(), distance_.end(),
                  std::numeric_limits<Weight>::infinity());
        std::fill(previous_node_.begin(), previous_node_.end(), -1);
        std::fill(previous_edge_.begin(), previous_edge_.end(), -1);
        std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<>>
            queue;
        distance_[static_cast<std::size_t>(to)] = Weight{0};
        queue.push({Weight{0}, to});
        if (stats != nullptr) {
          ++stats->dynamic_dijkstra_runs;
        }
        while (!queue.empty()) {
          const auto [queued_distance, node] = queue.top();
          queue.pop();
          if (queued_distance != distance_[static_cast<std::size_t>(node)]) {
            continue;
          }
          for (std::size_t residual_index = 0;
               residual_index <
               network_[static_cast<std::size_t>(node)].size();
               ++residual_index) {
            if (node >= row_base_ && node < column_base_ &&
                blocked_edges[static_cast<std::size_t>(node - row_base_)]
                             [residual_index] != 0) {
              continue;
            }
            const ResidualEdge& edge =
                network_[static_cast<std::size_t>(node)][residual_index];
            if (edge.capacity == 0) {
              continue;
            }
            const Weight raw_reduced =
                edge.cost + potential_[static_cast<std::size_t>(node)] -
                potential_[static_cast<std::size_t>(edge.destination)];
            const Weight scale =
                (std::max)({Weight{1}, std::fabs(edge.cost),
                            std::fabs(potential_[static_cast<std::size_t>(node)]),
                            std::fabs(potential_[static_cast<std::size_t>(
                                edge.destination)])});
            const Weight edge_tolerance =
                Weight{256} * std::numeric_limits<Weight>::epsilon() * scale;
            if (raw_reduced < -edge_tolerance) {
              return false;
            }
            const Weight reduced = (std::max)(Weight{0}, raw_reduced);
            const Weight next = queued_distance + reduced;
            if (next >=
                distance_[static_cast<std::size_t>(edge.destination)]) {
              continue;
            }
            distance_[static_cast<std::size_t>(edge.destination)] = next;
            previous_node_[static_cast<std::size_t>(edge.destination)] = node;
            previous_edge_[static_cast<std::size_t>(edge.destination)] =
                static_cast<int>(residual_index);
            queue.push({next, edge.destination});
          }
        }

        std::size_t chosen = pending.size();
        Weight best_cycle = Weight{0};
        for (std::size_t index = 0; index < pending.size(); ++index) {
          const PendingEdge& item = pending[index];
          const int from = row_base_ + static_cast<int>(item.row);
          const Weight path_reduced =
              distance_[static_cast<std::size_t>(from)];
          if (!std::isfinite(path_reduced)) {
            continue;
          }
          const ResidualEdge& inserted =
              network_[static_cast<std::size_t>(from)][item.edge_index];
          const Weight inserted_reduced =
              inserted.cost + potential_[static_cast<std::size_t>(from)] -
              potential_[static_cast<std::size_t>(to)];
          const Weight cycle_reduced = inserted_reduced + path_reduced;
          if (cycle_reduced < -item.tolerance &&
              (chosen == pending.size() || cycle_reduced < best_cycle)) {
            chosen = index;
            best_cycle = cycle_reduced;
          }
        }

        if (chosen == pending.size()) {
          Weight required_shift = Weight{0};
          for (const PendingEdge& item : pending) {
            const int from = row_base_ + static_cast<int>(item.row);
            const ResidualEdge& inserted =
                network_[static_cast<std::size_t>(from)][item.edge_index];
            const Weight inserted_reduced =
                inserted.cost + potential_[static_cast<std::size_t>(from)] -
                potential_[static_cast<std::size_t>(to)];
            required_shift =
                (std::max)(required_shift, -inserted_reduced);
          }
          for (std::size_t node = 0; node < network_.size(); ++node) {
            potential_[node] +=
                std::isfinite(distance_[node])
                    ? (std::min)(distance_[node], required_shift)
                    : required_shift;
          }
          for (const PendingEdge& item : pending) {
            blocked_edges[item.row][item.edge_index] = std::uint8_t{0};
          }
          pending.clear();
          continue;
        }

        const PendingEdge selected = pending[chosen];
        const int from = row_base_ + static_cast<int>(selected.row);
        const Weight path_reduced = distance_[static_cast<std::size_t>(from)];
        for (std::size_t node = 0; node < network_.size(); ++node) {
          potential_[node] +=
              std::isfinite(distance_[node])
                  ? (std::min)(distance_[node], path_reduced)
                  : path_reduced;
        }
        blocked_edges[selected.row][selected.edge_index] = std::uint8_t{0};
        send_flow(from, selected.edge_index);
        for (int node = from; node != to;
             node = previous_node_[static_cast<std::size_t>(node)]) {
          const int parent = previous_node_[static_cast<std::size_t>(node)];
          const int residual_index =
              previous_edge_[static_cast<std::size_t>(node)];
          if (parent < 0 || residual_index < 0) {
            return false;
          }
          send_flow(parent, static_cast<std::size_t>(residual_index));
        }
        pending.erase(pending.begin() + static_cast<std::ptrdiff_t>(chosen));
        if (stats != nullptr) {
          ++stats->augmentations;
        }
      }
    }
    return true;
  }

  MatchingResult current_matching() const { return matching(); }

 private:
  struct KnownEdge {
    std::uint32_t column = 0;
    std::size_t edge_index = 0;
  };

  void send_flow(int from, std::size_t edge_index) {
    ResidualEdge& edge =
        network_[static_cast<std::size_t>(from)][edge_index];
    --edge.capacity;
    ++network_[static_cast<std::size_t>(edge.destination)]
             [static_cast<std::size_t>(edge.reverse)]
                 .capacity;
  }

  MatchingResult matching() const {
    MatchingResult result(rows_);
    for (std::size_t row = 0; row < rows_; ++row) {
      const int node = row_base_ + static_cast<int>(row);
      for (const KnownEdge& known : known_[row]) {
        const ResidualEdge& edge =
            network_[static_cast<std::size_t>(node)][known.edge_index];
        if (edge.capacity == 0) {
          result.row_to_column[row] = static_cast<int>(known.column);
          result.saving -= edge.cost;
          break;
        }
      }
    }
    return result;
  }

  std::size_t rows_ = 0;
  std::size_t columns_ = 0;
  int row_base_ = 0;
  int column_base_ = 0;
  int sink_ = 0;
  std::vector<std::vector<ResidualEdge>> network_;
  std::vector<std::size_t> source_edges_;
  std::vector<std::size_t> sink_edges_;
  std::size_t virtual_edge_ = 0;
  std::vector<std::vector<KnownEdge>> known_;
  std::vector<Weight> distance_;
  std::vector<int> previous_node_;
  std::vector<int> previous_edge_;
  std::vector<Weight> potential_;
};

using TopEdgeKey = std::pair<Weight, std::int64_t>;
using TopEdgeHeap =
    std::priority_queue<TopEdgeKey, std::vector<TopEdgeKey>,
                        std::greater<TopEdgeKey>>;

void retain_top_edge(TopEdgeHeap& heap, std::size_t limit, Weight score,
                     std::size_t column) {
  const TopEdgeKey key{score, -static_cast<std::int64_t>(column)};
  if (heap.size() < limit) {
    heap.push(key);
  } else if (key > heap.top()) {
    heap.pop();
    heap.push(key);
  }
}

template <class Function>
void for_each_pricing_candidate(const PreparedDiagram& first,
                                const PreparedDiagram& second,
                                WassersteinMetric metric, std::size_t row,
                                bool full_scan, Function&& function) {
  if (full_scan) {
    for (std::size_t column = 0; column < second.finite_points().size();
         ++column) {
      function(column);
    }
    return;
  }
  const auto& first_u = first.finite_midpoints();
  const auto& first_v = first.finite_half_persistences();
  const auto& sorted_second_u = second.sorted_finite_midpoints();
  const auto& order = second.finite_midpoint_order();
  const double radius =
      metric == WassersteinMetric::w1_linf
          ? 2.0 * first_v[row]
          : std::sqrt(2.0 * first_v[row] *
                      second.max_finite_half_persistence());
  const auto begin = std::lower_bound(sorted_second_u.begin(), sorted_second_u.end(),
                                      first_u[row] - radius);
  const auto end = std::upper_bound(begin, sorted_second_u.end(),
                                    first_u[row] + radius);
  for (auto iterator = begin; iterator != end; ++iterator) {
    function(order[static_cast<std::size_t>(iterator - sorted_second_u.begin())]);
  }
}

void topk_seed_candidates(const PreparedDiagram& first,
                          const PreparedDiagram& second,
                          WassersteinMetric metric, std::size_t requested_top_k,
                          bool full_scan, WassersteinStats* stats,
                          CandidateRows& candidates) {
  const auto& first_points = first.finite_points();
  const auto& second_points = second.finite_points();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_v = second.finite_half_persistences();
  candidates.reset(first_points.size(), second_points.size());
  const std::size_t top_k =
      (std::max)(std::size_t{1}, (std::min)(requested_top_k, second_points.size()));
  std::uint64_t examined = 0;
  std::uint64_t positive_edges = 0;
  for (std::size_t row = 0; row < first_points.size(); ++row) {
    TopEdgeHeap heap;
    for_each_pricing_candidate(first, second, metric, row, full_scan,
                               [&](std::size_t column) {
      ++examined;
      const Weight saving = raw_saving(first_points[row], second_points[column],
                                       first_v[row], second_v[column], metric);
      if (saving <= Weight{0}) {
        return;
      }
      ++positive_edges;
      retain_top_edge(heap, top_k, saving, column);
    });
    auto& edges = candidates.edges[row];
    edges.reserve(heap.size());
    while (!heap.empty()) {
      const std::size_t column =
          static_cast<std::size_t>(-heap.top().second);
      edges.push_back({static_cast<std::uint32_t>(column),
                       static_cast<double>(heap.top().first)});
      heap.pop();
    }
    std::sort(edges.begin(), edges.end(),
              [](const Edge& left, const Edge& right) {
                return left.column < right.column;
              });
    candidates.edge_count += edges.size();
  }
  if (stats != nullptr) {
    stats->candidate_pairs += examined;
    stats->positive_edges += positive_edges;
    stats->peak_materialized_edges =
        (std::max)(stats->peak_materialized_edges,
                   static_cast<std::uint64_t>(candidates.edge_count));
  }
}

class PricingKdTree {
 public:
  explicit PricingKdTree(const PreparedDiagram& diagram) : diagram_(diagram) {
    const std::size_t columns = diagram.finite_points().size();
    while (base_ < columns) {
      base_ *= 2;
    }
    minimum_u_.assign(base_ * 2, std::numeric_limits<double>::infinity());
    maximum_u_.assign(base_ * 2, -std::numeric_limits<double>::infinity());
    maximum_v_.assign(base_ * 2, 0.0);
    minimum_beta_.assign(base_ * 2,
                         std::numeric_limits<Weight>::infinity());
    const auto& order = diagram.finite_midpoint_order();
    const auto& u = diagram.finite_midpoints();
    const auto& v = diagram.finite_half_persistences();
    for (std::size_t position = 0; position < columns; ++position) {
      const std::size_t column = order[position];
      const std::size_t node = base_ + position;
      minimum_u_[node] = u[column];
      maximum_u_[node] = u[column];
      maximum_v_[node] = v[column];
    }
    for (std::size_t node = base_; node-- > 1;) {
      const std::size_t left = node * 2;
      const std::size_t right = left + 1;
      minimum_u_[node] = (std::min)(minimum_u_[left], minimum_u_[right]);
      maximum_u_[node] = (std::max)(maximum_u_[left], maximum_u_[right]);
      maximum_v_[node] = (std::max)(maximum_v_[left], maximum_v_[right]);
    }
  }

  PricingKdTree(const PreparedDiagram& diagram,
                const std::vector<Weight>& column_potential)
      : PricingKdTree(diagram) {
    update_column_potential(column_potential);
  }

  void update_column_potential(
      const std::vector<Weight>& column_potential) {
    const std::size_t columns = diagram_.finite_points().size();
    const auto& order = diagram_.finite_midpoint_order();
    for (std::size_t position = 0; position < columns; ++position) {
      minimum_beta_[base_ + position] = column_potential[order[position]];
    }
    for (std::size_t node = base_; node-- > 1;) {
      minimum_beta_[node] =
          (std::min)(minimum_beta_[node * 2], minimum_beta_[node * 2 + 1]);
    }
  }

  template <class Function>
  void for_each_possible_violation(WassersteinMetric metric, double row_u,
                                   double row_v, Weight row_potential,
                                   Function&& function,
                                   std::vector<std::size_t>* scratch = nullptr) const {
    if (diagram_.finite_points().empty()) {
      return;
    }
    const auto& order = diagram_.finite_midpoint_order();
    std::vector<std::size_t> local_stack;
    std::vector<std::size_t>& stack =
        scratch == nullptr ? local_stack : *scratch;
    stack.clear();
    stack.push_back(1);
    while (!stack.empty()) {
      const std::size_t node = stack.back();
      stack.pop_back();
      if (!could_violate(node, metric, row_u, row_v, row_potential)) {
        continue;
      }
      if (node >= base_) {
        const std::size_t position = node - base_;
        if (position < order.size()) {
          function(order[position]);
        }
        continue;
      }
      stack.push_back(node * 2 + 1);
      stack.push_back(node * 2);
    }
  }

 private:
  bool could_violate(std::size_t node, WassersteinMetric metric, double row_u,
                     double row_v, Weight row_potential) const {
    if (!std::isfinite(minimum_u_[node])) {
      return false;
    }
    double minimum_du = 0.0;
    if (row_u < minimum_u_[node]) {
      minimum_du = minimum_u_[node] - row_u;
    } else if (row_u > maximum_u_[node]) {
      minimum_du = row_u - maximum_u_[node];
    }
    const Weight saving_upper =
        metric == WassersteinMetric::w1_linf
            ? static_cast<Weight>(
                  2.0 * (std::min)(row_v, maximum_v_[node]) - minimum_du)
            : static_cast<Weight>(4.0 * row_v * maximum_v_[node] -
                                  2.0 * minimum_du * minimum_du);
    const Weight reduced_upper =
        saving_upper - row_potential - minimum_beta_[node];
    const Weight scale =
        (std::max)({Weight{1}, std::fabs(saving_upper),
                    std::fabs(row_potential), std::fabs(minimum_beta_[node])});
    // The bound is assembled from binary64 node summaries. Keep a small
    // outward margin so a rounded-down upper bound cannot hide a genuinely
    // positive reduced saving.
    return reduced_upper >=
           -Weight{128} * std::numeric_limits<Weight>::epsilon() * scale;
  }

  const PreparedDiagram& diagram_;
  std::size_t base_ = 1;
  std::vector<double> minimum_u_;
  std::vector<double> maximum_u_;
  std::vector<double> maximum_v_;
  std::vector<Weight> minimum_beta_;
};

struct PricingRoundResult {
  std::size_t added = 0;
  std::uint64_t priced = 0;
  std::uint64_t violations = 0;
};

PricingRoundResult price_and_add_edges(
    const PreparedDiagram& first, const PreparedDiagram& second,
    WassersteinMetric metric, std::size_t requested_top_k, bool full_scan,
    bool kd_tree, bool simd_scan, const SparseDualState& dual,
    WassersteinStats* stats, CandidateRows& candidates,
    PricingKdTree* reusable_pricing_tree = nullptr,
    std::vector<std::size_t>* reusable_tree_stack = nullptr) {
  const auto& first_points = first.finite_points();
  const auto& second_points = second.finite_points();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_v = second.finite_half_persistences();
  const std::size_t top_k =
      (std::max)(std::size_t{1}, (std::min)(requested_top_k, second_points.size()));
  std::size_t added = 0;
  std::uint64_t priced = 0;
  std::uint64_t violations = 0;
  const bool use_simd = simd_scan && cpu_has_avx2();
  std::vector<double> simd_savings(use_simd ? second_points.size() : 0);
  std::optional<PricingKdTree> local_pricing_tree;
  PricingKdTree* pricing_tree = nullptr;
  if (kd_tree && dual.valid) {
    if (reusable_pricing_tree != nullptr) {
      reusable_pricing_tree->update_column_potential(dual.column_potential);
      pricing_tree = reusable_pricing_tree;
      if (stats != nullptr) {
        ++stats->pricing_kdtree_updates;
      }
    } else {
      local_pricing_tree.emplace(second, dual.column_potential);
      pricing_tree = &*local_pricing_tree;
      if (stats != nullptr) {
        ++stats->pricing_kdtree_builds;
      }
    }
  }
  for (std::size_t row = 0; row < first_points.size(); ++row) {
    auto& edges = candidates.edges[row];
    std::vector<Edge> pending;
    TopEdgeHeap heap;
    const auto inspect_saving = [&](std::size_t column, Weight saving) {
      const auto included = std::lower_bound(
          edges.begin(), edges.end(), column,
          [](const Edge& edge, std::size_t target) {
            return edge.column < target;
          });
      if (included != edges.end() && included->column == column) {
        return;
      }
      ++priced;
      if (saving <= Weight{0}) {
        return;
      }
      if (!dual.valid) {
        pending.push_back({static_cast<std::uint32_t>(column),
                           static_cast<double>(saving)});
        ++violations;
        return;
      }
      const Weight reduced = dual.row_potential[row] +
                             dual.column_potential[column] - saving;
      if (reduced < Weight{0}) {
        ++violations;
        retain_top_edge(heap, top_k, -reduced, column);
      }
    };
    const auto inspect = [&](std::size_t column) {
      inspect_saving(column,
                     raw_saving(first_points[row], second_points[column],
                                first_v[row], second_v[column], metric));
    };
    if (use_simd) {
#if defined(BOTTLENECK_HAVE_WASSERSTEIN_AVX2)
      if (metric == WassersteinMetric::w1_linf) {
        detail::fill_w1_savings_avx2(
            first.finite_midpoints()[row], first_v[row],
            second.finite_midpoints().data(), second_v.data(),
            second_points.size(), simd_savings.data());
      } else {
        detail::fill_w2_savings_avx2(
            first.finite_midpoints()[row], first_v[row],
            second.finite_midpoints().data(), second_v.data(),
            second_points.size(), simd_savings.data());
      }
      for (std::size_t column = 0; column < second_points.size(); ++column) {
        inspect_saving(column, static_cast<Weight>(simd_savings[column]));
      }
#endif
    } else if (pricing_tree != nullptr) {
      pricing_tree->for_each_possible_violation(
          metric, first.finite_midpoints()[row], first_v[row],
          dual.row_potential[row], inspect, reusable_tree_stack);
    } else {
      for_each_pricing_candidate(first, second, metric, row, full_scan,
                                 inspect);
    }
    while (!heap.empty()) {
      const std::size_t column =
          static_cast<std::size_t>(-heap.top().second);
      const Weight saving = raw_saving(first_points[row], second_points[column],
                                       first_v[row], second_v[column], metric);
      pending.push_back({static_cast<std::uint32_t>(column),
                         static_cast<double>(saving)});
      heap.pop();
    }
    if (!pending.empty()) {
      added += pending.size();
      edges.insert(edges.end(), pending.begin(), pending.end());
      std::sort(edges.begin(), edges.end(),
                [](const Edge& left, const Edge& right) {
                  return left.column < right.column;
                });
    }
  }
  candidates.edge_count += added;
  if (stats != nullptr) {
    ++stats->pricing_rounds;
    if (use_simd) {
      ++stats->pricing_simd_rounds;
    } else if (full_scan) {
      ++stats->pricing_full_scan_rounds;
    } else if (kd_tree) {
      ++stats->pricing_kdtree_rounds;
    } else {
      ++stats->pricing_sweep_rounds;
    }
    stats->priced_edges += priced;
    stats->pricing_violations += violations;
    stats->peak_materialized_edges =
        (std::max)(stats->peak_materialized_edges,
                   static_cast<std::uint64_t>(candidates.edge_count));
  }
  return {added, priced, violations};
}

MatchingResult solve_priced_restricted_graph(const WeightedGraph& graph,
                                             WassersteinStats* stats,
                                             SparseDualState& dual);

MatchingResult priced_topk_matching(const PreparedDiagram& first,
                                    const PreparedDiagram& second,
                                    WassersteinMetric metric, std::size_t top_k,
                                    bool full_scan, bool kd_tree, bool simd_scan,
                                    bool persistent_kd_tree, bool adaptive_pricing,
                                    bool incremental,
                                    WassersteinStats* stats,
                                    CandidateRows& candidates,
                                    WeightedGraph& graph) {
  const auto candidate_start = Clock::now();
  topk_seed_candidates(first, second, metric, top_k, full_scan, stats,
                       candidates);
  if (stats != nullptr) {
    stats->candidate_time_ns += elapsed_ns(candidate_start, Clock::now());
  }
  const double window_density =
      adaptive_pricing ? sampled_window_density(first, second, metric) : 0.0;
  std::uint64_t previous_priced = 0;
  std::uint64_t previous_violations = 0;
  std::optional<PricingKdTree> reusable_pricing_tree;
  std::vector<std::size_t> reusable_tree_stack;
  if (persistent_kd_tree) {
    const auto tree_start = Clock::now();
    reusable_pricing_tree.emplace(second);
    reusable_tree_stack.reserve(64);
    if (stats != nullptr) {
      ++stats->pricing_kdtree_builds;
      stats->pricing_time_ns += elapsed_ns(tree_start, Clock::now());
    }
  }
  MatchingResult matching(first.finite_points().size());
  bool initial_solve = true;
  while (true) {
    const auto graph_start = Clock::now();
    build_graph_into(candidates, WassersteinGraphStrategy::csr, stats, &first,
                     &second, metric, graph);
    if (stats != nullptr) {
      stats->graph_time_ns += elapsed_ns(graph_start, Clock::now());
    }

    SparseDualState dual;
    const auto solver_start = Clock::now();
    if (incremental && !initial_solve) {
      matching = sparse_cycle_reoptimize(graph, matching, stats, dual);
    } else {
      matching = solve_priced_restricted_graph(graph, stats, dual);
    }
    initial_solve = false;
    if (stats != nullptr) {
      stats->solver_time_ns += elapsed_ns(solver_start, Clock::now());
    }

    const auto pricing_start = Clock::now();
    bool round_full_scan = full_scan;
    bool round_kd_tree = kd_tree;
    bool round_simd_scan = simd_scan;
    if (adaptive_pricing) {
      round_full_scan = window_density >= 0.75;
      // Exact-guarded SIMD pricing is retained as an explicit E1 experiment.
      // Its terminal scalar certificate erased the local scan win, so the E4
      // router deliberately excludes it and combines only stable oracle paths.
      round_simd_scan = false;
      const double previous_hit_rate =
          previous_priced == 0
              ? 1.0
              : static_cast<double>(previous_violations) /
                    static_cast<double>(previous_priced);
      round_kd_tree = !round_full_scan &&
                      (std::max)(first.finite_points().size(),
                                 second.finite_points().size()) >= 256 &&
                      window_density >= 0.02 && previous_priced != 0 &&
                      previous_hit_rate < 0.02;
    }
    PricingRoundResult pricing = price_and_add_edges(
        first, second, metric, top_k, round_full_scan, round_kd_tree,
        round_simd_scan, dual, stats, candidates,
        reusable_pricing_tree.has_value() ? &*reusable_pricing_tree : nullptr,
        reusable_pricing_tree.has_value() ? &reusable_tree_stack : nullptr);
    if (pricing.added == 0 && round_simd_scan) {
      // SIMD saving evaluation is a fast binary64 filter, not the final exact
      // certificate. A scalar full scan on the terminal SIMD round prevents a
      // rounded-down saving from hiding the last violated edge.
      pricing = price_and_add_edges(first, second, metric, top_k, true, false,
                                    false, dual, stats, candidates);
    }
    previous_priced = pricing.priced;
    previous_violations = pricing.violations;
    if (stats != nullptr) {
      stats->pricing_time_ns += elapsed_ns(pricing_start, Clock::now());
    }
    if (pricing.added == 0) {
      return matching;
    }
  }
}

MatchingResult priced_topk_persistent_matching(
    const PreparedDiagram& first, const PreparedDiagram& second,
    WassersteinMetric metric, std::size_t top_k, WassersteinStats* stats,
    CandidateRows& candidates, WeightedGraph& graph) {
  const auto candidate_start = Clock::now();
  topk_seed_candidates(first, second, metric, top_k, false, stats, candidates);
  if (stats != nullptr) {
    stats->candidate_time_ns += elapsed_ns(candidate_start, Clock::now());
  }
  const auto graph_start = Clock::now();
  build_graph_into(candidates, WassersteinGraphStrategy::csr, stats, &first,
                   &second, metric, graph);
  if (stats != nullptr) {
    stats->graph_time_ns += elapsed_ns(graph_start, Clock::now());
  }
  SparseDualState dual;
  const auto solver_start = Clock::now();
  MatchingResult matching = solve_priced_restricted_graph(graph, stats, dual);
  PersistentCycleReoptimizer persistent(first.finite_points().size(),
                                        second.finite_points().size());
  persistent.initialize(candidates, matching);
  if (stats != nullptr) {
    stats->solver_time_ns += elapsed_ns(solver_start, Clock::now());
  }
  while (true) {
    const auto pricing_start = Clock::now();
    const PricingRoundResult pricing = price_and_add_edges(
        first, second, metric, top_k, false, false, false, dual, stats,
        candidates);
    if (stats != nullptr) {
      stats->pricing_time_ns += elapsed_ns(pricing_start, Clock::now());
    }
    if (pricing.added == 0) {
      return matching;
    }
    const auto incremental_start = Clock::now();
    persistent.append(candidates);
    matching = persistent.reoptimize(stats);
    const auto rebuild_start = Clock::now();
    build_graph_into(candidates, WassersteinGraphStrategy::csr, stats, &first,
                     &second, metric, graph);
    if (stats != nullptr) {
      stats->graph_time_ns += elapsed_ns(rebuild_start, Clock::now());
    }
    dual.valid = build_matching_dual(graph, matching, dual);
    if (stats != nullptr) {
      stats->solver_time_ns += elapsed_ns(incremental_start, Clock::now());
    }
  }
}

MatchingResult priced_topk_dynamic_matching(
    const PreparedDiagram& first, const PreparedDiagram& second,
    WassersteinMetric metric, std::size_t top_k, bool batched,
    WassersteinStats* stats, CandidateRows& candidates, WeightedGraph& graph) {
  const auto candidate_start = Clock::now();
  topk_seed_candidates(first, second, metric, top_k, false, stats, candidates);
  if (stats != nullptr) {
    stats->candidate_time_ns += elapsed_ns(candidate_start, Clock::now());
  }
  const auto graph_start = Clock::now();
  build_graph_into(candidates, WassersteinGraphStrategy::csr, stats, &first,
                   &second, metric, graph);
  if (stats != nullptr) {
    stats->graph_time_ns += elapsed_ns(graph_start, Clock::now());
  }
  SparseDualState dual;
  const auto solver_start = Clock::now();
  MatchingResult matching = solve_priced_restricted_graph(graph, stats, dual);
  PersistentCycleReoptimizer dynamic(first.finite_points().size(),
                                     second.finite_points().size());
  dynamic.initialize(candidates, matching, &dual);
  if (stats != nullptr) {
    stats->solver_time_ns += elapsed_ns(solver_start, Clock::now());
  }
  while (true) {
    const auto pricing_start = Clock::now();
    const PricingRoundResult pricing = price_and_add_edges(
        first, second, metric, top_k, false, false, false, dual, stats,
        candidates);
    if (stats != nullptr) {
      stats->pricing_time_ns += elapsed_ns(pricing_start, Clock::now());
    }
    if (pricing.added == 0) {
      return matching;
    }
    const auto incremental_start = Clock::now();
    const bool incremental_ok =
        dual.valid &&
        (batched ? dynamic.append_and_reoptimize_batched(candidates, stats)
                 : dynamic.append_and_reoptimize(candidates, stats));
    if (incremental_ok) {
      matching = dynamic.current_matching();
    }
    if (stats != nullptr) {
      stats->solver_time_ns += elapsed_ns(incremental_start, Clock::now());
    }

    const auto rebuild_start = Clock::now();
    build_graph_into(candidates, WassersteinGraphStrategy::csr, stats, &first,
                     &second, metric, graph);
    if (stats != nullptr) {
      stats->graph_time_ns += elapsed_ns(rebuild_start, Clock::now());
    }

    bool certified = false;
    if (incremental_ok) {
      const auto certificate_start = Clock::now();
      certified = build_matching_dual(graph, matching, dual);
      dual.valid = certified;
      if (stats != nullptr) {
        stats->solver_time_ns += elapsed_ns(certificate_start, Clock::now());
      }
    }
    if (!incremental_ok || !certified) {
      const auto fallback_start = Clock::now();
      if (stats != nullptr) {
        ++stats->sparse_fallbacks;
      }
      matching = solve_priced_restricted_graph(graph, stats, dual);
      dynamic = PersistentCycleReoptimizer(first.finite_points().size(),
                                           second.finite_points().size());
      dynamic.initialize(candidates, matching, &dual);
      if (stats != nullptr) {
        stats->solver_time_ns += elapsed_ns(fallback_start, Clock::now());
      }
    }
  }
}

std::optional<MatchingResult> certified_greedy(const WeightedGraph& graph,
                                               WassersteinWarmStart strategy,
                                               WassersteinStats* stats) {
  if (strategy == WassersteinWarmStart::none || graph.edge_count == 0) {
    return std::nullopt;
  }
  struct RankedEdge {
    std::size_t row = 0;
    std::size_t column = 0;
    double saving = 0.0;
  };
  std::vector<RankedEdge> edges;
  edges.reserve(graph.edge_count);
  std::vector<double> row_max(graph.rows, 0.0);
  std::vector<double> column_max(graph.columns, 0.0);
  for (std::size_t row = 0; row < graph.rows; ++row) {
    graph.for_each_edge(row, [&](std::size_t column, double saving) {
      edges.push_back({row, column, saving});
      row_max[row] = (std::max)(row_max[row], saving);
      column_max[column] = (std::max)(column_max[column], saving);
    });
  }
  if (strategy == WassersteinWarmStart::global_descending) {
    std::sort(edges.begin(), edges.end(), [](const RankedEdge& left,
                                             const RankedEdge& right) {
      return left.saving > right.saving;
    });
  } else {
    std::stable_sort(edges.begin(), edges.end(), [](const RankedEdge& left,
                                                    const RankedEdge& right) {
      if (left.row != right.row) {
        return left.row < right.row;
      }
      return left.saving > right.saving;
    });
  }
  std::vector<int> row_match(graph.rows, -1);
  std::vector<int> column_match(graph.columns, -1);
  std::vector<double> matched_weight(graph.rows, 0.0);
  Weight total = 0;
  for (const RankedEdge& edge : edges) {
    if (row_match[edge.row] >= 0 || column_match[edge.column] >= 0) {
      continue;
    }
    row_match[edge.row] = static_cast<int>(edge.column);
    column_match[edge.column] = static_cast<int>(edge.row);
    matched_weight[edge.row] = edge.saving;
    total += edge.saving;
    if (stats != nullptr) {
      ++stats->greedy_matches;
    }
  }
  bool reaches_row_upper = true;
  for (std::size_t row = 0; row < graph.rows; ++row) {
    if (row_max[row] > 0.0 && matched_weight[row] != row_max[row]) {
      reaches_row_upper = false;
      break;
    }
  }
  bool reaches_column_upper = true;
  for (std::size_t column = 0; column < graph.columns; ++column) {
    if (column_max[column] == 0.0) {
      continue;
    }
    const int row = column_match[column];
    if (row < 0 || matched_weight[static_cast<std::size_t>(row)] != column_max[column]) {
      reaches_column_upper = false;
      break;
    }
  }
  if (reaches_row_upper || reaches_column_upper) {
    if (stats != nullptr) {
      ++stats->warm_start_certificates;
    }
    MatchingResult result(graph.rows);
    result.saving = total;
    result.row_to_column = std::move(row_match);
    return result;
  }
  return std::nullopt;
}

MatchingResult solve_graph(const WeightedGraph& graph,
                           WassersteinMatcherStrategy requested,
                           WassersteinWarmStart warm_start,
                           WassersteinMetric metric,
                           WassersteinStats* stats,
                           WassersteinWorkspace::Impl* workspace = nullptr) {
  const double density = graph.rows == 0 || graph.columns == 0
                             ? 0.0
                             : static_cast<double>(graph.edge_count) /
                                   static_cast<double>(graph.rows * graph.columns);
  if ((std::max)(graph.rows, graph.columns) <= 32 || density <= 0.05) {
    if (const auto certified = certified_greedy(graph, warm_start, stats)) {
      return *certified;
    }
  }
  WassersteinMatcherStrategy strategy = requested;
  if (strategy == WassersteinMatcherStrategy::adaptive) {
    const bool use_dense =
        (std::max)(graph.rows, graph.columns) <= 24 || density >= 0.15;
    if (!use_dense) {
      strategy = WassersteinMatcherStrategy::sparse_sap;
    } else if (metric == WassersteinMetric::w1_linf &&
               (std::max)(graph.rows, graph.columns) >= 512) {
      strategy = WassersteinMatcherStrategy::dense_sap_row_reduction;
    } else {
      strategy = WassersteinMatcherStrategy::dense_sap;
    }
  }
  switch (strategy) {
    case WassersteinMatcherStrategy::dense_hungarian:
      return dense_assignment(graph, false, stats);
    case WassersteinMatcherStrategy::dense_sap:
      return dense_rectangular_sap(graph, true, stats, workspace, false);
    case WassersteinMatcherStrategy::dense_sap_row_reduction:
      return dense_rectangular_sap(graph, true, stats, workspace, true);
    case WassersteinMatcherStrategy::dense_sap_jv_reduction:
      return dense_rectangular_sap(graph, true, stats, workspace, false, true);
    case WassersteinMatcherStrategy::sparse_sap:
      return sparse_shortest_augmenting_path(graph, stats);
    case WassersteinMatcherStrategy::sparse_sap_arena:
      return sparse_shortest_augmenting_path_arena(graph, stats);
    case WassersteinMatcherStrategy::adaptive:
      break;
  }
  return dense_assignment(graph, false, stats);
}

struct Component {
  std::vector<std::size_t> rows;
  std::vector<std::size_t> columns;
};

std::vector<Component> graph_components(const WeightedGraph& graph,
                                        WassersteinStats* stats) {
  std::vector<std::vector<std::size_t>> reverse(graph.columns);
  for (std::size_t row = 0; row < graph.rows; ++row) {
    graph.for_each_edge(row, [&](std::size_t column, double) {
      reverse[column].push_back(row);
    });
  }
  std::vector<bool> visited_rows(graph.rows, false);
  std::vector<bool> visited_columns(graph.columns, false);
  std::vector<Component> result;
  for (std::size_t start = 0; start < graph.rows; ++start) {
    bool has_edge = false;
    graph.for_each_edge(start, [&](std::size_t, double) { has_edge = true; });
    if (!has_edge || visited_rows[start]) {
      continue;
    }
    Component component;
    std::queue<std::pair<bool, std::size_t>> queue;
    visited_rows[start] = true;
    queue.push({true, start});
    while (!queue.empty()) {
      const auto [is_row, index] = queue.front();
      queue.pop();
      if (is_row) {
        component.rows.push_back(index);
        graph.for_each_edge(index, [&](std::size_t column, double) {
          if (!visited_columns[column]) {
            visited_columns[column] = true;
            queue.push({false, column});
          }
        });
      } else {
        component.columns.push_back(index);
        for (std::size_t row : reverse[index]) {
          if (!visited_rows[row]) {
            visited_rows[row] = true;
            queue.push({true, row});
          }
        }
      }
    }
    if (stats != nullptr) {
      ++stats->component_count;
      stats->largest_component =
          (std::max)(stats->largest_component,
                     static_cast<std::uint64_t>(component.rows.size() +
                                                component.columns.size()));
    }
    result.push_back(std::move(component));
  }
  return result;
}

WeightedGraph component_graph(const WeightedGraph& graph, const Component& component) {
  std::vector<int> column_map(graph.columns, -1);
  for (std::size_t index = 0; index < component.columns.size(); ++index) {
    column_map[component.columns[index]] = static_cast<int>(index);
  }
  CandidateRows candidates(component.rows.size(), component.columns.size());
  for (std::size_t local_row = 0; local_row < component.rows.size(); ++local_row) {
    graph.for_each_edge(component.rows[local_row],
                        [&](std::size_t column, double saving) {
      const int local_column = column_map[column];
      if (local_column >= 0) {
        candidates.edges[local_row].push_back(
            {static_cast<std::uint32_t>(local_column), saving});
        ++candidates.edge_count;
      }
    });
  }
  return build_graph(candidates, WassersteinGraphStrategy::row_vectors, nullptr);
}

MatchingResult exhaustive_component(const WeightedGraph& graph) {
  MatchingResult best(graph.rows);
  std::vector<int> current(graph.rows, -1);
  std::vector<bool> used(graph.columns, false);
  const auto search = [&](auto&& self, std::size_t row, Weight saving) -> void {
    if (row == graph.rows) {
      if (saving > best.saving) {
        best.saving = saving;
        best.row_to_column = current;
      }
      return;
    }
    self(self, row + 1, saving);
    graph.for_each_edge(row, [&](std::size_t column, double edge_saving) {
      if (used[column]) {
        return;
      }
      used[column] = true;
      current[row] = static_cast<int>(column);
      self(self, row + 1, saving + edge_saving);
      current[row] = -1;
      used[column] = false;
    });
  };
  search(search, 0, 0);
  return best;
}

MatchingResult solve_priced_restricted_graph(const WeightedGraph& graph,
                                             WassersteinStats* stats,
                                             SparseDualState& dual) {
  MatchingResult matching(graph.rows);
  dual.row_potential.assign(graph.rows, Weight{0});
  dual.column_potential.assign(graph.columns, Weight{0});
  dual.valid = true;
  const std::vector<Component> components = graph_components(graph, stats);
  for (const Component& component : components) {
    WeightedGraph local = component_graph(graph, component);
    SparseDualState local_dual;
    MatchingResult local_matching;
    if (local.rows + local.columns <= tiny_component_limit) {
      if (stats != nullptr) {
        ++stats->tiny_components;
      }
      local_matching = exhaustive_component(local);
      local_dual.valid =
          build_matching_dual(local, local_matching, local_dual);
    } else {
      local_matching = sparse_shortest_augmenting_path(
          local, stats, true, &local_dual);
    }
    matching.saving += local_matching.saving;
    dual.valid = dual.valid && local_dual.valid;
    for (std::size_t local_row = 0; local_row < local.rows; ++local_row) {
      dual.row_potential[component.rows[local_row]] =
          local_dual.row_potential[local_row];
      const int local_column = local_matching.row_to_column[local_row];
      if (local_column >= 0) {
        matching.row_to_column[component.rows[local_row]] = static_cast<int>(
            component.columns[static_cast<std::size_t>(local_column)]);
      }
    }
    for (std::size_t local_column = 0; local_column < local.columns;
         ++local_column) {
      dual.column_potential[component.columns[local_column]] =
          local_dual.column_potential[local_column];
    }
  }
  return matching;
}

struct DuplicateGroup {
  Point point;
  double half_persistence = 0.0;
  int multiplicity = 0;
};

std::vector<DuplicateGroup> duplicate_groups(const PreparedDiagram& diagram) {
  const auto& points = diagram.finite_points();
  const auto& persistence = diagram.finite_half_persistences();
  const auto& representatives = diagram.finite_duplicate_representatives();
  const auto& multiplicities = diagram.finite_duplicate_multiplicities();
  std::vector<DuplicateGroup> groups;
  groups.reserve(representatives.size());
  for (std::size_t group = 0; group < representatives.size(); ++group) {
    const std::size_t index = representatives[group];
    groups.push_back({points[index], persistence[index],
                      static_cast<int>(multiplicities[group])});
  }
  return groups;
}

std::optional<double> duplicate_compressed_distance(
  const PreparedDiagram& first, const PreparedDiagram& second,
    WassersteinMetric metric, Weight essential,
    WassersteinDuplicateStrategy strategy, WassersteinStats* stats) {
  const auto candidate_start = Clock::now();
  const std::size_t first_group_count =
      first.finite_duplicate_representatives().size();
  const std::size_t second_group_count =
      second.finite_duplicate_representatives().size();
  const std::size_t removed =
      first.finite_points().size() + second.finite_points().size() -
      first_group_count - second_group_count;
  const std::uint64_t group_pairs =
      static_cast<std::uint64_t>(first_group_count) * second_group_count;
  const std::uint64_t original_pairs =
      static_cast<std::uint64_t>(first.finite_points().size()) *
      second.finite_points().size();
  if (stats != nullptr) {
    stats->duplicate_groups += first_group_count + second_group_count;
    stats->duplicate_points_removed += removed;
  }
  if (strategy == WassersteinDuplicateStrategy::adaptive &&
      (removed == 0 || group_pairs > original_pairs / 16)) {
    if (stats != nullptr) {
      stats->candidate_time_ns += elapsed_ns(candidate_start, Clock::now());
    }
    return std::nullopt;
  }
  const std::vector<DuplicateGroup> first_groups = duplicate_groups(first);
  const std::vector<DuplicateGroup> second_groups = duplicate_groups(second);

  const int source = 0;
  const int row_base = 1;
  const int column_base = row_base + static_cast<int>(first_groups.size());
  const int sink = column_base + static_cast<int>(second_groups.size());
  std::vector<std::vector<ResidualEdge>> network(static_cast<std::size_t>(sink + 1));
  for (std::size_t row = 0; row < first_groups.size(); ++row) {
    add_residual_edge(network, source, row_base + static_cast<int>(row), 0,
                      first_groups[row].multiplicity);
  }
  for (std::size_t column = 0; column < second_groups.size(); ++column) {
    add_residual_edge(network, column_base + static_cast<int>(column), sink, 0,
                      second_groups[column].multiplicity);
  }
  std::vector<Weight> maximum_column_saving(second_groups.size(), Weight{0});
  std::uint64_t positive_edges = 0;
  for (std::size_t row = 0; row < first_groups.size(); ++row) {
    for (std::size_t column = 0; column < second_groups.size(); ++column) {
      if (stats != nullptr) {
        ++stats->candidate_pairs;
      }
      const Weight saving =
          raw_saving(first_groups[row].point, second_groups[column].point,
                     first_groups[row].half_persistence,
                     second_groups[column].half_persistence, metric);
      if (saving <= Weight{0}) {
        continue;
      }
      ++positive_edges;
      maximum_column_saving[column] =
          (std::max)(maximum_column_saving[column], saving);
      add_residual_edge(
          network, row_base + static_cast<int>(row),
          column_base + static_cast<int>(column), -saving,
          (std::min)(first_groups[row].multiplicity,
                     second_groups[column].multiplicity));
    }
  }
  if (stats != nullptr) {
    stats->positive_edges += positive_edges;
    stats->pruned_pairs += group_pairs - positive_edges;
    stats->candidate_time_ns += elapsed_ns(candidate_start, Clock::now());
  }

  const auto solver_start = Clock::now();
  const std::size_t node_count = network.size();
  std::vector<Weight> potential(node_count, Weight{0});
  Weight minimum_column_potential = Weight{0};
  for (std::size_t column = 0; column < second_groups.size(); ++column) {
    const Weight value = -maximum_column_saving[column];
    potential[static_cast<std::size_t>(column_base) + column] = value;
    minimum_column_potential = (std::min)(minimum_column_potential, value);
  }
  potential[static_cast<std::size_t>(sink)] = minimum_column_potential;

  while (true) {
    std::vector<Weight> distance(node_count,
                                 std::numeric_limits<Weight>::infinity());
    std::vector<int> previous_node(node_count, -1);
    std::vector<int> previous_edge(node_count, -1);
    using QueueItem = std::pair<Weight, int>;
    std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<>> queue;
    distance[source] = 0;
    queue.push({0, source});
    while (!queue.empty()) {
      const auto [queued_distance, node] = queue.top();
      queue.pop();
      if (queued_distance != distance[static_cast<std::size_t>(node)]) {
        continue;
      }
      for (std::size_t edge_index = 0; edge_index < network[node].size();
           ++edge_index) {
        const ResidualEdge& edge = network[node][edge_index];
        if (edge.capacity == 0) {
          continue;
        }
        const Weight raw_reduced =
            edge.cost + potential[static_cast<std::size_t>(node)] -
            potential[static_cast<std::size_t>(edge.destination)];
        if (!std::isfinite(raw_reduced)) {
          if (stats != nullptr) {
            ++stats->sparse_fallbacks;
          }
          return std::nullopt;
        }
        const Weight reduced = (std::max)(Weight{0}, raw_reduced);
        const Weight candidate =
            distance[static_cast<std::size_t>(node)] + reduced;
        if (distance[static_cast<std::size_t>(edge.destination)] <= candidate) {
          continue;
        }
        distance[static_cast<std::size_t>(edge.destination)] = candidate;
        previous_node[static_cast<std::size_t>(edge.destination)] = node;
        previous_edge[static_cast<std::size_t>(edge.destination)] =
            static_cast<int>(edge_index);
        queue.push({candidate, edge.destination});
      }
    }
    if (previous_node[static_cast<std::size_t>(sink)] < 0) {
      break;
    }
    const Weight path_cost = distance[static_cast<std::size_t>(sink)] -
                             potential[static_cast<std::size_t>(source)] +
                             potential[static_cast<std::size_t>(sink)];
    if (path_cost >= Weight{0}) {
      break;
    }
    for (std::size_t node = 0; node < node_count; ++node) {
      if (std::isfinite(distance[node])) {
        potential[node] += distance[node];
      }
    }
    int flow = (std::numeric_limits<int>::max)();
    for (int node = sink; node != source;
         node = previous_node[static_cast<std::size_t>(node)]) {
      const int parent = previous_node[static_cast<std::size_t>(node)];
      const int edge_index = previous_edge[static_cast<std::size_t>(node)];
      flow = (std::min)(flow,
                        network[parent][static_cast<std::size_t>(edge_index)].capacity);
    }
    for (int node = sink; node != source;
         node = previous_node[static_cast<std::size_t>(node)]) {
      const int parent = previous_node[static_cast<std::size_t>(node)];
      const int edge_index = previous_edge[static_cast<std::size_t>(node)];
      ResidualEdge& edge = network[parent][static_cast<std::size_t>(edge_index)];
      edge.capacity -= flow;
      network[node][static_cast<std::size_t>(edge.reverse)].capacity += flow;
    }
    if (stats != nullptr) {
      ++stats->augmentations;
    }
  }

  Weight powered_distance = essential;
  std::vector<int> row_used(first_groups.size(), 0);
  std::vector<int> column_used(second_groups.size(), 0);
  for (std::size_t row = 0; row < first_groups.size(); ++row) {
    const int row_node = row_base + static_cast<int>(row);
    for (const ResidualEdge& edge : network[static_cast<std::size_t>(row_node)]) {
      if (edge.destination < column_base || edge.destination >= sink) {
        continue;
      }
      const int flow = network[static_cast<std::size_t>(edge.destination)]
                                   [static_cast<std::size_t>(edge.reverse)]
                                       .capacity;
      if (flow == 0) {
        continue;
      }
      const std::size_t column =
          static_cast<std::size_t>(edge.destination - column_base);
      row_used[row] += flow;
      column_used[column] += flow;
      powered_distance += static_cast<Weight>(flow) *
                          cross_power(first_groups[row].point,
                                      second_groups[column].point, metric);
    }
  }
  for (std::size_t row = 0; row < first_groups.size(); ++row) {
    powered_distance +=
        static_cast<Weight>(first_groups[row].multiplicity - row_used[row]) *
        diagonal_power(first_groups[row].half_persistence, metric);
  }
  for (std::size_t column = 0; column < second_groups.size(); ++column) {
    powered_distance +=
        static_cast<Weight>(second_groups[column].multiplicity - column_used[column]) *
        diagonal_power(second_groups[column].half_persistence, metric);
  }
  if (stats != nullptr) {
    stats->solver_time_ns += elapsed_ns(solver_start, Clock::now());
  }
  powered_distance = (std::max)(Weight{0}, powered_distance);
  return metric == WassersteinMetric::w1_linf
             ? static_cast<double>(powered_distance)
             : static_cast<double>(std::sqrt(powered_distance));
}

Weight essential_power(const PreparedDiagram& first, const PreparedDiagram& second,
                       WassersteinMetric metric) {
  if (first.fully_infinite_count() != second.fully_infinite_count() ||
      first.positive_infinite_births().size() !=
          second.positive_infinite_births().size() ||
      first.negative_infinite_deaths().size() !=
          second.negative_infinite_deaths().size()) {
    return std::numeric_limits<Weight>::infinity();
  }
  const auto accumulate = [&](const std::vector<double>& left,
                              const std::vector<double>& right) {
    Weight subtotal = 0;
    for (std::size_t index = 0; index < left.size(); ++index) {
      const Weight delta =
          std::fabs(static_cast<Weight>(left[index]) - right[index]);
      subtotal += metric == WassersteinMetric::w1_linf ? delta : delta * delta;
    }
    return subtotal;
  };
  return accumulate(first.positive_infinite_births(),
                     second.positive_infinite_births()) +
         accumulate(first.negative_infinite_deaths(),
                     second.negative_infinite_deaths());
}

double distance_from_matching(const PreparedDiagram& first,
                              const PreparedDiagram& second,
                              WassersteinMetric metric, Weight essential,
                              const MatchingResult& matching) {
  Weight powered_distance = essential;
  std::vector<bool> matched_columns(second.finite_points().size(), false);
  for (std::size_t row = 0; row < first.finite_points().size(); ++row) {
    const int column = matching.row_to_column[row];
    if (column < 0) {
      powered_distance +=
          diagonal_power(first.finite_half_persistences()[row], metric);
      continue;
    }
    const std::size_t matched_column = static_cast<std::size_t>(column);
    matched_columns[matched_column] = true;
    powered_distance += cross_power(first.finite_points()[row],
                                    second.finite_points()[matched_column], metric);
  }
  for (std::size_t column = 0; column < second.finite_points().size(); ++column) {
    if (!matched_columns[column]) {
      powered_distance +=
          diagonal_power(second.finite_half_persistences()[column], metric);
    }
  }
  powered_distance = (std::max)(Weight{0}, powered_distance);
  return metric == WassersteinMetric::w1_linf
             ? static_cast<double>(powered_distance)
             : static_cast<double>(std::sqrt(powered_distance));
}

}  // namespace

double wasserstein_distance_impl(const PreparedDiagram& first,
                                 const PreparedDiagram& second,
                                 const WassersteinConfig& config,
                                 WassersteinStats* stats,
                                 WassersteinWorkspace::Impl* workspace) {
  const std::size_t rows = first.finite_points().size();
  const std::size_t columns = second.finite_points().size();
  const std::uint64_t possible_pairs =
      static_cast<std::uint64_t>(rows) * columns;
  const std::uint64_t positive_before = stats == nullptr ? 0 : stats->positive_edges;
  if (stats != nullptr) {
    stats->possible_pairs += possible_pairs;
  }

  const Weight essential = essential_power(first, second, config.metric);
  if (std::isinf(essential)) {
    return std::numeric_limits<double>::infinity();
  }
  if (config.duplicates != WassersteinDuplicateStrategy::none) {
    if (const auto compressed = duplicate_compressed_distance(
            first, second, config.metric, essential, config.duplicates, stats)) {
      return *compressed;
    }
  }
  if (config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_full_scan ||
      config.candidates == WassersteinCandidateStrategy::topk_pricing_simd ||
      config.candidates == WassersteinCandidateStrategy::topk_pricing_sweep ||
      config.candidates == WassersteinCandidateStrategy::topk_pricing_kdtree ||
      config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_kdtree_persistent ||
      config.candidates == WassersteinCandidateStrategy::topk_pricing_adaptive ||
      config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_sweep_incremental ||
      config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_sweep_persistent ||
      config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_sweep_dynamic ||
      config.candidates == WassersteinCandidateStrategy::
                               topk_pricing_sweep_dynamic_batched) {
    CandidateRows local_candidates(rows, columns);
    CandidateRows& candidates =
        workspace == nullptr ? local_candidates : workspace->candidates;
    WeightedGraph local_graph;
    WeightedGraph& graph = workspace == nullptr ? local_graph : workspace->graph;
    MatchingResult matching;
    if (config.candidates ==
        WassersteinCandidateStrategy::topk_pricing_sweep_persistent) {
      matching = priced_topk_persistent_matching(first, second, config.metric,
                                                 config.top_k, stats, candidates,
                                                 graph);
    } else if (config.candidates ==
                   WassersteinCandidateStrategy::topk_pricing_sweep_dynamic ||
               config.candidates == WassersteinCandidateStrategy::
                                        topk_pricing_sweep_dynamic_batched) {
      matching = priced_topk_dynamic_matching(
          first, second, config.metric, config.top_k,
          config.candidates == WassersteinCandidateStrategy::
                                   topk_pricing_sweep_dynamic_batched,
          stats, candidates, graph);
    } else {
      matching = priced_topk_matching(
                  first, second, config.metric, config.top_k,
                  config.candidates ==
                          WassersteinCandidateStrategy::topk_pricing_full_scan ||
                      config.candidates ==
                          WassersteinCandidateStrategy::topk_pricing_simd,
                  config.candidates ==
                          WassersteinCandidateStrategy::topk_pricing_kdtree ||
                      config.candidates == WassersteinCandidateStrategy::
                                               topk_pricing_kdtree_persistent,
                  config.candidates ==
                      WassersteinCandidateStrategy::topk_pricing_simd,
                  config.candidates == WassersteinCandidateStrategy::
                                           topk_pricing_kdtree_persistent,
                  config.candidates ==
                      WassersteinCandidateStrategy::topk_pricing_adaptive,
                  config.candidates == WassersteinCandidateStrategy::
                                           topk_pricing_sweep_incremental,
                   stats, candidates, graph);
    }
    if (stats != nullptr) {
      const std::uint64_t call_positive = stats->positive_edges - positive_before;
      stats->pruned_pairs += possible_pairs - call_positive;
    }
    return distance_from_matching(first, second, config.metric, essential, matching);
  }
  const auto candidate_start = Clock::now();
  CandidateRows local_candidates(rows, columns);
  CandidateRows& candidates =
      workspace == nullptr ? local_candidates : workspace->candidates;
  generate_candidates(first, second, config, stats, candidates);
  const auto candidate_stop = Clock::now();
  if (stats != nullptr) {
    stats->candidate_time_ns += elapsed_ns(candidate_start, candidate_stop);
    const std::uint64_t call_positive = stats->positive_edges - positive_before;
    stats->pruned_pairs += possible_pairs - call_positive;
  }

  const auto graph_start = Clock::now();
  WeightedGraph local_graph;
  WeightedGraph& graph = workspace == nullptr ? local_graph : workspace->graph;
  build_graph_into(candidates, config.graph, stats, &first, &second, config.metric,
                   graph);
  const auto graph_stop = Clock::now();
  if (stats != nullptr) {
    stats->graph_time_ns += elapsed_ns(graph_start, graph_stop);
  }

  WassersteinComponentStrategy component_strategy = config.components;
  if (component_strategy == WassersteinComponentStrategy::adaptive) {
    const double density = possible_pairs == 0
                               ? 0.0
                               : static_cast<double>(graph.edge_count) /
                                     static_cast<double>(possible_pairs);
    if (density >= 0.15) {
      component_strategy = WassersteinComponentStrategy::none;
    }
  }

  MatchingResult matching(rows);
  if (component_strategy == WassersteinComponentStrategy::none) {
    const auto solver_start = Clock::now();
    matching =
        solve_graph(graph, config.matcher, config.warm_start, config.metric, stats,
                    workspace);
    const auto solver_stop = Clock::now();
    if (stats != nullptr) {
      stats->solver_time_ns += elapsed_ns(solver_start, solver_stop);
    }
  } else {
    const auto component_start = Clock::now();
    const std::vector<Component> components = graph_components(graph, stats);
    const auto component_stop = Clock::now();
    if (stats != nullptr) {
      stats->component_time_ns += elapsed_ns(component_start, component_stop);
    }
    const auto solver_start = Clock::now();
    if ((component_strategy == WassersteinComponentStrategy::parallel_dense ||
         component_strategy == WassersteinComponentStrategy::parallel_sparse) &&
        components.size() >= 2 && rows + columns >= 256) {
      struct ParallelComponentResult {
        MatchingResult matching;
        WassersteinStats stats;
      };
      std::vector<ParallelComponentResult> results(components.size());
      std::atomic<std::size_t> next_component{0};
      const std::size_t hardware =
          (std::max)(std::size_t{2}, static_cast<std::size_t>(
                                          std::thread::hardware_concurrency()));
      const std::size_t thread_count =
          (std::min)({std::size_t{8}, hardware, components.size()});
      std::vector<std::thread> workers;
      workers.reserve(thread_count);
      for (std::size_t worker = 0; worker < thread_count; ++worker) {
        workers.emplace_back([&] {
          while (true) {
            const std::size_t index = next_component.fetch_add(1);
            if (index >= components.size()) {
              break;
            }
            const Component& component = components[index];
            WeightedGraph local = component_graph(graph, component);
            ParallelComponentResult& result = results[index];
            if (local.rows + local.columns <= tiny_component_limit) {
              ++result.stats.tiny_components;
              result.matching = exhaustive_component(local);
            } else if (component_strategy ==
                       WassersteinComponentStrategy::parallel_dense) {
              result.matching = dense_rectangular_sap(
                  local, true, &result.stats, nullptr, false);
            } else {
              result.matching =
                  sparse_shortest_augmenting_path(local, &result.stats);
            }
          }
        });
      }
      for (std::thread& worker : workers) {
        worker.join();
      }
      for (std::size_t index = 0; index < components.size(); ++index) {
        const Component& component = components[index];
        const ParallelComponentResult& result = results[index];
        matching.saving += result.matching.saving;
        for (std::size_t local_row = 0; local_row < component.rows.size();
             ++local_row) {
          const int local_column = result.matching.row_to_column[local_row];
          if (local_column >= 0) {
            matching.row_to_column[component.rows[local_row]] = static_cast<int>(
                component.columns[static_cast<std::size_t>(local_column)]);
          }
        }
        if (stats != nullptr) {
          stats->active_rows += result.stats.active_rows;
          stats->active_columns += result.stats.active_columns;
          stats->augmentations += result.stats.augmentations;
          stats->sparse_fallbacks += result.stats.sparse_fallbacks;
          stats->tiny_components += result.stats.tiny_components;
        }
      }
      const auto solver_stop = Clock::now();
      if (stats != nullptr) {
        stats->solver_time_ns += elapsed_ns(solver_start, solver_stop);
      }
      return distance_from_matching(first, second, config.metric, essential,
                                    matching);
    }
    if (component_strategy == WassersteinComponentStrategy::parallel_dense) {
      component_strategy = WassersteinComponentStrategy::dense;
    } else if (component_strategy ==
               WassersteinComponentStrategy::parallel_sparse) {
      component_strategy = WassersteinComponentStrategy::tiny_sparse;
    }
    // Reuse the already measured decomposition through a local loop rather than
    // calling solve_components, which would discover the same components again.
    for (const Component& component : components) {
      if ((component_strategy == WassersteinComponentStrategy::tiny_sparse ||
           component_strategy == WassersteinComponentStrategy::adaptive) &&
          (component.rows.size() == 1 || component.columns.size() == 1)) {
        Weight best = 0;
        std::size_t best_row = graph.rows;
        std::size_t best_column = graph.columns;
        if (component.rows.size() == 1) {
          const std::size_t target_row = component.rows.front();
          graph.for_each_edge(target_row,
                              [&](std::size_t column, double saving) {
            if (static_cast<Weight>(saving) > best) {
              best = saving;
              best_row = target_row;
              best_column = column;
            }
          });
        } else {
          const std::size_t target_column = component.columns.front();
          for (std::size_t row : component.rows) {
            graph.for_each_edge(row, [&](std::size_t column, double saving) {
              if (column == target_column && static_cast<Weight>(saving) > best) {
                best = saving;
                best_row = row;
                best_column = column;
              }
            });
          }
        }
        if (stats != nullptr) {
          ++stats->tiny_components;
        }
        if (best_row != graph.rows) {
          matching.saving += best;
          matching.row_to_column[best_row] = static_cast<int>(best_column);
        }
        continue;
      }
      WeightedGraph local = component_graph(graph, component);
      const std::size_t vertex_count = local.rows + local.columns;
      const double density = local.rows == 0 || local.columns == 0
                                 ? 0.0
                                 : static_cast<double>(local.edge_count) /
                                       static_cast<double>(local.rows * local.columns);
      if ((component_strategy == WassersteinComponentStrategy::tiny_sparse ||
           component_strategy == WassersteinComponentStrategy::adaptive) &&
          vertex_count <= tiny_component_limit) {
        if (stats != nullptr) {
          ++stats->tiny_components;
        }
        MatchingResult local_matching = exhaustive_component(local);
        matching.saving += local_matching.saving;
        for (std::size_t local_row = 0; local_row < local.rows; ++local_row) {
          const int local_column = local_matching.row_to_column[local_row];
          if (local_column >= 0) {
            matching.row_to_column[component.rows[local_row]] = static_cast<int>(
                component.columns[static_cast<std::size_t>(local_column)]);
          }
        }
        continue;
      }
      WassersteinMatcherStrategy local_matcher = config.matcher;
      if (component_strategy == WassersteinComponentStrategy::dense) {
        local_matcher = WassersteinMatcherStrategy::dense_sap;
      } else if (component_strategy == WassersteinComponentStrategy::sparse ||
                 component_strategy == WassersteinComponentStrategy::tiny_sparse ||
                 component_strategy ==
                     WassersteinComponentStrategy::parallel_sparse) {
        local_matcher = WassersteinMatcherStrategy::sparse_sap;
      } else if (component_strategy == WassersteinComponentStrategy::adaptive) {
        local_matcher = (std::max)(local.rows, local.columns) <= 24 || density >= 0.20
                            ? WassersteinMatcherStrategy::dense_sap
                            : WassersteinMatcherStrategy::sparse_sap;
      }
      MatchingResult local_matching =
          solve_graph(local, local_matcher,
                      local_matcher == WassersteinMatcherStrategy::dense_sap
                          ? WassersteinWarmStart::none
                          : config.warm_start,
                      config.metric, stats, workspace);
      matching.saving += local_matching.saving;
      for (std::size_t local_row = 0; local_row < local.rows; ++local_row) {
        const int local_column = local_matching.row_to_column[local_row];
        if (local_column >= 0) {
          matching.row_to_column[component.rows[local_row]] = static_cast<int>(
              component.columns[static_cast<std::size_t>(local_column)]);
        }
      }
    }
    const auto solver_stop = Clock::now();
    if (stats != nullptr) {
      stats->solver_time_ns += elapsed_ns(solver_start, solver_stop);
    }
  }

  return distance_from_matching(first, second, config.metric, essential, matching);
}

double wasserstein_distance(const PreparedDiagram& first,
                            const PreparedDiagram& second,
                            const WassersteinConfig& config,
                            WassersteinStats* stats) {
  return wasserstein_distance_impl(first, second, config, stats, nullptr);
}

double wasserstein_distance(const PreparedDiagram& first,
                            const PreparedDiagram& second,
                            WassersteinWorkspace& workspace,
                            const WassersteinConfig& config,
                            WassersteinStats* stats) {
  return wasserstein_distance_impl(first, second, config, stats,
                                   workspace.impl_.get());
}

double wasserstein_distance(const Diagram& first, const Diagram& second,
                            const WassersteinConfig& config, WassersteinStats* stats) {
  const auto start = Clock::now();
  const PreparedDiagram prepared_first(first);
  const PreparedDiagram prepared_second(second);
  const auto stop = Clock::now();
  if (stats != nullptr) {
    stats->prepare_time_ns += elapsed_ns(start, stop);
  }
  return wasserstein_distance(prepared_first, prepared_second, config, stats);
}

WassersteinWorkspace::WassersteinWorkspace() : impl_(std::make_unique<Impl>()) {}
WassersteinWorkspace::~WassersteinWorkspace() = default;
WassersteinWorkspace::WassersteinWorkspace(WassersteinWorkspace&&) noexcept = default;
WassersteinWorkspace& WassersteinWorkspace::operator=(
    WassersteinWorkspace&&) noexcept = default;

void wasserstein_distances(const PreparedDiagram& query,
                           std::span<const PreparedDiagram> diagrams,
                           std::span<double> output,
                           const WassersteinConfig& config,
                           WassersteinStats* stats) {
  WassersteinWorkspace workspace;
  wasserstein_distances(query, diagrams, output, workspace, config, stats);
}

void wasserstein_distances(const PreparedDiagram& query,
                           std::span<const PreparedDiagram> diagrams,
                           std::span<double> output,
                           WassersteinWorkspace& workspace,
                           const WassersteinConfig& config,
                           WassersteinStats* stats) {
  if (output.size() != diagrams.size()) {
    throw std::invalid_argument("output size must match diagram count");
  }
  for (std::size_t index = 0; index < diagrams.size(); ++index) {
    output[index] =
        wasserstein_distance(query, diagrams[index], workspace, config, stats);
  }
}

std::vector<double> wasserstein_distances(
    const PreparedDiagram& query, std::span<const PreparedDiagram> diagrams,
    const WassersteinConfig& config, WassersteinStats* stats) {
  std::vector<double> result(diagrams.size());
  wasserstein_distances(query, diagrams, result, config, stats);
  return result;
}

const char* to_string(WassersteinMetric metric) noexcept {
  switch (metric) {
    case WassersteinMetric::w1_linf:
      return "w1_linf";
    case WassersteinMetric::w2_l2:
      return "w2_l2";
  }
  return "unknown";
}

const char* to_string(WassersteinCandidateStrategy strategy) noexcept {
  switch (strategy) {
    case WassersteinCandidateStrategy::dense_scalar:
      return "dense_scalar";
    case WassersteinCandidateStrategy::dense_blocked:
      return "dense_blocked";
    case WassersteinCandidateStrategy::dense_avx2:
      return "dense_avx2";
    case WassersteinCandidateStrategy::dense_parallel:
      return "dense_parallel";
    case WassersteinCandidateStrategy::sweep_binary:
      return "sweep_binary";
    case WassersteinCandidateStrategy::sweep_two_pointer:
      return "sweep_two_pointer";
    case WassersteinCandidateStrategy::topk_pricing_full_scan:
      return "topk_pricing_full_scan";
    case WassersteinCandidateStrategy::topk_pricing_simd:
      return "topk_pricing_simd";
    case WassersteinCandidateStrategy::topk_pricing_sweep:
      return "topk_pricing_sweep";
    case WassersteinCandidateStrategy::topk_pricing_kdtree:
      return "topk_pricing_kdtree";
    case WassersteinCandidateStrategy::topk_pricing_kdtree_persistent:
      return "topk_pricing_kdtree_persistent";
    case WassersteinCandidateStrategy::topk_pricing_adaptive:
      return "topk_pricing_adaptive";
    case WassersteinCandidateStrategy::topk_pricing_sweep_incremental:
      return "topk_pricing_sweep_incremental";
    case WassersteinCandidateStrategy::topk_pricing_sweep_persistent:
      return "topk_pricing_sweep_persistent";
    case WassersteinCandidateStrategy::topk_pricing_sweep_dynamic:
      return "topk_pricing_sweep_dynamic";
    case WassersteinCandidateStrategy::topk_pricing_sweep_dynamic_batched:
      return "topk_pricing_sweep_dynamic_batched";
    case WassersteinCandidateStrategy::adaptive:
      return "adaptive";
  }
  return "unknown";
}

const char* to_string(WassersteinGraphStrategy strategy) noexcept {
  switch (strategy) {
    case WassersteinGraphStrategy::dense_matrix:
      return "dense_matrix";
    case WassersteinGraphStrategy::csr:
      return "csr";
    case WassersteinGraphStrategy::row_vectors:
      return "row_vectors";
    case WassersteinGraphStrategy::fixed_degree_4:
      return "fixed_degree_4";
    case WassersteinGraphStrategy::fixed_degree_8:
      return "fixed_degree_8";
    case WassersteinGraphStrategy::fixed_degree_16:
      return "fixed_degree_16";
    case WassersteinGraphStrategy::fixed_degree_32:
      return "fixed_degree_32";
    case WassersteinGraphStrategy::bitmask_lazy:
      return "bitmask_lazy";
    case WassersteinGraphStrategy::block_sparse_16:
      return "block_sparse_16";
    case WassersteinGraphStrategy::block_sparse_32:
      return "block_sparse_32";
    case WassersteinGraphStrategy::adaptive:
      return "adaptive";
  }
  return "unknown";
}

const char* to_string(WassersteinMatcherStrategy strategy) noexcept {
  switch (strategy) {
    case WassersteinMatcherStrategy::dense_hungarian:
      return "dense_hungarian";
    case WassersteinMatcherStrategy::dense_sap:
      return "dense_sap";
    case WassersteinMatcherStrategy::dense_sap_row_reduction:
      return "dense_sap_row_reduction";
    case WassersteinMatcherStrategy::dense_sap_jv_reduction:
      return "dense_sap_jv_reduction";
    case WassersteinMatcherStrategy::sparse_sap:
      return "sparse_sap";
    case WassersteinMatcherStrategy::sparse_sap_arena:
      return "sparse_sap_arena";
    case WassersteinMatcherStrategy::adaptive:
      return "adaptive";
  }
  return "unknown";
}

const char* to_string(WassersteinComponentStrategy strategy) noexcept {
  switch (strategy) {
    case WassersteinComponentStrategy::none:
      return "none";
    case WassersteinComponentStrategy::dense:
      return "dense";
    case WassersteinComponentStrategy::sparse:
      return "sparse";
    case WassersteinComponentStrategy::tiny_sparse:
      return "tiny_sparse";
    case WassersteinComponentStrategy::parallel_dense:
      return "parallel_dense";
    case WassersteinComponentStrategy::parallel_sparse:
      return "parallel_sparse";
    case WassersteinComponentStrategy::adaptive:
      return "adaptive";
  }
  return "unknown";
}

const char* to_string(WassersteinWarmStart strategy) noexcept {
  switch (strategy) {
    case WassersteinWarmStart::none:
      return "none";
    case WassersteinWarmStart::row_max:
      return "row_max";
    case WassersteinWarmStart::global_descending:
      return "global_descending";
  }
  return "unknown";
}

const char* to_string(WassersteinDuplicateStrategy strategy) noexcept {
  switch (strategy) {
    case WassersteinDuplicateStrategy::none:
      return "none";
    case WassersteinDuplicateStrategy::exact:
      return "exact";
    case WassersteinDuplicateStrategy::adaptive:
      return "adaptive";
  }
  return "unknown";
}

}  // namespace bottleneck
