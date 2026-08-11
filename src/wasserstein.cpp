#include <bottleneck/wasserstein.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
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

bool cpu_has_avx2() noexcept {
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

struct CandidateRows {
  CandidateRows(std::size_t row_count, std::size_t column_count)
      : rows(row_count), columns(column_count), edges(row_count) {}

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

enum class StoredGraphKind { dense, csr, rows };

struct WeightedGraph {
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
    for (const Edge& edge : row_edges[row]) {
      function(static_cast<std::size_t>(edge.column), edge.saving);
    }
  }

  [[nodiscard]] double saving(std::size_t row, std::size_t column) const {
    if (kind == StoredGraphKind::dense) {
      return dense[row * columns + column];
    }
    double result = 0.0;
    for_each_edge(row, [&](std::size_t candidate, double value) {
      if (candidate == column) {
        result = value;
      }
    });
    return result;
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
};

Weight diagonal_power(double half_persistence, WassersteinMetric metric) {
  const Weight persistence = static_cast<Weight>(half_persistence);
  return metric == WassersteinMetric::w1_linf
             ? persistence
             : Weight{2} * persistence * persistence;
}

Weight raw_saving(const Point& first, const Point& second,
                  double first_half_persistence, double second_half_persistence,
                  WassersteinMetric metric) {
  const Weight first_diagonal = diagonal_power(first_half_persistence, metric);
  const Weight second_diagonal = diagonal_power(second_half_persistence, metric);
  const Weight birth_delta =
      std::fabs(static_cast<Weight>(first.birth) - second.birth);
  const Weight death_delta =
      std::fabs(static_cast<Weight>(first.death) - second.death);
  if (metric == WassersteinMetric::w1_linf) {
    return first_diagonal + second_diagonal - (std::max)(birth_delta, death_delta);
  }
  return first_diagonal + second_diagonal -
         (birth_delta * birth_delta + death_delta * death_delta);
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

CandidateRows dense_scalar_candidates(const PreparedDiagram& first,
                                      const PreparedDiagram& second,
                                      WassersteinMetric metric,
                                      WassersteinStats* stats) {
  const auto& first_points = first.finite_points();
  const auto& second_points = second.finite_points();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_v = second.finite_half_persistences();
  CandidateRows result(first_points.size(), second_points.size());
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
  return result;
}

CandidateRows dense_blocked_candidates(const PreparedDiagram& first,
                                       const PreparedDiagram& second,
                                       WassersteinMetric metric,
                                       WassersteinStats* stats) {
  const auto& first_u = first.finite_midpoints();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_u = second.finite_midpoints();
  const auto& second_v = second.finite_half_persistences();
  CandidateRows result(first_u.size(), second_u.size());
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
  return result;
}

CandidateRows dense_avx2_candidates(const PreparedDiagram& first,
                                    const PreparedDiagram& second,
                                    WassersteinMetric metric,
                                    WassersteinStats* stats) {
#if defined(BOTTLENECK_HAVE_WASSERSTEIN_AVX2)
  if (cpu_has_avx2()) {
    const auto& first_u = first.finite_midpoints();
    const auto& first_v = first.finite_half_persistences();
    const auto& second_u = second.finite_midpoints();
    const auto& second_v = second.finite_half_persistences();
    CandidateRows result(first_u.size(), second_u.size());
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
    return result;
  }
#endif
  return dense_blocked_candidates(first, second, metric, stats);
}

CandidateRows sweep_binary_candidates(const PreparedDiagram& first,
                                      const PreparedDiagram& second,
                                      WassersteinMetric metric,
                                      WassersteinStats* stats) {
  const auto& first_u = first.finite_midpoints();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_u = second.finite_midpoints();
  const auto& second_v = second.finite_half_persistences();
  const auto& sorted_u = second.sorted_finite_midpoints();
  const auto& order = second.finite_midpoint_order();
  CandidateRows result(first_u.size(), second_u.size());
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
  return result;
}

CandidateRows sweep_two_pointer_candidates(const PreparedDiagram& first,
                                           const PreparedDiagram& second,
                                           WassersteinMetric metric,
                                           WassersteinStats* stats) {
  if (metric == WassersteinMetric::w2_l2) {
    return sweep_binary_candidates(first, second, metric, stats);
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

  CandidateRows result(first_u.size(), second_u.size());
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
  return result;
}

double sampled_positive_density(const PreparedDiagram& first,
                                const PreparedDiagram& second,
                                WassersteinMetric metric) {
  const std::size_t sample_rows = (std::min)(std::size_t{8}, first.finite_points().size());
  if (sample_rows == 0 || second.finite_points().empty()) {
    return 0.0;
  }
  std::size_t positive = 0;
  std::size_t sampled = 0;
  const auto& first_points = first.finite_points();
  const auto& second_points = second.finite_points();
  const auto& first_v = first.finite_half_persistences();
  const auto& second_v = second.finite_half_persistences();
  for (std::size_t sample = 0; sample < sample_rows; ++sample) {
    const std::size_t row = sample * first_points.size() / sample_rows;
    for (std::size_t column = 0; column < second_points.size(); ++column) {
      positive += static_cast<std::size_t>(
          raw_saving(first_points[row], second_points[column], first_v[row],
                     second_v[column], metric) > Weight{0});
      ++sampled;
    }
  }
  return static_cast<double>(positive) / static_cast<double>(sampled);
}

CandidateRows generate_candidates(const PreparedDiagram& first,
                                  const PreparedDiagram& second,
                                  WassersteinConfig config,
                                  WassersteinStats* stats) {
  WassersteinCandidateStrategy strategy = config.candidates;
  if (strategy == WassersteinCandidateStrategy::adaptive) {
    const std::size_t pairs =
        first.finite_points().size() * second.finite_points().size();
    if (pairs <= 1024) {
      strategy = WassersteinCandidateStrategy::dense_scalar;
    } else if (sampled_positive_density(first, second, config.metric) >= 0.35) {
      strategy = cpu_has_avx2() ? WassersteinCandidateStrategy::dense_avx2
                                : WassersteinCandidateStrategy::dense_blocked;
    } else {
      strategy = WassersteinCandidateStrategy::sweep_binary;
    }
  }
  switch (strategy) {
    case WassersteinCandidateStrategy::dense_scalar:
      return dense_scalar_candidates(first, second, config.metric, stats);
    case WassersteinCandidateStrategy::dense_blocked:
      return dense_blocked_candidates(first, second, config.metric, stats);
    case WassersteinCandidateStrategy::dense_avx2:
      return dense_avx2_candidates(first, second, config.metric, stats);
    case WassersteinCandidateStrategy::sweep_binary:
      return sweep_binary_candidates(first, second, config.metric, stats);
    case WassersteinCandidateStrategy::sweep_two_pointer:
      return sweep_two_pointer_candidates(first, second, config.metric, stats);
    case WassersteinCandidateStrategy::adaptive:
      break;
  }
  return dense_scalar_candidates(first, second, config.metric, stats);
}

WeightedGraph build_graph(const CandidateRows& candidates,
                          WassersteinGraphStrategy requested,
                          WassersteinStats* stats) {
  WassersteinGraphStrategy strategy = requested;
  const double density = candidates.rows == 0 || candidates.columns == 0
                             ? 0.0
                             : static_cast<double>(candidates.edge_count) /
                                   static_cast<double>(candidates.rows * candidates.columns);
  if (strategy == WassersteinGraphStrategy::adaptive) {
    strategy = candidates.rows * candidates.columns <= 1024 || density >= 0.25
                   ? WassersteinGraphStrategy::dense_matrix
                   : WassersteinGraphStrategy::csr;
  }

  WeightedGraph graph;
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
  } else {
    graph.kind = StoredGraphKind::rows;
    graph.row_edges = candidates.edges;
    if (stats != nullptr) {
      stats->graph_bytes += graph.row_edges.size() * sizeof(std::vector<Edge>) +
                            graph.edge_count * sizeof(Edge);
    }
  }
  return graph;
}

struct ActiveVertices {
  std::vector<std::size_t> rows;
  std::vector<std::size_t> columns;
};

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

Weight dense_assignment(const WeightedGraph& graph, bool compact,
                        WassersteinStats* stats) {
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
    return Weight{0};
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
          saving = graph.saving(active.rows[current_row - 1],
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

  Weight result = 0;
  for (std::size_t column = 1; column <= active.columns.size(); ++column) {
    const std::size_t row = matched_row[column];
    if (row != 0 && row <= active.rows.size()) {
      result += graph.saving(active.rows[row - 1], active.columns[column - 1]);
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
                       int destination, Weight cost) {
  const int source_reverse = static_cast<int>(network[destination].size());
  const int destination_reverse = static_cast<int>(network[source].size());
  network[source].push_back({destination, source_reverse, 1, cost});
  network[destination].push_back({source, destination_reverse, 0, -cost});
}

Weight sparse_shortest_augmenting_path(const WeightedGraph& graph,
                                       WassersteinStats* stats) {
  const ActiveVertices active = active_vertices(graph);
  if (stats != nullptr) {
    stats->active_rows += active.rows.size();
    stats->active_columns += active.columns.size();
  }
  if (active.rows.empty() || active.columns.empty()) {
    return Weight{0};
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
  while (true) {
    std::vector<Weight> distance(node_count,
                                 std::numeric_limits<Weight>::infinity());
    std::vector<int> previous_node(node_count, -1);
    std::vector<int> previous_edge(node_count, -1);
    std::vector<bool> queued(node_count, false);
    std::deque<int> queue;
    distance[source] = 0;
    queue.push_back(source);
    queued[source] = true;
    while (!queue.empty()) {
      const int node = queue.front();
      queue.pop_front();
      queued[static_cast<std::size_t>(node)] = false;
      for (std::size_t edge_index = 0; edge_index < network[node].size(); ++edge_index) {
        const ResidualEdge& edge = network[node][edge_index];
        if (edge.capacity == 0 ||
            distance[static_cast<std::size_t>(edge.destination)] <=
                distance[static_cast<std::size_t>(node)] + edge.cost) {
          continue;
        }
        distance[static_cast<std::size_t>(edge.destination)] =
            distance[static_cast<std::size_t>(node)] + edge.cost;
        previous_node[static_cast<std::size_t>(edge.destination)] = node;
        previous_edge[static_cast<std::size_t>(edge.destination)] =
            static_cast<int>(edge_index);
        if (!queued[static_cast<std::size_t>(edge.destination)]) {
          queued[static_cast<std::size_t>(edge.destination)] = true;
          queue.push_back(edge.destination);
        }
      }
    }
    if (previous_node[static_cast<std::size_t>(sink)] < 0 ||
        distance[static_cast<std::size_t>(sink)] >= Weight{0}) {
      break;
    }
    for (int node = sink; node != source;
         node = previous_node[static_cast<std::size_t>(node)]) {
      const int parent = previous_node[static_cast<std::size_t>(node)];
      const int edge_index = previous_edge[static_cast<std::size_t>(node)];
      ResidualEdge& edge = network[parent][static_cast<std::size_t>(edge_index)];
      --edge.capacity;
      ++network[node][static_cast<std::size_t>(edge.reverse)].capacity;
    }
    total_saving -= distance[static_cast<std::size_t>(sink)];
    if (stats != nullptr) {
      ++stats->augmentations;
    }
  }
  return total_saving;
}

std::optional<Weight> certified_greedy(const WeightedGraph& graph,
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
    return total;
  }
  return std::nullopt;
}

Weight solve_graph(const WeightedGraph& graph, WassersteinMatcherStrategy requested,
                   WassersteinWarmStart warm_start, WassersteinStats* stats) {
  if (const auto certified = certified_greedy(graph, warm_start, stats)) {
    return *certified;
  }
  WassersteinMatcherStrategy strategy = requested;
  if (strategy == WassersteinMatcherStrategy::adaptive) {
    const double density = graph.rows == 0 || graph.columns == 0
                               ? 0.0
                               : static_cast<double>(graph.edge_count) /
                                     static_cast<double>(graph.rows * graph.columns);
    strategy = (std::max)(graph.rows, graph.columns) <= 24 || density >= 0.20
                   ? WassersteinMatcherStrategy::dense_sap
                   : WassersteinMatcherStrategy::sparse_sap;
  }
  switch (strategy) {
    case WassersteinMatcherStrategy::dense_hungarian:
      return dense_assignment(graph, false, stats);
    case WassersteinMatcherStrategy::dense_sap:
      return dense_assignment(graph, true, stats);
    case WassersteinMatcherStrategy::sparse_sap:
      return sparse_shortest_augmenting_path(graph, stats);
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

Weight exhaustive_component(const WeightedGraph& graph) {
  Weight best = 0;
  std::vector<bool> used(graph.columns, false);
  const auto search = [&](auto&& self, std::size_t row, Weight saving) -> void {
    if (row == graph.rows) {
      best = (std::max)(best, saving);
      return;
    }
    self(self, row + 1, saving);
    graph.for_each_edge(row, [&](std::size_t column, double edge_saving) {
      if (used[column]) {
        return;
      }
      used[column] = true;
      self(self, row + 1, saving + edge_saving);
      used[column] = false;
    });
  };
  search(search, 0, 0);
  return best;
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

}  // namespace

double wasserstein_distance(const PreparedDiagram& first, const PreparedDiagram& second,
                            const WassersteinConfig& config, WassersteinStats* stats) {
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
  Weight diagonal_baseline = essential;
  for (double persistence : first.finite_half_persistences()) {
    diagonal_baseline += diagonal_power(persistence, config.metric);
  }
  for (double persistence : second.finite_half_persistences()) {
    diagonal_baseline += diagonal_power(persistence, config.metric);
  }

  const auto candidate_start = Clock::now();
  const CandidateRows candidates = generate_candidates(first, second, config, stats);
  const auto candidate_stop = Clock::now();
  if (stats != nullptr) {
    stats->candidate_time_ns += elapsed_ns(candidate_start, candidate_stop);
    const std::uint64_t call_positive = stats->positive_edges - positive_before;
    stats->pruned_pairs += possible_pairs - call_positive;
  }

  const auto graph_start = Clock::now();
  const WeightedGraph graph = build_graph(candidates, config.graph, stats);
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
    if (density >= 0.50) {
      component_strategy = WassersteinComponentStrategy::none;
    }
  }

  Weight best_saving = 0;
  const auto solver_start = Clock::now();
  if (component_strategy == WassersteinComponentStrategy::none) {
    best_saving = solve_graph(graph, config.matcher, config.warm_start, stats);
  } else {
    const auto component_start = Clock::now();
    const std::vector<Component> components = graph_components(graph, stats);
    const auto component_stop = Clock::now();
    if (stats != nullptr) {
      stats->component_time_ns += elapsed_ns(component_start, component_stop);
    }
    // Reuse the already measured decomposition through a local loop rather than
    // calling solve_components, which would discover the same components again.
    for (const Component& component : components) {
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
        best_saving += exhaustive_component(local);
        continue;
      }
      WassersteinMatcherStrategy local_matcher = config.matcher;
      if (component_strategy == WassersteinComponentStrategy::dense) {
        local_matcher = WassersteinMatcherStrategy::dense_sap;
      } else if (component_strategy == WassersteinComponentStrategy::sparse ||
                 component_strategy == WassersteinComponentStrategy::tiny_sparse) {
        local_matcher = WassersteinMatcherStrategy::sparse_sap;
      } else if (component_strategy == WassersteinComponentStrategy::adaptive) {
        local_matcher = (std::max)(local.rows, local.columns) <= 24 || density >= 0.20
                            ? WassersteinMatcherStrategy::dense_sap
                            : WassersteinMatcherStrategy::sparse_sap;
      }
      best_saving += solve_graph(local, local_matcher, config.warm_start, stats);
    }
  }
  const auto solver_stop = Clock::now();
  if (stats != nullptr) {
    stats->solver_time_ns += elapsed_ns(solver_start, solver_stop);
  }

  const Weight powered_distance =
      (std::max)(Weight{0}, diagonal_baseline - best_saving);
  return config.metric == WassersteinMetric::w1_linf
             ? static_cast<double>(powered_distance)
             : static_cast<double>(std::sqrt(powered_distance));
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
    case WassersteinCandidateStrategy::sweep_binary:
      return "sweep_binary";
    case WassersteinCandidateStrategy::sweep_two_pointer:
      return "sweep_two_pointer";
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
    case WassersteinMatcherStrategy::sparse_sap:
      return "sparse_sap";
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

}  // namespace bottleneck
