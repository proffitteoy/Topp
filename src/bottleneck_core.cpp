#include <bottleneck/core.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#if defined(BOTTLENECK_HAVE_AVX2_KERNEL) && defined(_MSC_VER)
#include <intrin.h>
#endif

#if defined(BOTTLENECK_HAVE_AVX2_KERNEL)
namespace bottleneck::detail {
void fill_cross_row_avx2(double first_birth, double first_death, const double* second_births,
                         const double* second_deaths, std::size_t size, double* output) noexcept;
}
#endif

namespace bottleneck {
namespace {

constexpr int unmatched = -1;

bool cpu_has_avx2() noexcept {
#if defined(BOTTLENECK_HAVE_AVX2_KERNEL) && defined(_MSC_VER) && defined(_M_X64)
  static const bool available = [] {
    int registers[4]{};
    __cpuid(registers, 1);
    constexpr int osxsave = 1 << 27;
    constexpr int avx = 1 << 28;
    if ((registers[2] & (osxsave | avx)) != (osxsave | avx) || (_xgetbv(0) & 0x6) != 0x6) {
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

double diagonal_distance(const Point& point) {
  const double midpoint = (point.birth + point.death) / 2.0;
  return (std::max)(std::fabs(point.birth - midpoint), std::fabs(point.death - midpoint));
}

double point_distance(const Point& first, const Point& second) {
  return (std::max)(std::fabs(first.birth - second.birth),
                    std::fabs(first.death - second.death));
}

struct DistanceTable {
  explicit DistanceTable(const PreparedDiagram& first_prepared,
                         const PreparedDiagram& second_prepared, DistanceStrategy strategy)
      : first(first_prepared.finite_points()), second(second_prepared.finite_points()),
        first_births(first_prepared.finite_births()), first_deaths(first_prepared.finite_deaths()),
        second_births(second_prepared.finite_births()), second_deaths(second_prepared.finite_deaths()),
        second_sorted_births(second_prepared.sorted_finite_births()),
        second_birth_order(second_prepared.finite_birth_order()),
        strategy(strategy), first_diagonal(first_prepared.finite_diagonal_distances()),
        second_diagonal(second_prepared.finite_diagonal_distances()),
        first_max_diagonal(first_prepared.max_finite_diagonal_distance()),
        second_max_diagonal(second_prepared.max_finite_diagonal_distance()) {
    if (strategy == DistanceStrategy::dense_aos || strategy == DistanceStrategy::dense_soa ||
        strategy == DistanceStrategy::dense_soa_avx2) {
      cross.resize(first.size() * second.size());
      for (std::size_t i = 0; i < first.size(); ++i) {
#if defined(BOTTLENECK_HAVE_AVX2_KERNEL)
        if (strategy == DistanceStrategy::dense_soa_avx2 && cpu_has_avx2()) {
          detail::fill_cross_row_avx2(first_births[i], first_deaths[i], second_births.data(),
                                      second_deaths.data(), second.size(),
                                      cross.data() + i * second.size());
          continue;
        }
#endif
        for (std::size_t j = 0; j < second.size(); ++j) {
          cross[i * second.size() + j] = compute_cross(i, j);
        }
      }
    }
  }

  [[nodiscard]] std::size_t side_size() const noexcept {
    return first.size() + second.size();
  }

  [[nodiscard]] double cross_distance(std::size_t first_index, std::size_t second_index) const {
    if (!cross.empty()) {
      return cross[first_index * second.size() + second_index];
    }
    return compute_cross(first_index, second_index);
  }

  [[nodiscard]] bool edge_allowed(std::size_t left, std::size_t right, double threshold) const {
    const std::size_t n = first.size();
    const std::size_t m = second.size();
    if (left < n) {
      if (right < m) {
        return cross_distance(left, right) <= threshold;
      }
      return right == m + left && first_diagonal[left] <= threshold;
    }

    const std::size_t projected_second = left - n;
    if (right < m) {
      return right == projected_second && second_diagonal[projected_second] <= threshold;
    }
    return true;
  }

  [[nodiscard]] double compute_cross(std::size_t first_index, std::size_t second_index) const {
    if (strategy == DistanceStrategy::dense_soa || strategy == DistanceStrategy::dense_soa_avx2 ||
        strategy == DistanceStrategy::recompute_soa) {
      return (std::max)(std::fabs(first_births[first_index] - second_births[second_index]),
                        std::fabs(first_deaths[first_index] - second_deaths[second_index]));
    }
    return point_distance(first[first_index], second[second_index]);
  }

  const Diagram& first;
  const Diagram& second;
  const std::vector<double>& first_births;
  const std::vector<double>& first_deaths;
  const std::vector<double>& second_births;
  const std::vector<double>& second_deaths;
  const std::vector<double>& second_sorted_births;
  const std::vector<std::size_t>& second_birth_order;
  DistanceStrategy strategy;
  const std::vector<double>& first_diagonal;
  const std::vector<double>& second_diagonal;
  double first_max_diagonal;
  double second_max_diagonal;
  std::vector<double> cross;
};

class ThresholdGraph {
 public:
  struct Cursor {
    std::size_t left = 0;
    std::size_t next = 0;
    std::size_t word_index = 0;
    std::size_t end = 0;
    std::uint64_t word = 0;
  };

  ThresholdGraph(const DistanceTable& table, double threshold, AdjacencyStrategy strategy,
                 SolverStats* stats)
      : table_(table), threshold_(threshold), strategy_(strategy), side_size_(table.side_size()) {
    if (strategy_ == AdjacencyStrategy::adaptive) {
      strategy_ = side_size_ >= 128 ? AdjacencyStrategy::x_sweep_csr
                                    : AdjacencyStrategy::bitset_auto;
    }
    if (strategy_ == AdjacencyStrategy::bitset_auto) {
      strategy_ = side_size_ <= 64    ? AdjacencyStrategy::bitset64
                  : side_size_ <= 128 ? AdjacencyStrategy::bitset128
                  : side_size_ <= 256 ? AdjacencyStrategy::bitset256
                                      : AdjacencyStrategy::dynamic_bitset;
    }
    if (strategy_ == AdjacencyStrategy::bitset64 && side_size_ > 64) {
      strategy_ = AdjacencyStrategy::dynamic_bitset;
    } else if (strategy_ == AdjacencyStrategy::bitset128 && side_size_ > 128) {
      strategy_ = AdjacencyStrategy::dynamic_bitset;
    } else if (strategy_ == AdjacencyStrategy::bitset256 && side_size_ > 256) {
      strategy_ = AdjacencyStrategy::dynamic_bitset;
    }
    if (strategy_ == AdjacencyStrategy::dense_byte) {
      dense_.resize(side_size_ * side_size_);
      for (std::size_t left = 0; left < side_size_; ++left) {
        for_each_neighbor(left, stats, [&](std::size_t right) {
          dense_[left * side_size_ + right] = 1;
        });
      }
    } else if (strategy_ == AdjacencyStrategy::sparse_csr ||
               strategy_ == AdjacencyStrategy::x_sweep_csr) {
      offsets_.reserve(side_size_ + 1);
      offsets_.push_back(0);
      for (std::size_t left = 0; left < side_size_; ++left) {
        const auto append = [&](std::size_t right) {
          sparse_neighbors_.push_back(static_cast<int>(right));
        };
        if (strategy_ == AdjacencyStrategy::x_sweep_csr) {
          for_each_x_sweep_neighbor(left, stats, append);
        } else {
          for_each_neighbor(left, stats, append);
        }
        offsets_.push_back(sparse_neighbors_.size());
      }
    } else if (strategy_ == AdjacencyStrategy::bitset64) {
      bits64_.assign(side_size_, 0);
      for (std::size_t left = 0; left < side_size_; ++left) {
        for_each_neighbor(left, stats, [&](std::size_t right) {
          bits64_[left] |= std::uint64_t{1} << right;
        });
      }
    } else if (strategy_ == AdjacencyStrategy::bitset128) {
      bits128_.resize(side_size_);
      for (std::size_t left = 0; left < side_size_; ++left) {
        for_each_neighbor(left, stats, [&](std::size_t right) {
          bits128_[left][right / 64U] |= std::uint64_t{1} << (right % 64U);
        });
      }
    } else if (strategy_ == AdjacencyStrategy::bitset256) {
      bits256_.resize(side_size_);
      for (std::size_t left = 0; left < side_size_; ++left) {
        for_each_neighbor(left, stats, [&](std::size_t right) {
          bits256_[left][right / 64U] |= std::uint64_t{1} << (right % 64U);
        });
      }
    } else if (strategy_ == AdjacencyStrategy::dynamic_bitset) {
      words_per_row_ = (side_size_ + 63U) / 64U;
      bits_.assign(side_size_ * words_per_row_, 0);
      for (std::size_t left = 0; left < side_size_; ++left) {
        for_each_neighbor(left, stats, [&](std::size_t right) {
          bits_[left * words_per_row_ + right / 64U] |= std::uint64_t{1} << (right % 64U);
        });
      }
    }
  }

  [[nodiscard]] std::size_t size() const noexcept { return side_size_; }

  [[nodiscard]] Cursor cursor(std::size_t left) const {
    Cursor result;
    result.left = left;
    if (strategy_ == AdjacencyStrategy::sparse_csr ||
        strategy_ == AdjacencyStrategy::x_sweep_csr) {
      result.next = offsets_[left];
      result.end = offsets_[left + 1];
    } else if (strategy_ == AdjacencyStrategy::bitset64) {
      result.word = bits64_[left];
    } else if (strategy_ == AdjacencyStrategy::bitset128) {
      result.word = bits128_[left][0];
      result.end = 2;
    } else if (strategy_ == AdjacencyStrategy::bitset256) {
      result.word = bits256_[left][0];
      result.end = 4;
    } else if (strategy_ == AdjacencyStrategy::dynamic_bitset && words_per_row_ != 0) {
      result.word = bits_[left * words_per_row_];
    }
    return result;
  }

  bool next(Cursor& cursor, int& right, SolverStats* stats) const {
    if (strategy_ == AdjacencyStrategy::sparse_csr ||
        strategy_ == AdjacencyStrategy::x_sweep_csr) {
      if (cursor.next == cursor.end) {
        return false;
      }
      right = sparse_neighbors_[cursor.next++];
      emit(stats);
      return true;
    }
    if (strategy_ == AdjacencyStrategy::bitset64) {
      if (cursor.word == 0) {
        return false;
      }
      const auto bit = static_cast<std::size_t>(std::countr_zero(cursor.word));
      cursor.word &= cursor.word - 1;
      right = static_cast<int>(bit);
      emit(stats);
      return true;
    }
    if (strategy_ == AdjacencyStrategy::bitset128 ||
        strategy_ == AdjacencyStrategy::bitset256) {
      while (cursor.word_index < cursor.end) {
        if (cursor.word != 0) {
          const auto bit = static_cast<std::size_t>(std::countr_zero(cursor.word));
          cursor.word &= cursor.word - 1;
          right = static_cast<int>(cursor.word_index * 64U + bit);
          emit(stats);
          return true;
        }
        ++cursor.word_index;
        if (cursor.word_index < cursor.end) {
          cursor.word = strategy_ == AdjacencyStrategy::bitset128
                            ? bits128_[cursor.left][cursor.word_index]
                            : bits256_[cursor.left][cursor.word_index];
        }
      }
      return false;
    }
    if (strategy_ == AdjacencyStrategy::dynamic_bitset) {
      while (cursor.word_index < words_per_row_) {
        if (cursor.word != 0) {
          const auto bit = static_cast<std::size_t>(std::countr_zero(cursor.word));
          cursor.word &= cursor.word - 1;
          const std::size_t index = cursor.word_index * 64U + bit;
          if (index < side_size_) {
            right = static_cast<int>(index);
            emit(stats);
            return true;
          }
        } else {
          ++cursor.word_index;
          if (cursor.word_index < words_per_row_) {
            cursor.word = bits_[cursor.left * words_per_row_ + cursor.word_index];
          }
        }
      }
      return false;
    }

    while (cursor.next < side_size_) {
      const std::size_t candidate = cursor.next++;
      bool allowed = false;
      if (strategy_ == AdjacencyStrategy::dense_byte) {
        allowed = dense_[cursor.left * side_size_ + candidate] != 0;
      } else {
        allowed = check(cursor.left, candidate, stats);
      }
      if (allowed) {
        right = static_cast<int>(candidate);
        emit(stats);
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] std::size_t degree(std::size_t left, SolverStats* stats) const {
    if (strategy_ == AdjacencyStrategy::sparse_csr ||
        strategy_ == AdjacencyStrategy::x_sweep_csr) {
      return offsets_[left + 1] - offsets_[left];
    }
    if (strategy_ == AdjacencyStrategy::bitset64) {
      return std::popcount(bits64_[left]);
    }
    if (strategy_ == AdjacencyStrategy::bitset128) {
      return std::popcount(bits128_[left][0]) + std::popcount(bits128_[left][1]);
    }
    if (strategy_ == AdjacencyStrategy::bitset256) {
      std::size_t result = 0;
      for (std::uint64_t word : bits256_[left]) {
        result += std::popcount(word);
      }
      return result;
    }
    if (strategy_ == AdjacencyStrategy::dynamic_bitset) {
      std::size_t result = 0;
      for (std::size_t word = 0; word < words_per_row_; ++word) {
        result += std::popcount(bits_[left * words_per_row_ + word]);
      }
      return result;
    }
    if (strategy_ == AdjacencyStrategy::dense_byte) {
      return static_cast<std::size_t>(std::count(
          dense_.begin() + static_cast<std::ptrdiff_t>(left * side_size_),
          dense_.begin() + static_cast<std::ptrdiff_t>((left + 1) * side_size_),
          std::uint8_t{1}));
    }

    std::size_t result = 0;
    for (std::size_t right = 0; right < side_size_; ++right) {
      result += static_cast<std::size_t>(check(left, right, stats));
    }
    return result;
  }

 private:
  template <class Function>
  void for_each_neighbor(std::size_t left, SolverStats* stats, Function&& function) const {
    const std::size_t n = table_.first.size();
    const std::size_t m = table_.second.size();
    if (left < n) {
      for (std::size_t right = 0; right < m; ++right) {
        if (check(left, right, stats)) {
          function(right);
        }
      }
      const std::size_t diagonal_right = m + left;
      if (check(left, diagonal_right, stats)) {
        function(diagonal_right);
      }
      return;
    }

    const std::size_t projected_second = left - n;
    if (check(left, projected_second, stats)) {
      function(projected_second);
    }
    for (std::size_t projected_first = 0; projected_first < n; ++projected_first) {
      function(m + projected_first);
    }
  }

  template <class Function>
  void for_each_x_sweep_neighbor(std::size_t left, SolverStats* stats,
                                 Function&& function) const {
    const std::size_t n = table_.first.size();
    const std::size_t m = table_.second.size();
    if (left < n) {
      const double birth = table_.first_births[left];
      const double lower = std::nextafter(birth - threshold_,
                                          -std::numeric_limits<double>::infinity());
      const double upper = std::nextafter(birth + threshold_,
                                          std::numeric_limits<double>::infinity());
      auto begin = std::lower_bound(table_.second_sorted_births.begin(),
                                    table_.second_sorted_births.end(), lower);
      const auto end = std::upper_bound(begin, table_.second_sorted_births.end(), upper);
      for (; begin != end; ++begin) {
        const std::size_t position =
            static_cast<std::size_t>(begin - table_.second_sorted_births.begin());
        const std::size_t right = table_.second_birth_order[position];
        if (stats != nullptr) {
          ++stats->adjacency_checks;
        }
        if (std::fabs(table_.first_births[left] - table_.second_births[right]) <= threshold_ &&
            std::fabs(table_.first_deaths[left] - table_.second_deaths[right]) <= threshold_) {
          function(right);
        }
      }
      const std::size_t diagonal_right = m + left;
      if (check(left, diagonal_right, stats)) {
        function(diagonal_right);
      }
      return;
    }

    const std::size_t projected_second = left - n;
    if (check(left, projected_second, stats)) {
      function(projected_second);
    }
    for (std::size_t projected_first = 0; projected_first < n; ++projected_first) {
      function(m + projected_first);
    }
  }

  bool check(std::size_t left, std::size_t right, SolverStats* stats) const {
    if (stats != nullptr) {
      ++stats->adjacency_checks;
    }
    return table_.edge_allowed(left, right, threshold_);
  }

  static void emit(SolverStats* stats) {
    if (stats != nullptr) {
      ++stats->emitted_edges;
    }
  }

  const DistanceTable& table_;
  double threshold_;
  AdjacencyStrategy strategy_;
  std::size_t side_size_;
  std::size_t words_per_row_ = 0;
  std::vector<std::uint8_t> dense_;
  std::vector<std::size_t> offsets_;
  std::vector<int> sparse_neighbors_;
  std::vector<std::uint64_t> bits64_;
  std::vector<std::array<std::uint64_t, 2>> bits128_;
  std::vector<std::array<std::uint64_t, 4>> bits256_;
  std::vector<std::uint64_t> bits_;
};

void fill_vertex_order(const ThresholdGraph& graph, VertexOrder order, SolverStats* stats,
                       std::vector<int>& vertices, std::vector<std::size_t>& degrees) {
  vertices.resize(graph.size());
  std::iota(vertices.begin(), vertices.end(), 0);
  if (order == VertexOrder::degree_ascending) {
    degrees.resize(graph.size());
    for (std::size_t vertex = 0; vertex < graph.size(); ++vertex) {
      degrees[vertex] = graph.degree(vertex, stats);
    }
    std::stable_sort(vertices.begin(), vertices.end(), [&](int first, int second) {
      return degrees[static_cast<std::size_t>(first)] < degrees[static_cast<std::size_t>(second)];
    });
  }
}

std::vector<int> vertex_order(const ThresholdGraph& graph, VertexOrder order, SolverStats* stats) {
  std::vector<int> vertices;
  std::vector<std::size_t> degrees;
  fill_vertex_order(graph, order, stats, vertices, degrees);
  return vertices;
}

bool kuhn_dfs(int left, const ThresholdGraph& graph, std::vector<int>& match_left,
              std::vector<int>& match_right, std::vector<std::uint32_t>& seen,
              std::uint32_t generation, SolverStats* stats) {
  auto cursor = graph.cursor(static_cast<std::size_t>(left));
  int right = 0;
  while (graph.next(cursor, right, stats)) {
    if (seen[static_cast<std::size_t>(right)] == generation) {
      continue;
    }
    seen[static_cast<std::size_t>(right)] = generation;
    const int previous = match_right[static_cast<std::size_t>(right)];
    if (previous == unmatched ||
        kuhn_dfs(previous, graph, match_left, match_right, seen, generation, stats)) {
      match_left[static_cast<std::size_t>(left)] = right;
      match_right[static_cast<std::size_t>(right)] = left;
      return true;
    }
  }
  return false;
}

void greedy_initialize(const ThresholdGraph& graph, const std::vector<int>& vertices,
                       std::vector<int>& match_left, std::vector<int>& match_right,
                       SolverStats* stats) {
  for (int left : vertices) {
    auto cursor = graph.cursor(static_cast<std::size_t>(left));
    int right = 0;
    while (graph.next(cursor, right, stats)) {
      if (match_right[static_cast<std::size_t>(right)] == unmatched) {
        match_left[static_cast<std::size_t>(left)] = right;
        match_right[static_cast<std::size_t>(right)] = left;
        break;
      }
    }
  }
}

struct KuhnWorkspace {
  std::vector<int> vertices;
  std::vector<std::size_t> degrees;
  std::vector<int> match_left;
  std::vector<int> match_right;
  std::vector<std::uint32_t> seen;
};

bool perfect_matching_kuhn(const ThresholdGraph& graph, VertexOrder order, bool greedy, bool reuse,
                           SolverStats* stats) {
  KuhnWorkspace local_workspace;
  static thread_local KuhnWorkspace reusable_workspace;
  KuhnWorkspace& workspace = reuse ? reusable_workspace : local_workspace;
  fill_vertex_order(graph, order, stats, workspace.vertices, workspace.degrees);
  workspace.match_left.assign(graph.size(), unmatched);
  workspace.match_right.assign(graph.size(), unmatched);
  if (greedy) {
    greedy_initialize(graph, workspace.vertices, workspace.match_left, workspace.match_right, stats);
  }
  workspace.seen.assign(graph.size(), 0);
  std::uint32_t generation = 0;
  for (int left : workspace.vertices) {
    if (workspace.match_left[static_cast<std::size_t>(left)] != unmatched) {
      continue;
    }
    if (stats != nullptr) {
      ++stats->augment_searches;
    }
    ++generation;
    if (generation == 0) {
      std::fill(workspace.seen.begin(), workspace.seen.end(), 0);
      generation = 1;
    }
    if (!kuhn_dfs(left, graph, workspace.match_left, workspace.match_right, workspace.seen,
                  generation, stats)) {
      return false;
    }
  }
  return true;
}

bool fixed_kuhn_dfs(int left, const ThresholdGraph& graph, std::array<int, 256>& match_left,
                    std::array<int, 256>& match_right,
                    std::array<std::uint32_t, 256>& seen, std::uint32_t generation,
                    SolverStats* stats) {
  auto cursor = graph.cursor(static_cast<std::size_t>(left));
  int right = 0;
  while (graph.next(cursor, right, stats)) {
    const auto right_index = static_cast<std::size_t>(right);
    if (seen[right_index] == generation) {
      continue;
    }
    seen[right_index] = generation;
    const int previous = match_right[right_index];
    if (previous == unmatched ||
        fixed_kuhn_dfs(previous, graph, match_left, match_right, seen, generation, stats)) {
      match_left[static_cast<std::size_t>(left)] = right;
      match_right[right_index] = left;
      return true;
    }
  }
  return false;
}

bool perfect_matching_fixed_greedy_kuhn(const ThresholdGraph& graph, VertexOrder order,
                                        SolverStats* stats) {
  if (graph.size() > 256) {
    return perfect_matching_kuhn(graph, order, true, true, stats);
  }
  std::array<int, 256> vertices{};
  std::array<std::size_t, 256> degrees{};
  std::array<int, 256> match_left{};
  std::array<int, 256> match_right{};
  std::array<std::uint32_t, 256> seen{};
  std::fill_n(match_left.begin(), graph.size(), unmatched);
  std::fill_n(match_right.begin(), graph.size(), unmatched);
  for (std::size_t index = 0; index < graph.size(); ++index) {
    vertices[index] = static_cast<int>(index);
    if (order == VertexOrder::degree_ascending) {
      degrees[index] = graph.degree(index, stats);
    }
  }
  if (order == VertexOrder::degree_ascending) {
    std::stable_sort(vertices.begin(), vertices.begin() + static_cast<std::ptrdiff_t>(graph.size()),
                     [&](int first, int second) {
      return degrees[static_cast<std::size_t>(first)] <
             degrees[static_cast<std::size_t>(second)];
    });
  }
  for (std::size_t index = 0; index < graph.size(); ++index) {
    const int left = vertices[index];
    auto cursor = graph.cursor(static_cast<std::size_t>(left));
    int right = 0;
    while (graph.next(cursor, right, stats)) {
      if (match_right[static_cast<std::size_t>(right)] == unmatched) {
        match_left[static_cast<std::size_t>(left)] = right;
        match_right[static_cast<std::size_t>(right)] = left;
        break;
      }
    }
  }
  std::uint32_t generation = 0;
  for (std::size_t index = 0; index < graph.size(); ++index) {
    const int left = vertices[index];
    if (match_left[static_cast<std::size_t>(left)] != unmatched) {
      continue;
    }
    if (stats != nullptr) {
      ++stats->augment_searches;
    }
    ++generation;
    if (!fixed_kuhn_dfs(left, graph, match_left, match_right, seen, generation, stats)) {
      return false;
    }
  }
  return true;
}

struct MaterializedBitGraph {
  explicit MaterializedBitGraph(const ThresholdGraph& graph, SolverStats* stats)
      : size(graph.size()), words((size + 63U) / 64U), forward(size * words, 0),
        reverse(size * words, 0), degree_left(size, 0), degree_right(size, 0) {
    for (std::size_t left = 0; left < size; ++left) {
      auto cursor = graph.cursor(left);
      int right_value = 0;
      while (graph.next(cursor, right_value, stats)) {
        const auto right = static_cast<std::size_t>(right_value);
        forward[left * words + right / 64U] |= std::uint64_t{1} << (right % 64U);
        reverse[right * words + left / 64U] |= std::uint64_t{1} << (left % 64U);
        ++degree_left[left];
        ++degree_right[right];
      }
    }
  }

  template <class Function>
  void for_each(const std::vector<std::uint64_t>& rows, std::size_t row,
                const std::vector<std::uint64_t>& active, Function&& function) const {
    for (std::size_t word_index = 0; word_index < words; ++word_index) {
      std::uint64_t word = rows[row * words + word_index] & active[word_index];
      while (word != 0) {
        const std::size_t bit = static_cast<std::size_t>(std::countr_zero(word));
        function(word_index * 64U + bit);
        word &= word - 1;
      }
    }
  }

  [[nodiscard]] int first(const std::vector<std::uint64_t>& rows, std::size_t row,
                          const std::vector<std::uint64_t>& active) const {
    for (std::size_t word_index = 0; word_index < words; ++word_index) {
      const std::uint64_t word = rows[row * words + word_index] & active[word_index];
      if (word != 0) {
        return static_cast<int>(word_index * 64U +
                                static_cast<std::size_t>(std::countr_zero(word)));
      }
    }
    return unmatched;
  }

  std::size_t size;
  std::size_t words;
  std::vector<std::uint64_t> forward;
  std::vector<std::uint64_t> reverse;
  std::vector<std::size_t> degree_left;
  std::vector<std::size_t> degree_right;
};

std::vector<std::uint64_t> full_mask(std::size_t size) {
  std::vector<std::uint64_t> result((size + 63U) / 64U, ~std::uint64_t{0});
  if (!result.empty() && size % 64U != 0) {
    result.back() = (std::uint64_t{1} << (size % 64U)) - 1;
  }
  return result;
}

bool mask_contains(const std::vector<std::uint64_t>& mask, std::size_t index) {
  return (mask[index / 64U] & (std::uint64_t{1} << (index % 64U))) != 0;
}

void mask_remove(std::vector<std::uint64_t>& mask, std::size_t index) {
  mask[index / 64U] &= ~(std::uint64_t{1} << (index % 64U));
}

bool reduced_bit_kuhn_dfs(int left, const MaterializedBitGraph& graph,
                          const std::vector<std::uint64_t>& active_left,
                          const std::vector<std::uint64_t>& active_right,
                          std::vector<int>& match_right, std::vector<std::uint32_t>& seen,
                          std::uint32_t generation, SolverStats* stats) {
  const std::size_t left_index = static_cast<std::size_t>(left);
  for (std::size_t word_index = 0; word_index < graph.words; ++word_index) {
    std::uint64_t word = graph.forward[left_index * graph.words + word_index] &
                         active_right[word_index];
    while (word != 0) {
      const std::size_t bit = static_cast<std::size_t>(std::countr_zero(word));
      const std::size_t right = word_index * 64U + bit;
      word &= word - 1;
      if (stats != nullptr) {
        ++stats->emitted_edges;
      }
      if (seen[right] == generation) {
        continue;
      }
      seen[right] = generation;
      const int previous = match_right[right];
      if (previous == unmatched ||
          (mask_contains(active_left, static_cast<std::size_t>(previous)) &&
           reduced_bit_kuhn_dfs(previous, graph, active_left, active_right, match_right, seen,
                                generation, stats))) {
        match_right[right] = left;
        return true;
      }
    }
  }
  return false;
}

bool perfect_matching_constraint_kuhn(const ThresholdGraph& threshold_graph, VertexOrder order,
                                      SolverStats* stats) {
  MaterializedBitGraph graph(threshold_graph, stats);
  std::vector<std::uint64_t> active_left = full_mask(graph.size);
  std::vector<std::uint64_t> active_right = full_mask(graph.size);
  std::vector<std::size_t> queue;
  queue.reserve(graph.size * 2);
  for (std::size_t index = 0; index < graph.size; ++index) {
    if (graph.degree_left[index] == 0 || graph.degree_right[index] == 0) {
      return false;
    }
    if (graph.degree_left[index] == 1) {
      queue.push_back(index);
    }
    if (graph.degree_right[index] == 1) {
      queue.push_back(graph.size + index);
    }
  }

  std::size_t active_count = graph.size;
  for (std::size_t head = 0; head < queue.size(); ++head) {
    if (stats != nullptr) {
      ++stats->constraint_rounds;
    }
    const bool is_right = queue[head] >= graph.size;
    const std::size_t vertex = is_right ? queue[head] - graph.size : queue[head];
    if (is_right) {
      if (!mask_contains(active_right, vertex) || graph.degree_right[vertex] != 1) {
        continue;
      }
    } else if (!mask_contains(active_left, vertex) || graph.degree_left[vertex] != 1) {
      continue;
    }

    const int left_value = is_right ? graph.first(graph.reverse, vertex, active_left)
                                    : static_cast<int>(vertex);
    const int right_value = is_right ? static_cast<int>(vertex)
                                     : graph.first(graph.forward, vertex, active_right);
    if (left_value == unmatched || right_value == unmatched) {
      return false;
    }
    const auto left = static_cast<std::size_t>(left_value);
    const auto right = static_cast<std::size_t>(right_value);
    mask_remove(active_left, left);
    mask_remove(active_right, right);
    graph.degree_left[left] = 0;
    graph.degree_right[right] = 0;
    --active_count;
    if (stats != nullptr) {
      ++stats->forced_matches;
    }

    graph.for_each(graph.forward, left, active_right, [&](std::size_t neighbor_right) {
      if (--graph.degree_right[neighbor_right] <= 1) {
        queue.push_back(graph.size + neighbor_right);
      }
    });
    graph.for_each(graph.reverse, right, active_left, [&](std::size_t neighbor_left) {
      if (--graph.degree_left[neighbor_left] <= 1) {
        queue.push_back(neighbor_left);
      }
    });
  }

  if (active_count == 0) {
    return true;
  }
  for (std::size_t index = 0; index < graph.size; ++index) {
    if ((mask_contains(active_left, index) && graph.degree_left[index] == 0) ||
        (mask_contains(active_right, index) && graph.degree_right[index] == 0)) {
      return false;
    }
  }

  std::vector<int> vertices;
  vertices.reserve(active_count);
  for (std::size_t left = 0; left < graph.size; ++left) {
    if (mask_contains(active_left, left)) {
      vertices.push_back(static_cast<int>(left));
    }
  }
  if (order == VertexOrder::degree_ascending) {
    std::stable_sort(vertices.begin(), vertices.end(), [&](int first, int second) {
      return graph.degree_left[static_cast<std::size_t>(first)] <
             graph.degree_left[static_cast<std::size_t>(second)];
    });
  }

  std::vector<int> match_right(graph.size, unmatched);
  std::vector<std::uint32_t> seen(graph.size, 0);
  std::uint32_t generation = 0;
  for (int left : vertices) {
    if (stats != nullptr) {
      ++stats->augment_searches;
    }
    ++generation;
    if (!reduced_bit_kuhn_dfs(left, graph, active_left, active_right, match_right, seen,
                              generation, stats)) {
      return false;
    }
  }
  return true;
}

bool perfect_matching_component_kuhn(const ThresholdGraph& threshold_graph, VertexOrder order,
                                     SolverStats* stats) {
  MaterializedBitGraph graph(threshold_graph, stats);
  const auto all_vertices = full_mask(graph.size);
  std::vector<int> left_component(graph.size, unmatched);
  std::vector<int> right_component(graph.size, unmatched);
  std::vector<std::vector<int>> component_left_vertices;
  std::vector<std::size_t> queue;
  queue.reserve(graph.size * 2);

  for (std::size_t index = 0; index < graph.size; ++index) {
    if (graph.degree_left[index] == 0 || graph.degree_right[index] == 0) {
      if (stats != nullptr) {
        ++stats->component_rejects;
      }
      return false;
    }
  }

  for (std::size_t seed = 0; seed < graph.size; ++seed) {
    if (left_component[seed] != unmatched) {
      continue;
    }
    const int component = static_cast<int>(component_left_vertices.size());
    component_left_vertices.emplace_back();
    queue.clear();
    queue.push_back(seed);
    left_component[seed] = component;
    std::size_t left_count = 0;
    std::size_t right_count = 0;
    for (std::size_t head = 0; head < queue.size(); ++head) {
      const std::size_t encoded = queue[head];
      if (encoded < graph.size) {
        const std::size_t left = encoded;
        ++left_count;
        component_left_vertices.back().push_back(static_cast<int>(left));
        graph.for_each(graph.forward, left, all_vertices, [&](std::size_t right) {
          if (right_component[right] == unmatched) {
            right_component[right] = component;
            queue.push_back(graph.size + right);
          }
        });
      } else {
        const std::size_t right = encoded - graph.size;
        ++right_count;
        graph.for_each(graph.reverse, right, all_vertices, [&](std::size_t left) {
          if (left_component[left] == unmatched) {
            left_component[left] = component;
            queue.push_back(left);
          }
        });
      }
    }
    if (stats != nullptr) {
      ++stats->component_count;
    }
    if (left_count != right_count) {
      if (stats != nullptr) {
        ++stats->component_rejects;
      }
      return false;
    }
  }

  std::vector<int> match_right(graph.size, unmatched);
  std::vector<std::uint32_t> seen(graph.size, 0);
  std::uint32_t generation = 0;
  for (auto& vertices : component_left_vertices) {
    if (order == VertexOrder::degree_ascending) {
      std::stable_sort(vertices.begin(), vertices.end(), [&](int first, int second) {
        return graph.degree_left[static_cast<std::size_t>(first)] <
               graph.degree_left[static_cast<std::size_t>(second)];
      });
    }
    for (int left : vertices) {
      if (stats != nullptr) {
        ++stats->augment_searches;
      }
      ++generation;
      if (!reduced_bit_kuhn_dfs(left, graph, all_vertices, all_vertices, match_right, seen,
                                generation, stats)) {
        return false;
      }
    }
  }
  return true;
}

bool hk_bfs(const ThresholdGraph& graph, const std::vector<int>& match_left,
            const std::vector<int>& match_right, std::vector<int>& distance,
            const std::vector<int>& order, SolverStats* stats) {
  std::vector<int> queue;
  queue.reserve(graph.size());
  for (int left : order) {
    if (match_left[static_cast<std::size_t>(left)] == unmatched) {
      distance[static_cast<std::size_t>(left)] = 0;
      queue.push_back(left);
    } else {
      distance[static_cast<std::size_t>(left)] = unmatched;
    }
  }

  bool found_free_right = false;
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const int left = queue[head];
    auto cursor = graph.cursor(static_cast<std::size_t>(left));
    int right = 0;
    while (graph.next(cursor, right, stats)) {
      const int next_left = match_right[static_cast<std::size_t>(right)];
      if (next_left == unmatched) {
        found_free_right = true;
      } else if (distance[static_cast<std::size_t>(next_left)] == unmatched) {
        distance[static_cast<std::size_t>(next_left)] =
            distance[static_cast<std::size_t>(left)] + 1;
        queue.push_back(next_left);
      }
    }
  }
  return found_free_right;
}

bool hk_dfs(int left, const ThresholdGraph& graph, std::vector<int>& match_left,
            std::vector<int>& match_right, std::vector<int>& distance, SolverStats* stats) {
  auto cursor = graph.cursor(static_cast<std::size_t>(left));
  int right = 0;
  while (graph.next(cursor, right, stats)) {
    const int next_left = match_right[static_cast<std::size_t>(right)];
    if (next_left == unmatched ||
        (distance[static_cast<std::size_t>(next_left)] == distance[static_cast<std::size_t>(left)] + 1 &&
         hk_dfs(next_left, graph, match_left, match_right, distance, stats))) {
      match_left[static_cast<std::size_t>(left)] = right;
      match_right[static_cast<std::size_t>(right)] = left;
      return true;
    }
  }
  distance[static_cast<std::size_t>(left)] = unmatched;
  return false;
}

bool perfect_matching_hopcroft_karp(const ThresholdGraph& graph, VertexOrder order, bool greedy,
                                    SolverStats* stats) {
  const auto vertices = vertex_order(graph, order, stats);
  std::vector<int> match_left(graph.size(), unmatched);
  std::vector<int> match_right(graph.size(), unmatched);
  std::vector<int> distance(graph.size(), unmatched);
  if (greedy) {
    greedy_initialize(graph, vertices, match_left, match_right, stats);
  }
  std::size_t matching_size = static_cast<std::size_t>(
      std::count_if(match_left.begin(), match_left.end(), [](int value) { return value != unmatched; }));

  while (hk_bfs(graph, match_left, match_right, distance, vertices, stats)) {
    if (stats != nullptr) {
      ++stats->bfs_phases;
    }
    bool augmented = false;
    for (int left : vertices) {
      if (match_left[static_cast<std::size_t>(left)] == unmatched) {
        if (stats != nullptr) {
          ++stats->augment_searches;
        }
        if (hk_dfs(left, graph, match_left, match_right, distance, stats)) {
          ++matching_size;
          augmented = true;
        }
      }
    }
    if (!augmented) {
      break;
    }
  }
  return matching_size == graph.size();
}

struct FlowEdge {
  int target;
  int reverse;
  int capacity;
};

class Dinic {
 public:
  explicit Dinic(std::size_t size) : edges_(size), level_(size), cursor_(size) {}

  void add_edge(int source, int target, int capacity) {
    const int source_reverse = static_cast<int>(edges_[static_cast<std::size_t>(target)].size());
    const int target_reverse = static_cast<int>(edges_[static_cast<std::size_t>(source)].size());
    edges_[static_cast<std::size_t>(source)].push_back({target, source_reverse, capacity});
    edges_[static_cast<std::size_t>(target)].push_back({source, target_reverse, 0});
  }

  int max_flow(int source, int sink, SolverStats* stats) {
    int result = 0;
    while (build_levels(source, sink, stats)) {
      std::fill(cursor_.begin(), cursor_.end(), 0);
      while (const int pushed = push(source, sink, std::numeric_limits<int>::max(), stats)) {
        result += pushed;
      }
    }
    return result;
  }

 private:
  bool build_levels(int source, int sink, SolverStats* stats) {
    std::fill(level_.begin(), level_.end(), unmatched);
    std::vector<int> queue;
    queue.reserve(edges_.size());
    level_[static_cast<std::size_t>(source)] = 0;
    queue.push_back(source);
    for (std::size_t head = 0; head < queue.size(); ++head) {
      const int vertex = queue[head];
      for (const FlowEdge& edge : edges_[static_cast<std::size_t>(vertex)]) {
        if (edge.capacity > 0 && level_[static_cast<std::size_t>(edge.target)] == unmatched) {
          level_[static_cast<std::size_t>(edge.target)] =
              level_[static_cast<std::size_t>(vertex)] + 1;
          queue.push_back(edge.target);
        }
      }
    }
    if (stats != nullptr) {
      ++stats->bfs_phases;
    }
    return level_[static_cast<std::size_t>(sink)] != unmatched;
  }

  int push(int vertex, int sink, int flow, SolverStats* stats) {
    if (vertex == sink) {
      return flow;
    }
    auto& adjacency = edges_[static_cast<std::size_t>(vertex)];
    for (std::size_t& index = cursor_[static_cast<std::size_t>(vertex)]; index < adjacency.size();
         ++index) {
      FlowEdge& edge = adjacency[index];
      if (edge.capacity == 0 ||
          level_[static_cast<std::size_t>(edge.target)] !=
              level_[static_cast<std::size_t>(vertex)] + 1) {
        continue;
      }
      if (stats != nullptr) {
        ++stats->emitted_edges;
      }
      const int pushed = push(edge.target, sink, (std::min)(flow, edge.capacity), stats);
      if (pushed != 0) {
        edge.capacity -= pushed;
        edges_[static_cast<std::size_t>(edge.target)][static_cast<std::size_t>(edge.reverse)]
            .capacity += pushed;
        return pushed;
      }
    }
    return 0;
  }

  std::vector<std::vector<FlowEdge>> edges_;
  std::vector<int> level_;
  std::vector<std::size_t> cursor_;
};

bool mandatory_flow_within(const DistanceTable& table, double threshold, SolverStats* stats) {
  const std::size_t n = table.first.size();
  const std::size_t m = table.second.size();
  const int source = static_cast<int>(n + m);
  const int sink = source + 1;
  const int super_source = sink + 1;
  const int super_sink = super_source + 1;
  Dinic flow(static_cast<std::size_t>(super_sink + 1));
  std::vector<int> demand(static_cast<std::size_t>(super_sink + 1), 0);

  const auto add_bounded_edge = [&](int from, int to, int lower, int upper) {
    flow.add_edge(from, to, upper - lower);
    demand[static_cast<std::size_t>(from)] -= lower;
    demand[static_cast<std::size_t>(to)] += lower;
  };

  for (std::size_t left = 0; left < n; ++left) {
    const int mandatory = table.first_diagonal[left] > threshold ? 1 : 0;
    add_bounded_edge(source, static_cast<int>(left), mandatory, 1);
  }
  for (std::size_t right = 0; right < m; ++right) {
    const int mandatory = table.second_diagonal[right] > threshold ? 1 : 0;
    add_bounded_edge(static_cast<int>(n + right), sink, mandatory, 1);
  }
  for (std::size_t left = 0; left < n; ++left) {
    for (std::size_t right = 0; right < m; ++right) {
      if (stats != nullptr) {
        ++stats->adjacency_checks;
      }
      if (table.cross_distance(left, right) <= threshold) {
        flow.add_edge(static_cast<int>(left), static_cast<int>(n + right), 1);
      }
    }
  }
  flow.add_edge(sink, source, static_cast<int>(n + m));

  int required = 0;
  for (int vertex = 0; vertex <= sink; ++vertex) {
    if (demand[static_cast<std::size_t>(vertex)] > 0) {
      flow.add_edge(super_source, vertex, demand[static_cast<std::size_t>(vertex)]);
      required += demand[static_cast<std::size_t>(vertex)];
    } else if (demand[static_cast<std::size_t>(vertex)] < 0) {
      flow.add_edge(vertex, super_sink, -demand[static_cast<std::size_t>(vertex)]);
    }
  }
  return flow.max_flow(super_source, super_sink, stats) == required;
}

double sampled_cross_density(const DistanceTable& table, double threshold, SolverStats* stats) {
  const std::size_t total = table.first.size() * table.second.size();
  if (total == 0) {
    return 0.0;
  }
  const std::size_t sample_count = (std::min)(std::size_t{64}, total);
  std::size_t allowed = 0;
  for (std::size_t sample = 0; sample < sample_count; ++sample) {
    const std::size_t flat = sample * total / sample_count;
    const std::size_t left = flat / table.second.size();
    const std::size_t right = flat % table.second.size();
    if (stats != nullptr) {
      ++stats->adjacency_checks;
    }
    allowed += static_cast<std::size_t>(table.cross_distance(left, right) <= threshold);
  }
  return static_cast<double>(allowed) / static_cast<double>(sample_count);
}

bool finite_within(const DistanceTable& table, double threshold, const SolverConfig& config,
                    SolverStats* stats) {
  if (table.side_size() == 0) {
    return true;
  }
  if (std::fabs(table.first_max_diagonal - table.second_max_diagonal) > threshold) {
    if (stats != nullptr) {
      ++stats->lower_bound_rejects;
    }
    return false;
  }
  MatcherStrategy matcher = config.matcher;
  if (matcher == MatcherStrategy::adaptive) {
    const bool large = table.side_size() >= 384;
    matcher = large && sampled_cross_density(table, threshold, stats) <= 0.03
                  ? MatcherStrategy::mandatory_flow
                  : MatcherStrategy::reusable_greedy_kuhn;
  }
  if (matcher == MatcherStrategy::mandatory_flow) {
    return mandatory_flow_within(table, threshold, stats);
  }
  ThresholdGraph graph(table, threshold, config.adjacency, stats);
  if (matcher == MatcherStrategy::hopcroft_karp ||
      matcher == MatcherStrategy::greedy_hopcroft_karp) {
    return perfect_matching_hopcroft_karp(
        graph, config.vertex_order, matcher == MatcherStrategy::greedy_hopcroft_karp, stats);
  }
  if (matcher == MatcherStrategy::constraint_kuhn) {
    return perfect_matching_constraint_kuhn(graph, config.vertex_order, stats);
  }
  if (matcher == MatcherStrategy::component_kuhn) {
    return perfect_matching_component_kuhn(graph, config.vertex_order, stats);
  }
  if (matcher == MatcherStrategy::fixed_greedy_kuhn) {
    return perfect_matching_fixed_greedy_kuhn(graph, config.vertex_order, stats);
  }
  return perfect_matching_kuhn(
      graph, config.vertex_order,
      matcher == MatcherStrategy::greedy_kuhn || matcher == MatcherStrategy::reusable_greedy_kuhn,
      matcher == MatcherStrategy::reusable_greedy_kuhn, stats);
}

std::vector<double> candidates(const DistanceTable& table, CandidateStrategy strategy,
                               bool ordered, SolverStats* stats) {
  std::vector<double> result;
  result.reserve(1 + table.first_diagonal.size() + table.second_diagonal.size() +
                 table.first.size() * table.second.size());
  result.push_back(0.0);
  result.insert(result.end(), table.first_diagonal.begin(), table.first_diagonal.end());
  result.insert(result.end(), table.second_diagonal.begin(), table.second_diagonal.end());
  std::uint64_t clipped = 0;
  if (strategy == CandidateStrategy::sort_unique_clipped ||
      strategy == CandidateStrategy::sort_unique_greedy_clipped) {
    double upper_bound = (std::max)(table.first_max_diagonal, table.second_max_diagonal);
    if (strategy == CandidateStrategy::sort_unique_greedy_clipped) {
      struct CrossCandidate {
        double distance;
        std::size_t first;
        std::size_t second;
      };
      std::vector<CrossCandidate> cross;
      cross.reserve(table.first.size() * table.second.size());
      for (std::size_t i = 0; i < table.first.size(); ++i) {
        for (std::size_t j = 0; j < table.second.size(); ++j) {
          cross.push_back({table.cross_distance(i, j), i, j});
        }
      }
      std::sort(cross.begin(), cross.end(), [](const CrossCandidate& first,
                                               const CrossCandidate& second) {
        return first.distance < second.distance;
      });
      std::vector<std::uint8_t> matched_first(table.first.size(), 0);
      std::vector<std::uint8_t> matched_second(table.second.size(), 0);
      double matched_max = 0.0;
      for (const CrossCandidate& candidate : cross) {
        if (candidate.distance > upper_bound) {
          break;
        }
        if (matched_first[candidate.first] == 0 && matched_second[candidate.second] == 0) {
          matched_first[candidate.first] = 1;
          matched_second[candidate.second] = 1;
          matched_max = (std::max)(matched_max, candidate.distance);
        }
      }
      double greedy_upper_bound = matched_max;
      for (std::size_t i = 0; i < table.first.size(); ++i) {
        if (matched_first[i] == 0) {
          greedy_upper_bound = (std::max)(greedy_upper_bound, table.first_diagonal[i]);
        }
      }
      for (std::size_t j = 0; j < table.second.size(); ++j) {
        if (matched_second[j] == 0) {
          greedy_upper_bound = (std::max)(greedy_upper_bound, table.second_diagonal[j]);
        }
      }
      upper_bound = (std::min)(upper_bound, greedy_upper_bound);
      const auto retained_end = std::remove_if(result.begin(), result.end(),
                                               [&](double value) {
        return value > upper_bound;
      });
      clipped += static_cast<std::uint64_t>(result.end() - retained_end);
      result.erase(retained_end, result.end());
      for (const CrossCandidate& candidate : cross) {
        if (candidate.distance <= upper_bound) {
          result.push_back(candidate.distance);
        } else {
          ++clipped;
        }
      }
    } else {
      for (std::size_t i = 0; i < table.first.size(); ++i) {
        for (std::size_t j = 0; j < table.second.size(); ++j) {
          const double value = table.cross_distance(i, j);
          if (value <= upper_bound) {
            result.push_back(value);
          } else {
            ++clipped;
          }
        }
      }
    }
  } else {
    for (std::size_t i = 0; i < table.first.size(); ++i) {
      for (std::size_t j = 0; j < table.second.size(); ++j) {
        result.push_back(table.cross_distance(i, j));
      }
    }
  }
  if (stats != nullptr) {
    stats->raw_candidates += result.size() + clipped;
    stats->clipped_candidates += clipped;
  }
  if (ordered) {
    std::sort(result.begin(), result.end());
    if (strategy != CandidateStrategy::sort_all) {
      result.erase(std::unique(result.begin(), result.end()), result.end());
    }
  }
  if (stats != nullptr) {
    stats->retained_candidates += result.size();
  }
  return result;
}

struct WeightedEdge {
  double weight;
  int left;
  int right;
};

bool incremental_dfs(int left, const std::vector<std::vector<int>>& adjacency,
                     std::vector<int>& match_left, std::vector<int>& match_right,
                     std::vector<std::uint32_t>& seen, std::uint32_t generation,
                     SolverStats* stats) {
  for (int right : adjacency[static_cast<std::size_t>(left)]) {
    if (stats != nullptr) {
      ++stats->emitted_edges;
    }
    if (seen[static_cast<std::size_t>(right)] == generation) {
      continue;
    }
    seen[static_cast<std::size_t>(right)] = generation;
    const int previous = match_right[static_cast<std::size_t>(right)];
    if (previous == unmatched ||
        incremental_dfs(previous, adjacency, match_left, match_right, seen, generation, stats)) {
      match_left[static_cast<std::size_t>(left)] = right;
      match_right[static_cast<std::size_t>(right)] = left;
      return true;
    }
  }
  return false;
}

double finite_distance_incremental(const DistanceTable& table, const SolverConfig& config,
                                   bool blocked, SolverStats* stats) {
  const std::size_t n = table.first.size();
  const std::size_t m = table.second.size();
  const std::size_t side_size = table.side_size();
  if (side_size == 0) {
    return 0.0;
  }

  const double upper_bound = (std::max)(table.first_max_diagonal, table.second_max_diagonal);

  std::vector<WeightedEdge> edges;
  edges.reserve(2 * n * m + n + m);
  for (std::size_t i = 0; i < n; ++i) {
    edges.push_back({table.first_diagonal[i], static_cast<int>(i), static_cast<int>(m + i)});
  }
  for (std::size_t j = 0; j < m; ++j) {
    edges.push_back({table.second_diagonal[j], static_cast<int>(n + j), static_cast<int>(j)});
  }
  for (std::size_t j = 0; j < m; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      edges.push_back({0.0, static_cast<int>(n + j), static_cast<int>(m + i)});
    }
  }

  std::uint64_t clipped = 0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
      const double weight = table.cross_distance(i, j);
      if (weight <= upper_bound) {
        edges.push_back({weight, static_cast<int>(i), static_cast<int>(j)});
      } else {
        ++clipped;
      }
    }
  }

  std::sort(edges.begin(), edges.end(), [](const WeightedEdge& first, const WeightedEdge& second) {
    if (first.weight != second.weight) {
      return first.weight < second.weight;
    }
    if (first.left != second.left) {
      return first.left < second.left;
    }
    return first.right < second.right;
  });

  std::vector<std::size_t> group_starts;
  group_starts.reserve(edges.size() + 1);
  for (std::size_t index = 0; index < edges.size(); ++index) {
    if (index == 0 || edges[index].weight != edges[index - 1].weight) {
      group_starts.push_back(index);
    }
  }
  group_starts.push_back(edges.size());

  if (stats != nullptr) {
    stats->raw_candidates += 1 + n + m + n * m;
    stats->clipped_candidates += clipped;
    stats->retained_candidates += group_starts.size() - 1;
  }

  std::vector<std::vector<int>> adjacency(side_size);
  std::vector<int> match_left(side_size, unmatched);
  std::vector<int> match_right(side_size, unmatched);
  std::vector<std::uint32_t> seen(side_size, 0);
  std::uint32_t generation = 0;
  std::size_t matching_size = 0;

  const auto add_groups = [&](std::size_t first_group, std::size_t end_group) {
    for (std::size_t index = group_starts[first_group]; index < group_starts[end_group]; ++index) {
      adjacency[static_cast<std::size_t>(edges[index].left)].push_back(edges[index].right);
    }
  };

  const auto augment = [&]() {
    if (stats != nullptr) {
      ++stats->threshold_decisions;
    }
    std::vector<int> unmatched_left;
    unmatched_left.reserve(side_size - matching_size);
    for (std::size_t left = 0; left < side_size; ++left) {
      if (match_left[left] == unmatched) {
        unmatched_left.push_back(static_cast<int>(left));
      }
    }
    if (config.vertex_order == VertexOrder::degree_ascending) {
      std::stable_sort(unmatched_left.begin(), unmatched_left.end(), [&](int first, int second) {
        return adjacency[static_cast<std::size_t>(first)].size() <
               adjacency[static_cast<std::size_t>(second)].size();
      });
    }

    for (int left : unmatched_left) {
      if (stats != nullptr) {
        ++stats->augment_searches;
      }
      ++generation;
      if (generation == 0) {
        std::fill(seen.begin(), seen.end(), 0);
        generation = 1;
      }
      if (incremental_dfs(left, adjacency, match_left, match_right, seen, generation, stats)) {
        ++matching_size;
      }
    }
    return matching_size == side_size;
  };

  const std::size_t group_count = group_starts.size() - 1;
  if (!blocked) {
    for (std::size_t group = 0; group < group_count; ++group) {
      add_groups(group, group + 1);
      if (augment()) {
        return edges[group_starts[group]].weight;
      }
    }
    return upper_bound;
  }

  const std::size_t block_size = (std::max)(
      std::size_t{1}, static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(group_count)))));
  for (std::size_t block_begin = 0; block_begin < group_count; block_begin += block_size) {
    const std::size_t block_end = (std::min)(group_count, block_begin + block_size);
    std::vector<std::size_t> row_sizes(side_size);
    for (std::size_t left = 0; left < side_size; ++left) {
      row_sizes[left] = adjacency[left].size();
    }
    const std::vector<int> saved_match_left = match_left;
    const std::vector<int> saved_match_right = match_right;
    const std::size_t saved_matching_size = matching_size;

    add_groups(block_begin, block_end);
    if (!augment()) {
      continue;
    }

    for (std::size_t left = 0; left < side_size; ++left) {
      adjacency[left].resize(row_sizes[left]);
    }
    match_left = saved_match_left;
    match_right = saved_match_right;
    matching_size = saved_matching_size;
    for (std::size_t group = block_begin; group < block_end; ++group) {
      add_groups(group, group + 1);
      if (augment()) {
        return edges[group_starts[group]].weight;
      }
    }
    return upper_bound;
  }
  return upper_bound;
}

double finite_distance(const PreparedDiagram& first, const PreparedDiagram& second,
                       const SolverConfig& config, SolverStats* stats) {
  const bool use_x_sweep = config.adjacency == AdjacencyStrategy::x_sweep_csr ||
                           (config.adjacency == AdjacencyStrategy::adaptive &&
                            first.finite_points().size() + second.finite_points().size() >= 128);
  if (use_x_sweep &&
      first.finite_points().size() > second.finite_points().size()) {
    return finite_distance(second, first, config, stats);
  }
  const DistanceTable table(first, second, config.distance);
  if (config.threshold == ThresholdStrategy::incremental ||
      config.threshold == ThresholdStrategy::incremental_blocked) {
    return finite_distance_incremental(table, config,
                                       config.threshold == ThresholdStrategy::incremental_blocked,
                                       stats);
  }
  const bool quickselect = config.threshold == ThresholdStrategy::quickselect;
  auto radii = candidates(table, config.candidates, !quickselect, stats);

  const auto decision = [&](std::size_t index) {
    if (stats != nullptr) {
      ++stats->threshold_decisions;
    }
    return finite_within(table, radii[index], config, stats);
  };

  const double lower_bound =
      std::fabs(table.first_max_diagonal - table.second_max_diagonal);
  if (quickselect) {
    radii.erase(std::remove_if(radii.begin(), radii.end(),
                               [&](double value) { return value < lower_bound; }),
                radii.end());
    double best = (std::max)(table.first_max_diagonal, table.second_max_diagonal);
    auto begin = radii.begin();
    auto end = radii.end();
    while (begin != end) {
      auto middle = begin + (end - begin) / 2;
      std::nth_element(begin, middle, end);
      const double pivot = *middle;
      auto less_end = std::partition(begin, end, [&](double value) { return value < pivot; });
      auto equal_end =
          std::partition(less_end, end, [&](double value) { return value == pivot; });
      bool feasible = pivot >= best;
      if (!feasible) {
        if (stats != nullptr) {
          ++stats->threshold_decisions;
        }
        feasible = finite_within(table, pivot, config, stats);
      }
      if (feasible) {
        best = (std::min)(best, pivot);
        end = less_end;
      } else {
        begin = equal_end;
      }
    }
    return best;
  }
  std::size_t lower = static_cast<std::size_t>(
      std::lower_bound(radii.begin(), radii.end(), lower_bound) - radii.begin());
  std::size_t upper = radii.size() - 1;
  if (config.threshold == ThresholdStrategy::exponential && lower < upper) {
    if (decision(lower)) {
      return radii[lower];
    }
    const std::size_t origin = lower;
    std::size_t last_false = lower;
    std::size_t probe = lower + 1;
    while (probe < upper && !decision(probe)) {
      last_false = probe;
      probe = (std::min)(upper, origin + (probe - origin) * 2);
    }
    lower = last_false + 1;
    upper = probe;
  }

  while (lower < upper) {
    std::size_t middle = lower + (upper - lower) / 2;
    if (config.threshold == ThresholdStrategy::gudhi_alpha && table.side_size() > 1) {
      const double alpha = std::pow(static_cast<double>(table.side_size()), 1.0 / 5.0);
      middle = lower + static_cast<std::size_t>(
                           static_cast<double>(upper - lower - 1) / alpha);
    }
    if (decision(middle)) {
      upper = middle;
    } else {
      lower = middle + 1;
    }
  }
  return radii[lower];
}

double essential_distance(const PreparedDiagram& first, const PreparedDiagram& second) {
  if (first.fully_infinite_count() != second.fully_infinite_count() ||
      first.positive_infinite_births().size() != second.positive_infinite_births().size() ||
      first.negative_infinite_deaths().size() != second.negative_infinite_deaths().size()) {
    return std::numeric_limits<double>::infinity();
  }

  double result = 0.0;
  for (std::size_t i = 0; i < first.positive_infinite_births().size(); ++i) {
    result = (std::max)(result, std::fabs(first.positive_infinite_births()[i] -
                                         second.positive_infinite_births()[i]));
  }
  for (std::size_t i = 0; i < first.negative_infinite_deaths().size(); ++i) {
    result = (std::max)(result, std::fabs(first.negative_infinite_deaths()[i] -
                                         second.negative_infinite_deaths()[i]));
  }
  return result;
}

}  // namespace

PreparedDiagram::PreparedDiagram(const Diagram& diagram) {
  const double infinity = std::numeric_limits<double>::infinity();
  const double negative_infinity = -infinity;
  finite_points_.reserve(diagram.size());
  finite_births_.reserve(diagram.size());
  finite_deaths_.reserve(diagram.size());
  finite_diagonal_distances_.reserve(diagram.size());
  for (const Point& point : diagram) {
    // Match GUDHI e=0 preprocessing semantics for the exact oracle.
    if (point.birth == infinity || point.death == negative_infinity) {
      continue;
    }
    if (point.birth == negative_infinity && point.death == infinity) {
      ++fully_infinite_count_;
    } else if (point.birth == negative_infinity) {
      negative_infinite_deaths_.push_back(point.death);
    } else if (point.death == infinity) {
      positive_infinite_births_.push_back(point.birth);
    } else if (point.death - point.birth > 0.0) {
      finite_points_.push_back(point);
      finite_births_.push_back(point.birth);
      finite_deaths_.push_back(point.death);
      const double distance = diagonal_distance(point);
      finite_diagonal_distances_.push_back(distance);
      max_finite_diagonal_distance_ = (std::max)(max_finite_diagonal_distance_, distance);
    }
  }
  finite_birth_order_.resize(finite_births_.size());
  std::iota(finite_birth_order_.begin(), finite_birth_order_.end(), std::size_t{0});
  std::stable_sort(finite_birth_order_.begin(), finite_birth_order_.end(),
                   [&](std::size_t first, std::size_t second) {
    return finite_births_[first] < finite_births_[second];
  });
  sorted_finite_births_.reserve(finite_birth_order_.size());
  for (std::size_t index : finite_birth_order_) {
    sorted_finite_births_.push_back(finite_births_[index]);
  }
  std::sort(positive_infinite_births_.begin(), positive_infinite_births_.end());
  std::sort(negative_infinite_deaths_.begin(), negative_infinite_deaths_.end());
}

const Diagram& PreparedDiagram::finite_points() const noexcept { return finite_points_; }

const std::vector<double>& PreparedDiagram::finite_births() const noexcept { return finite_births_; }

const std::vector<double>& PreparedDiagram::finite_deaths() const noexcept { return finite_deaths_; }

const std::vector<double>& PreparedDiagram::finite_diagonal_distances() const noexcept {
  return finite_diagonal_distances_;
}

const std::vector<double>& PreparedDiagram::sorted_finite_births() const noexcept {
  return sorted_finite_births_;
}

const std::vector<std::size_t>& PreparedDiagram::finite_birth_order() const noexcept {
  return finite_birth_order_;
}

double PreparedDiagram::max_finite_diagonal_distance() const noexcept {
  return max_finite_diagonal_distance_;
}

const std::vector<double>& PreparedDiagram::positive_infinite_births() const noexcept {
  return positive_infinite_births_;
}

const std::vector<double>& PreparedDiagram::negative_infinite_deaths() const noexcept {
  return negative_infinite_deaths_;
}

std::size_t PreparedDiagram::fully_infinite_count() const noexcept { return fully_infinite_count_; }

double bottleneck_distance(const PreparedDiagram& first, const PreparedDiagram& second,
                           const SolverConfig& config, SolverStats* stats) {
  const double alive = essential_distance(first, second);
  if (std::isinf(alive)) {
    return alive;
  }
  const double finite = finite_distance(first, second, config, stats);
  return (std::max)(alive, finite);
}

double bottleneck_distance(const Diagram& first, const Diagram& second, const SolverConfig& config,
                           SolverStats* stats) {
  return bottleneck_distance(PreparedDiagram(first), PreparedDiagram(second), config, stats);
}

bool bottleneck_within(const PreparedDiagram& first, const PreparedDiagram& second, double threshold,
                       const SolverConfig& config, SolverStats* stats) {
  if (threshold < 0.0 || std::isnan(threshold)) {
    return false;
  }
  const double alive = essential_distance(first, second);
  if (alive > threshold) {
    return false;
  }
  if (stats != nullptr) {
    ++stats->threshold_decisions;
  }
  const bool use_x_sweep = config.adjacency == AdjacencyStrategy::x_sweep_csr ||
                           (config.adjacency == AdjacencyStrategy::adaptive &&
                            first.finite_points().size() + second.finite_points().size() >= 128);
  if (use_x_sweep &&
      first.finite_points().size() > second.finite_points().size()) {
    const DistanceTable table(second, first, config.distance);
    return finite_within(table, threshold, config, stats);
  }
  const DistanceTable table(first, second, config.distance);
  return finite_within(table, threshold, config, stats);
}

void bottleneck_distances(const PreparedDiagram& query,
                          std::span<const PreparedDiagram> diagrams,
                          std::span<double> output, const SolverConfig& config,
                          SolverStats* stats) {
  if (output.size() != diagrams.size()) {
    throw std::invalid_argument("output size must match diagram count");
  }
  for (std::size_t index = 0; index < diagrams.size(); ++index) {
    output[index] = bottleneck_distance(query, diagrams[index], config, stats);
  }
}

std::vector<double> bottleneck_distances(const PreparedDiagram& query,
                                         std::span<const PreparedDiagram> diagrams,
                                         const SolverConfig& config, SolverStats* stats) {
  std::vector<double> result(diagrams.size());
  bottleneck_distances(query, diagrams, result, config, stats);
  return result;
}

const char* to_string(CandidateStrategy strategy) noexcept {
  switch (strategy) {
    case CandidateStrategy::sort_all:
      return "sort_all";
    case CandidateStrategy::sort_unique:
      return "sort_unique";
    case CandidateStrategy::sort_unique_clipped:
      return "sort_unique_clipped";
    case CandidateStrategy::sort_unique_greedy_clipped:
      return "sort_unique_greedy_clipped";
  }
  return "unknown";
}

const char* to_string(ThresholdStrategy strategy) noexcept {
  switch (strategy) {
    case ThresholdStrategy::binary:
      return "binary";
    case ThresholdStrategy::gudhi_alpha:
      return "gudhi_alpha";
    case ThresholdStrategy::exponential:
      return "exponential";
    case ThresholdStrategy::quickselect:
      return "quickselect";
    case ThresholdStrategy::incremental:
      return "incremental";
    case ThresholdStrategy::incremental_blocked:
      return "incremental_blocked";
  }
  return "unknown";
}

const char* to_string(DistanceStrategy strategy) noexcept {
  switch (strategy) {
    case DistanceStrategy::dense_aos:
      return "dense_aos";
    case DistanceStrategy::dense_soa:
      return "dense_soa";
    case DistanceStrategy::dense_soa_avx2:
      return "dense_soa_avx2";
    case DistanceStrategy::recompute_aos:
      return "recompute_aos";
    case DistanceStrategy::recompute_soa:
      return "recompute_soa";
  }
  return "unknown";
}

const char* to_string(AdjacencyStrategy strategy) noexcept {
  switch (strategy) {
    case AdjacencyStrategy::on_demand:
      return "on_demand";
    case AdjacencyStrategy::dense_byte:
      return "dense_byte";
    case AdjacencyStrategy::sparse_csr:
      return "sparse_csr";
    case AdjacencyStrategy::x_sweep_csr:
      return "x_sweep_csr";
    case AdjacencyStrategy::bitset64:
      return "bitset64";
    case AdjacencyStrategy::bitset128:
      return "bitset128";
    case AdjacencyStrategy::bitset256:
      return "bitset256";
    case AdjacencyStrategy::bitset_auto:
      return "bitset_auto";
    case AdjacencyStrategy::adaptive:
      return "adaptive";
    case AdjacencyStrategy::dynamic_bitset:
      return "dynamic_bitset";
  }
  return "unknown";
}

const char* to_string(MatcherStrategy strategy) noexcept {
  switch (strategy) {
    case MatcherStrategy::kuhn:
      return "kuhn";
    case MatcherStrategy::greedy_kuhn:
      return "greedy_kuhn";
    case MatcherStrategy::reusable_greedy_kuhn:
      return "reusable_greedy_kuhn";
    case MatcherStrategy::fixed_greedy_kuhn:
      return "fixed_greedy_kuhn";
    case MatcherStrategy::constraint_kuhn:
      return "constraint_kuhn";
    case MatcherStrategy::component_kuhn:
      return "component_kuhn";
    case MatcherStrategy::mandatory_flow:
      return "mandatory_flow";
    case MatcherStrategy::adaptive:
      return "adaptive";
    case MatcherStrategy::hopcroft_karp:
      return "hopcroft_karp";
    case MatcherStrategy::greedy_hopcroft_karp:
      return "greedy_hopcroft_karp";
  }
  return "unknown";
}

const char* to_string(VertexOrder order) noexcept {
  return order == VertexOrder::degree_ascending ? "degree_ascending" : "natural";
}

}  // namespace bottleneck
