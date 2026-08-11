#include "geometric_backend.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace bottleneck::detail {
namespace {

constexpr int unmatched = -1;

double point_distance(const Point& first, const Point& second) noexcept {
  return (std::max)(std::fabs(first.birth - second.birth),
                    std::fabs(first.death - second.death));
}

class KdRangeIndex {
 public:
  struct ActiveState {
    std::vector<std::uint8_t> active_points;
    std::vector<std::size_t> remaining;
  };

  explicit KdRangeIndex(const Diagram& points) : points_(points) {
    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    nodes_.reserve(points.size());
    root_ = build(indices, 0, indices.size(), 0);
  }

  template <class Function>
  bool for_each_in_square(const Point& query, double radius, SolverStats* stats,
                          Function&& function) const {
    if (stats != nullptr) {
      ++stats->geometric_queries;
    }
    const double min_birth = std::nextafter(query.birth - radius,
                                            -std::numeric_limits<double>::infinity());
    const double max_birth = std::nextafter(query.birth + radius,
                                            std::numeric_limits<double>::infinity());
    const double min_death = std::nextafter(query.death - radius,
                                            -std::numeric_limits<double>::infinity());
    const double max_death = std::nextafter(query.death + radius,
                                            std::numeric_limits<double>::infinity());
    return visit(root_, min_birth, max_birth, min_death, max_death, stats, function);
  }

  ActiveState active_state(const std::vector<int>& match_right) const {
    ActiveState state;
    state.active_points.resize(points_.size(), 0);
    for (std::size_t point = 0; point < points_.size(); ++point) {
      state.active_points[point] =
          static_cast<std::uint8_t>(match_right[point] == unmatched);
    }
    state.remaining.resize(nodes_.size(), 0);
    fill_remaining(root_, state);
    return state;
  }

  int take_any_in_square(const Point& query, double radius, ActiveState& state,
                         SolverStats* stats) const {
    if (stats != nullptr) {
      ++stats->geometric_queries;
    }
    return take(root_, query, radius,
                std::nextafter(query.birth - radius,
                               -std::numeric_limits<double>::infinity()),
                std::nextafter(query.birth + radius,
                               std::numeric_limits<double>::infinity()),
                std::nextafter(query.death - radius,
                               -std::numeric_limits<double>::infinity()),
                std::nextafter(query.death + radius,
                               std::numeric_limits<double>::infinity()),
                state, stats);
  }

 private:
  struct Node {
    std::size_t point = 0;
    int left = unmatched;
    int right = unmatched;
    double min_birth = 0.0;
    double max_birth = 0.0;
    double min_death = 0.0;
    double max_death = 0.0;
  };

  int build(std::vector<std::size_t>& indices, std::size_t begin, std::size_t end,
            unsigned depth) {
    if (begin == end) {
      return unmatched;
    }
    const std::size_t middle = begin + (end - begin) / 2;
    const bool split_birth = (depth % 2U) == 0;
    std::nth_element(indices.begin() + static_cast<std::ptrdiff_t>(begin),
                     indices.begin() + static_cast<std::ptrdiff_t>(middle),
                     indices.begin() + static_cast<std::ptrdiff_t>(end),
                     [&](std::size_t first, std::size_t second) {
                       return split_birth ? points_[first].birth < points_[second].birth
                                          : points_[first].death < points_[second].death;
                     });
    const int node_index = static_cast<int>(nodes_.size());
    nodes_.push_back({});
    const std::size_t point_index = indices[middle];
    const int left = build(indices, begin, middle, depth + 1);
    const int right = build(indices, middle + 1, end, depth + 1);
    Node& node = nodes_[static_cast<std::size_t>(node_index)];
    node.point = point_index;
    node.left = left;
    node.right = right;
    node.min_birth = node.max_birth = points_[point_index].birth;
    node.min_death = node.max_death = points_[point_index].death;
    include_bounds(node, left);
    include_bounds(node, right);
    return node_index;
  }

  void include_bounds(Node& node, int child_index) {
    if (child_index == unmatched) {
      return;
    }
    const Node& child = nodes_[static_cast<std::size_t>(child_index)];
    node.min_birth = (std::min)(node.min_birth, child.min_birth);
    node.max_birth = (std::max)(node.max_birth, child.max_birth);
    node.min_death = (std::min)(node.min_death, child.min_death);
    node.max_death = (std::max)(node.max_death, child.max_death);
  }

  template <class Function>
  bool visit(int node_index, double min_birth, double max_birth, double min_death,
             double max_death, SolverStats* stats, Function& function) const {
    if (node_index == unmatched) {
      return true;
    }
    const Node& node = nodes_[static_cast<std::size_t>(node_index)];
    if (node.max_birth < min_birth || node.min_birth > max_birth ||
        node.max_death < min_death || node.min_death > max_death) {
      return true;
    }
    if (stats != nullptr) {
      ++stats->kd_nodes_visited;
    }
    const Point& point = points_[node.point];
    if (point.birth >= min_birth && point.birth <= max_birth &&
        point.death >= min_death && point.death <= max_death) {
      if (stats != nullptr) {
        ++stats->adjacency_checks;
        ++stats->emitted_edges;
      }
      if (!function(node.point)) {
        return false;
      }
    }
    return visit(node.left, min_birth, max_birth, min_death, max_death, stats, function) &&
           visit(node.right, min_birth, max_birth, min_death, max_death, stats, function);
  }

  std::size_t fill_remaining(int node_index, ActiveState& state) const {
    if (node_index == unmatched) {
      return 0;
    }
    const Node& node = nodes_[static_cast<std::size_t>(node_index)];
    const std::size_t count =
        static_cast<std::size_t>(state.active_points[node.point]) +
        fill_remaining(node.left, state) + fill_remaining(node.right, state);
    state.remaining[static_cast<std::size_t>(node_index)] = count;
    return count;
  }

  int take(int node_index, const Point& query, double radius, double min_birth,
           double max_birth, double min_death, double max_death, ActiveState& state,
           SolverStats* stats) const {
    if (node_index == unmatched || state.remaining[static_cast<std::size_t>(node_index)] == 0) {
      return unmatched;
    }
    const Node& node = nodes_[static_cast<std::size_t>(node_index)];
    if (node.max_birth < min_birth || node.min_birth > max_birth ||
        node.max_death < min_death || node.min_death > max_death) {
      return unmatched;
    }
    if (stats != nullptr) {
      ++stats->kd_nodes_visited;
    }
    int result = unmatched;
    if (state.active_points[node.point] != 0 &&
        point_distance(query, points_[node.point]) <= radius) {
      state.active_points[node.point] = 0;
      result = static_cast<int>(node.point);
      if (stats != nullptr) {
        ++stats->adjacency_checks;
        ++stats->emitted_edges;
      }
    } else {
      result = take(node.left, query, radius, min_birth, max_birth, min_death, max_death,
                    state, stats);
      if (result == unmatched) {
        result = take(node.right, query, radius, min_birth, max_birth, min_death, max_death,
                      state, stats);
      }
    }
    if (result != unmatched) {
      --state.remaining[static_cast<std::size_t>(node_index)];
    }
    return result;
  }

  const Diagram& points_;
  std::vector<Node> nodes_;
  int root_ = unmatched;
};

class GeometricMatchOracle {
 public:
  GeometricMatchOracle(const PreparedDiagram& first, const PreparedDiagram& second)
      : first_(first), second_(second), second_index_(second.finite_points()),
        side_size_(first.finite_points().size() + second.finite_points().size()),
        match_left_(side_size_, unmatched), match_right_(side_size_, unmatched),
        distance_(side_size_, unmatched) {}

  bool within(double radius, SolverStats* stats) {
    if (radius < 0.0 || std::isnan(radius)) {
      return false;
    }
    if (std::fabs(first_.max_finite_diagonal_distance() -
                  second_.max_finite_diagonal_distance()) > radius) {
      if (stats != nullptr) {
        ++stats->lower_bound_rejects;
      }
      return false;
    }
    if (side_size_ == 0) {
      return true;
    }
    if (radius < current_radius_) {
      trim(radius);
    }
    current_radius_ = radius;
    greedy_fill(radius, stats);
    std::size_t matching_size = static_cast<std::size_t>(std::count_if(
        match_left_.begin(), match_left_.end(), [](int right) { return right != unmatched; }));
    if (matching_size == side_size_) {
      return true;
    }
    while (bfs(radius, stats)) {
      if (stats != nullptr) {
        ++stats->bfs_phases;
      }
      bool augmented = false;
      for (std::size_t left = 0; left < side_size_; ++left) {
        if (match_left_[left] == unmatched) {
          if (stats != nullptr) {
            ++stats->augment_searches;
          }
          if (dfs(static_cast<int>(left), radius, stats)) {
            ++matching_size;
            augmented = true;
          }
        }
      }
      if (!augmented || matching_size == side_size_) {
        break;
      }
    }
    return matching_size == side_size_;
  }

 private:
  bool edge_allowed(std::size_t left, std::size_t right, double radius) const {
    const std::size_t n = first_.finite_points().size();
    const std::size_t m = second_.finite_points().size();
    if (left < n) {
      if (right < m) {
        return point_distance(first_.finite_points()[left], second_.finite_points()[right]) <= radius;
      }
      return right == m + left && first_.finite_diagonal_distances()[left] <= radius;
    }
    const std::size_t projected_second = left - n;
    if (right < m) {
      return right == projected_second &&
             second_.finite_diagonal_distances()[projected_second] <= radius;
    }
    return true;
  }

  void trim(double radius) {
    for (std::size_t left = 0; left < side_size_; ++left) {
      const int right = match_left_[left];
      if (right != unmatched && !edge_allowed(left, static_cast<std::size_t>(right), radius)) {
        match_left_[left] = unmatched;
        match_right_[static_cast<std::size_t>(right)] = unmatched;
      }
    }
  }

  template <class Function>
  bool for_each_neighbor(std::size_t left, double radius, SolverStats* stats,
                         Function&& function) const {
    const std::size_t n = first_.finite_points().size();
    const std::size_t m = second_.finite_points().size();
    if (left < n) {
      auto cross = [&](std::size_t right) {
        if (point_distance(first_.finite_points()[left], second_.finite_points()[right]) > radius) {
          return true;
        }
        return function(static_cast<int>(right));
      };
      if (!second_index_.for_each_in_square(first_.finite_points()[left], radius, stats, cross)) {
        return false;
      }
      if (first_.finite_diagonal_distances()[left] <= radius) {
        if (stats != nullptr) {
          ++stats->adjacency_checks;
          ++stats->emitted_edges;
        }
        return function(static_cast<int>(m + left));
      }
      return true;
    }
    const std::size_t projected_second = left - n;
    if (second_.finite_diagonal_distances()[projected_second] <= radius) {
      if (stats != nullptr) {
        ++stats->adjacency_checks;
        ++stats->emitted_edges;
      }
      if (!function(static_cast<int>(projected_second))) {
        return false;
      }
    }
    for (std::size_t projected_first = 0; projected_first < n; ++projected_first) {
      if (stats != nullptr) {
        ++stats->emitted_edges;
      }
      if (!function(static_cast<int>(m + projected_first))) {
        return false;
      }
    }
    return true;
  }

  void greedy_fill(double radius, SolverStats* stats) {
    const std::size_t n = first_.finite_points().size();
    const std::size_t m = second_.finite_points().size();
    auto active_cross = second_index_.active_state(match_right_);
    for (std::size_t left = 0; left < n; ++left) {
      if (match_left_[left] != unmatched) {
        continue;
      }
      const int cross_right = second_index_.take_any_in_square(
          first_.finite_points()[left], radius, active_cross, stats);
      if (cross_right != unmatched) {
        match_left_[left] = cross_right;
        match_right_[static_cast<std::size_t>(cross_right)] = static_cast<int>(left);
      } else if (first_.finite_diagonal_distances()[left] <= radius &&
                 match_right_[m + left] == unmatched) {
        const int right = static_cast<int>(m + left);
        match_left_[left] = right;
        match_right_[static_cast<std::size_t>(right)] = static_cast<int>(left);
      }
    }

    std::vector<int> free_diagonal_right;
    free_diagonal_right.reserve(n);
    for (std::size_t right = m; right < side_size_; ++right) {
      if (match_right_[right] == unmatched) {
        free_diagonal_right.push_back(static_cast<int>(right));
      }
    }
    std::size_t next_diagonal = 0;
    for (std::size_t projected_second = 0; projected_second < m; ++projected_second) {
      const std::size_t left = n + projected_second;
      if (match_left_[left] != unmatched) {
        continue;
      }
      if (second_.finite_diagonal_distances()[projected_second] <= radius &&
          match_right_[projected_second] == unmatched) {
        match_left_[left] = static_cast<int>(projected_second);
        match_right_[projected_second] = static_cast<int>(left);
        continue;
      }
      while (next_diagonal < free_diagonal_right.size() &&
             match_right_[static_cast<std::size_t>(free_diagonal_right[next_diagonal])] != unmatched) {
        ++next_diagonal;
      }
      if (next_diagonal < free_diagonal_right.size()) {
        const int right = free_diagonal_right[next_diagonal++];
        match_left_[left] = right;
        match_right_[static_cast<std::size_t>(right)] = static_cast<int>(left);
      }
    }
  }

  bool bfs(double radius, SolverStats* stats) {
    std::vector<int> queue;
    queue.reserve(side_size_);
    for (std::size_t left = 0; left < side_size_; ++left) {
      if (match_left_[left] == unmatched) {
        distance_[left] = 0;
        queue.push_back(static_cast<int>(left));
      } else {
        distance_[left] = unmatched;
      }
    }
    bool found_free_right = false;
    for (std::size_t head = 0; head < queue.size(); ++head) {
      const int left = queue[head];
      for_each_neighbor(static_cast<std::size_t>(left), radius, stats, [&](int right) {
        const int next_left = match_right_[static_cast<std::size_t>(right)];
        if (next_left == unmatched) {
          found_free_right = true;
        } else if (distance_[static_cast<std::size_t>(next_left)] == unmatched) {
          distance_[static_cast<std::size_t>(next_left)] =
              distance_[static_cast<std::size_t>(left)] + 1;
          queue.push_back(next_left);
        }
        return true;
      });
    }
    return found_free_right;
  }

  bool dfs(int left, double radius, SolverStats* stats) {
    bool matched = false;
    for_each_neighbor(static_cast<std::size_t>(left), radius, stats, [&](int right) {
      const int next_left = match_right_[static_cast<std::size_t>(right)];
      if (next_left == unmatched ||
          (distance_[static_cast<std::size_t>(next_left)] ==
               distance_[static_cast<std::size_t>(left)] + 1 &&
           dfs(next_left, radius, stats))) {
        match_left_[static_cast<std::size_t>(left)] = right;
        match_right_[static_cast<std::size_t>(right)] = left;
        matched = true;
        return false;
      }
      return true;
    });
    if (!matched) {
      distance_[static_cast<std::size_t>(left)] = unmatched;
    }
    return matched;
  }

  const PreparedDiagram& first_;
  const PreparedDiagram& second_;
  KdRangeIndex second_index_;
  std::size_t side_size_;
  std::vector<int> match_left_;
  std::vector<int> match_right_;
  std::vector<int> distance_;
  double current_radius_ = std::numeric_limits<double>::infinity();
};

std::vector<double> annulus_candidates(const PreparedDiagram& first,
                                       const PreparedDiagram& second, double lower,
                                       double upper, SolverStats* stats) {
  std::vector<double> result;
  result.reserve(first.finite_points().size() + second.finite_points().size());
  const auto retain = [&](double value) {
    if (value > lower && value <= upper) {
      result.push_back(value);
    }
  };
  for (double value : first.finite_diagonal_distances()) {
    retain(value);
  }
  for (double value : second.finite_diagonal_distances()) {
    retain(value);
  }

  const auto& sorted_births = second.sorted_finite_births();
  const auto& birth_order = second.finite_birth_order();
  for (const Point& point : first.finite_points()) {
    const double min_birth = std::nextafter(point.birth - upper,
                                            -std::numeric_limits<double>::infinity());
    const double max_birth = std::nextafter(point.birth + upper,
                                            std::numeric_limits<double>::infinity());
    auto begin = std::lower_bound(sorted_births.begin(), sorted_births.end(), min_birth);
    const auto end = std::upper_bound(begin, sorted_births.end(), max_birth);
    for (; begin != end; ++begin) {
      const std::size_t position = static_cast<std::size_t>(begin - sorted_births.begin());
      const Point& other = second.finite_points()[birth_order[position]];
      if (stats != nullptr) {
        ++stats->adjacency_checks;
      }
      const double value = point_distance(point, other);
      retain(value);
    }
  }
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  if (stats != nullptr) {
    const std::uint64_t raw = 1 + first.finite_points().size() + second.finite_points().size() +
                              first.finite_points().size() * second.finite_points().size();
    stats->raw_candidates += raw;
    stats->retained_candidates += result.size();
    stats->clipped_candidates += raw - (std::min)(raw, static_cast<std::uint64_t>(result.size()));
  }
  return result;
}

}  // namespace

bool geometric_bottleneck_within(const PreparedDiagram& first, const PreparedDiagram& second,
                                 double threshold, SolverStats* stats) {
  GeometricMatchOracle oracle(first, second);
  return oracle.within(threshold, stats);
}

double geometric_refined_bottleneck_distance(const PreparedDiagram& first_input,
                                             const PreparedDiagram& second_input,
                                             SolverStats* stats) {
  if (first_input.finite_points().size() > second_input.finite_points().size()) {
    return geometric_refined_bottleneck_distance(second_input, first_input, stats);
  }
  const PreparedDiagram& first = first_input;
  const PreparedDiagram& second = second_input;
  if (first.finite_points().empty() && second.finite_points().empty()) {
    return 0.0;
  }

  GeometricMatchOracle oracle(first, second);
  const auto decide = [&](double radius) {
    if (stats != nullptr) {
      ++stats->threshold_decisions;
      ++stats->refinement_rounds;
    }
    return oracle.within(radius, stats);
  };

  double lower = std::fabs(first.max_finite_diagonal_distance() -
                           second.max_finite_diagonal_distance());
  double upper = (std::max)(first.max_finite_diagonal_distance(),
                            second.max_finite_diagonal_distance());
  if (decide(lower)) {
    return lower;
  }
  if (lower == 0.0) {
    double probe = upper / 2.0;
    while (probe > 0.0 && decide(probe)) {
      upper = probe;
      probe /= 2.0;
    }
    lower = probe;
  }

  constexpr double relative_width = 0.01;
  while (upper > lower && upper - lower > upper * relative_width) {
    const double middle = lower + (upper - lower) / 2.0;
    if (middle == lower || middle == upper) {
      break;
    }
    if (decide(middle)) {
      upper = middle;
    } else {
      lower = middle;
    }
  }

  std::vector<double> radii = annulus_candidates(first, second, lower, upper, stats);
  std::size_t begin = 0;
  std::size_t end = radii.size();
  while (begin < end) {
    const std::size_t middle = begin + (end - begin) / 2;
    if (decide(radii[middle])) {
      end = middle;
    } else {
      begin = middle + 1;
    }
  }
  if (begin < radii.size()) {
    return radii[begin];
  }

  // Defensive exact fallback: a valid infeasible/feasible bracket must contain an exact
  // candidate, but retain a full-candidate path so a future range-index regression cannot
  // silently turn the public exact result into the approximate upper endpoint.
  radii.clear();
  radii.reserve(1 + first.finite_points().size() + second.finite_points().size() +
                first.finite_points().size() * second.finite_points().size());
  radii.push_back(0.0);
  radii.insert(radii.end(), first.finite_diagonal_distances().begin(),
               first.finite_diagonal_distances().end());
  radii.insert(radii.end(), second.finite_diagonal_distances().begin(),
               second.finite_diagonal_distances().end());
  for (const Point& point : first.finite_points()) {
    for (const Point& other : second.finite_points()) {
      radii.push_back(point_distance(point, other));
    }
  }
  std::sort(radii.begin(), radii.end());
  radii.erase(std::unique(radii.begin(), radii.end()), radii.end());
  begin = 0;
  end = radii.size();
  while (begin < end) {
    const std::size_t middle = begin + (end - begin) / 2;
    if (decide(radii[middle])) {
      end = middle;
    } else {
      begin = middle + 1;
    }
  }
  return radii[begin];
}

}  // namespace bottleneck::detail
