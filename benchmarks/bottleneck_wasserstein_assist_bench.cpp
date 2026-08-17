#include <bottleneck/core.hpp>
#include <bottleneck/wasserstein.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using bottleneck::Diagram;
using bottleneck::PreparedDiagram;

enum class Pattern {
  uniform,
  clustered,
  near_diagonal,
  duplicate_heavy,
  separated,
  threshold_shell
};

struct Options {
  std::size_t repetitions = 1;
  std::size_t rounds = 7;
  std::size_t max_points = 128;
  std::string_view pattern = "all";
};

struct Workload {
  Pattern pattern;
  std::size_t first_size;
  std::size_t second_size;
};

struct MatchingResult {
  double distance = 0.0;
  double upper = 0.0;
  std::size_t summands = 0;
  std::vector<std::pair<std::size_t, std::size_t>> cross_matches;
  std::vector<double> edge_costs;
};

struct GreedyResult {
  MatchingResult matching;
  double upper_at_five_percent = 0.0;
  double upper_at_ten_percent = 0.0;
};

struct WeightedEdge {
  std::size_t row = 0;
  std::size_t column = 0;
  long double saving = 0.0L;
};

const char* pattern_name(Pattern pattern) {
  switch (pattern) {
    case Pattern::uniform:
      return "uniform";
    case Pattern::clustered:
      return "clustered";
    case Pattern::near_diagonal:
      return "near_diagonal";
    case Pattern::duplicate_heavy:
      return "duplicate_heavy";
    case Pattern::separated:
      return "separated";
    case Pattern::threshold_shell:
      return "threshold_shell";
  }
  return "unknown";
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; index += 2) {
    if (index + 1 >= argc) {
      throw std::invalid_argument("missing benchmark option value");
    }
    const std::string_view name(argv[index]);
    if (name == "--repetitions") {
      options.repetitions = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--rounds") {
      options.rounds = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--max-points") {
      options.max_points = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--pattern") {
      options.pattern = argv[index + 1];
    } else {
      throw std::invalid_argument("unknown benchmark option");
    }
  }
  if (options.repetitions == 0 || options.rounds == 0 || options.max_points == 0) {
    throw std::invalid_argument("benchmark sizes must be positive");
  }
  return options;
}

Diagram generate_diagram(std::mt19937_64& generator, Pattern pattern, std::size_t size,
                         bool second) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 0.15);
  const double shell_radius = std::ldexp(1.0, -10);
  Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    double birth = 0.0;
    double persistence = 0.0;
    switch (pattern) {
      case Pattern::uniform:
        birth = -4.0 + 8.0 * unit(generator);
        persistence = 0.05 + 1.95 * unit(generator);
        break;
      case Pattern::clustered: {
        const double center = static_cast<double>(static_cast<int>(index % 4) - 2) * 1.5;
        birth = center + normal(generator) + (second ? 0.01 : 0.0);
        persistence = 0.4 + 0.8 * unit(generator);
        break;
      }
      case Pattern::near_diagonal:
        birth = -4.0 + 8.0 * unit(generator);
        persistence = 0.0001 + 0.04 * unit(generator);
        break;
      case Pattern::duplicate_heavy:
        birth = static_cast<double>(static_cast<int>(index % 8) - 4) * 0.5 +
                (second ? 0.01 : 0.0);
        persistence = 0.25 + 0.125 * static_cast<double>(index % 4);
        break;
      case Pattern::separated:
        birth = (second ? 10.0 : -12.0) + 2.0 * unit(generator);
        persistence = 0.05 + 1.95 * unit(generator);
        break;
      case Pattern::threshold_shell:
        birth = static_cast<double>(index) * shell_radius /
                static_cast<double>((std::max)(std::size_t{1}, size));
        persistence = 2.0 + (second ? shell_radius : 0.0) - birth;
        break;
    }
    result.push_back({birth, birth + persistence});
  }
  return result;
}

long double power_cost(double value, unsigned q) {
  long double result = 1.0L;
  const long double base = static_cast<long double>(value);
  for (unsigned exponent = 0; exponent < q; ++exponent) {
    result *= base;
  }
  return result;
}

double cross_cost(const bottleneck::Point& first, const bottleneck::Point& second) {
  return (std::max)(std::fabs(first.birth - second.birth),
                    std::fabs(first.death - second.death));
}

std::vector<long double> savings_matrix(const PreparedDiagram& first,
                                        const PreparedDiagram& second, unsigned q,
                                        std::size_t side) {
  std::vector<long double> weights(side * side, 0.0L);
  for (std::size_t row = 0; row < first.finite_points().size(); ++row) {
    for (std::size_t column = 0; column < second.finite_points().size(); ++column) {
      const long double saving =
          power_cost(first.finite_diagonal_distances()[row], q) +
          power_cost(second.finite_diagonal_distances()[column], q) -
          power_cost(cross_cost(first.finite_points()[row], second.finite_points()[column]), q);
      weights[row * side + column] = (std::max)(0.0L, saving);
    }
  }
  return weights;
}

MatchingResult materialize_matching(const PreparedDiagram& first,
                                    const PreparedDiagram& second, unsigned q,
                                    const std::vector<long double>& weights,
                                    const std::vector<std::size_t>& assignment,
                                    std::size_t side) {
  std::vector<std::uint8_t> matched_first(first.finite_points().size(), 0);
  std::vector<std::uint8_t> matched_second(second.finite_points().size(), 0);
  MatchingResult result;
  long double total_power = 0.0L;
  for (std::size_t row = 0; row < first.finite_points().size(); ++row) {
    const std::size_t column = assignment[row];
    if (column < second.finite_points().size() && weights[row * side + column] > 0.0L) {
      const double cost = cross_cost(first.finite_points()[row], second.finite_points()[column]);
      result.cross_matches.emplace_back(row, column);
      result.edge_costs.push_back(cost);
      total_power += power_cost(cost, q);
      matched_first[row] = 1;
      matched_second[column] = 1;
    }
  }
  for (std::size_t row = 0; row < matched_first.size(); ++row) {
    if (matched_first[row] == 0) {
      const double cost = first.finite_diagonal_distances()[row];
      result.edge_costs.push_back(cost);
      total_power += power_cost(cost, q);
    }
  }
  for (std::size_t column = 0; column < matched_second.size(); ++column) {
    if (matched_second[column] == 0) {
      const double cost = second.finite_diagonal_distances()[column];
      result.edge_costs.push_back(cost);
      total_power += power_cost(cost, q);
    }
  }
  result.summands = result.edge_costs.size();
  for (double cost : result.edge_costs) {
    result.upper = (std::max)(result.upper, cost);
  }
  result.distance = total_power <= 0.0L
                        ? 0.0
                        : static_cast<double>(std::pow(total_power, 1.0L / q));
  return result;
}

MatchingResult optimal_power_matching(const PreparedDiagram& first,
                                      const PreparedDiagram& second, unsigned q) {
  const std::size_t side = (std::max)(first.finite_points().size(),
                                      second.finite_points().size());
  if (side == 0) {
    return {};
  }
  const std::vector<long double> weights = savings_matrix(first, second, q, side);
  std::vector<long double> row_potential(side + 1, 0.0L);
  std::vector<long double> column_potential(side + 1, 0.0L);
  std::vector<std::size_t> matched_row(side + 1, 0);
  std::vector<std::size_t> predecessor(side + 1, 0);
  for (std::size_t row = 1; row <= side; ++row) {
    matched_row[0] = row;
    std::size_t column0 = 0;
    std::vector<long double> minimum(side + 1,
                                     std::numeric_limits<long double>::infinity());
    std::vector<std::uint8_t> used(side + 1, 0);
    do {
      used[column0] = 1;
      const std::size_t row0 = matched_row[column0];
      long double delta = std::numeric_limits<long double>::infinity();
      std::size_t column1 = 0;
      for (std::size_t column = 1; column <= side; ++column) {
        if (used[column] != 0) {
          continue;
        }
        const long double reduced =
            -weights[(row0 - 1) * side + (column - 1)] - row_potential[row0] -
            column_potential[column];
        if (reduced < minimum[column]) {
          minimum[column] = reduced;
          predecessor[column] = column0;
        }
        if (minimum[column] < delta) {
          delta = minimum[column];
          column1 = column;
        }
      }
      for (std::size_t column = 0; column <= side; ++column) {
        if (used[column] != 0) {
          row_potential[matched_row[column]] += delta;
          column_potential[column] -= delta;
        } else {
          minimum[column] -= delta;
        }
      }
      column0 = column1;
    } while (matched_row[column0] != 0);
    do {
      const std::size_t column1 = predecessor[column0];
      matched_row[column0] = matched_row[column1];
      column0 = column1;
    } while (column0 != 0);
  }
  std::vector<std::size_t> assignment(side, side);
  for (std::size_t column = 1; column <= side; ++column) {
    if (matched_row[column] != 0) {
      assignment[matched_row[column] - 1] = column - 1;
    }
  }
  return materialize_matching(first, second, q, weights, assignment, side);
}

double greedy_upper_from_prefix(const PreparedDiagram& first,
                                const PreparedDiagram& second,
                                const std::vector<WeightedEdge>& edges,
                                std::size_t prefix) {
  std::vector<std::uint8_t> matched_first(first.finite_points().size(), 0);
  std::vector<std::uint8_t> matched_second(second.finite_points().size(), 0);
  double upper = 0.0;
  for (std::size_t index = 0; index < prefix; ++index) {
    const WeightedEdge& edge = edges[index];
    if (matched_first[edge.row] == 0 && matched_second[edge.column] == 0) {
      matched_first[edge.row] = 1;
      matched_second[edge.column] = 1;
      upper = (std::max)(upper,
                         cross_cost(first.finite_points()[edge.row],
                                    second.finite_points()[edge.column]));
    }
  }
  for (std::size_t row = 0; row < matched_first.size(); ++row) {
    if (matched_first[row] == 0) {
      upper = (std::max)(upper, first.finite_diagonal_distances()[row]);
    }
  }
  for (std::size_t column = 0; column < matched_second.size(); ++column) {
    if (matched_second[column] == 0) {
      upper = (std::max)(upper, second.finite_diagonal_distances()[column]);
    }
  }
  return upper;
}

GreedyResult greedy_power_matching(const PreparedDiagram& first,
                                   const PreparedDiagram& second, unsigned q) {
  const std::size_t side = (std::max)(first.finite_points().size(),
                                      second.finite_points().size());
  const std::vector<long double> weights = savings_matrix(first, second, q, side);
  std::vector<WeightedEdge> edges;
  edges.reserve(first.finite_points().size() * second.finite_points().size());
  for (std::size_t row = 0; row < first.finite_points().size(); ++row) {
    for (std::size_t column = 0; column < second.finite_points().size(); ++column) {
      const long double saving = weights[row * side + column];
      if (saving > 0.0L) {
        edges.push_back({row, column, saving});
      }
    }
  }
  std::sort(edges.begin(), edges.end(), [](const WeightedEdge& left,
                                           const WeightedEdge& right) {
    return left.saving > right.saving;
  });
  std::vector<std::uint8_t> matched_first(first.finite_points().size(), 0);
  std::vector<std::uint8_t> matched_second(second.finite_points().size(), 0);
  std::vector<std::size_t> assignment(side, side);
  for (const WeightedEdge& edge : edges) {
    if (matched_first[edge.row] == 0 && matched_second[edge.column] == 0) {
      matched_first[edge.row] = 1;
      matched_second[edge.column] = 1;
      assignment[edge.row] = edge.column;
    }
  }
  GreedyResult result;
  result.matching = materialize_matching(first, second, q, weights, assignment, side);
  result.upper_at_five_percent = greedy_upper_from_prefix(
      first, second, edges, (edges.size() + 19) / 20);
  result.upper_at_ten_percent = greedy_upper_from_prefix(
      first, second, edges, (edges.size() + 9) / 10);
  return result;
}

double percentile(std::vector<double> values, double quantile) {
  std::sort(values.begin(), values.end());
  const std::size_t index = static_cast<std::size_t>(
      std::ceil(quantile * static_cast<double>(values.size())) - 1.0);
  return values[(std::min)(index, values.size() - 1)];
}

std::vector<double> bottleneck_candidates(const PreparedDiagram& first,
                                          const PreparedDiagram& second) {
  std::vector<double> result{0.0};
  result.insert(result.end(), first.finite_diagonal_distances().begin(),
                first.finite_diagonal_distances().end());
  result.insert(result.end(), second.finite_diagonal_distances().begin(),
                second.finite_diagonal_distances().end());
  for (const auto& left : first.finite_points()) {
    for (const auto& right : second.finite_points()) {
      result.push_back(cross_cost(left, right));
    }
  }
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

double candidate_reduction(const std::vector<double>& candidates, double lower, double upper) {
  const auto begin = std::lower_bound(candidates.begin(), candidates.end(), lower);
  const auto end = std::upper_bound(begin, candidates.end(), upper);
  return 1.0 - static_cast<double>(end - begin) /
                   static_cast<double>(candidates.size());
}

double matching_survival(const MatchingResult& matching, double threshold) {
  if (matching.edge_costs.empty()) {
    return 1.0;
  }
  const std::size_t surviving = static_cast<std::size_t>(
      std::count_if(matching.edge_costs.begin(), matching.edge_costs.end(),
                    [&](double cost) { return cost <= threshold; }));
  return static_cast<double>(surviving) /
         static_cast<double>(matching.edge_costs.size());
}

double cross_overlap(const MatchingResult& current, const MatchingResult& previous) {
  if (current.cross_matches.empty() && previous.cross_matches.empty()) {
    return 1.0;
  }
  std::size_t intersection = 0;
  for (const auto& edge : current.cross_matches) {
    intersection += static_cast<std::size_t>(
        std::find(previous.cross_matches.begin(), previous.cross_matches.end(), edge) !=
        previous.cross_matches.end());
  }
  return static_cast<double>(intersection) /
         static_cast<double>((std::max)(current.cross_matches.size(),
                                        previous.cross_matches.size()));
}

bool close(double first, double second) {
  const double scale = (std::max)({1.0, std::fabs(first), std::fabs(second)});
  return std::fabs(first - second) <= 1e-10 * scale;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    const std::vector<unsigned> exponents{1, 2, 4, 8, 16};
    const std::vector<Pattern> patterns{Pattern::uniform, Pattern::clustered,
                                        Pattern::near_diagonal, Pattern::duplicate_heavy,
                                        Pattern::separated, Pattern::threshold_shell};
    const std::vector<std::size_t> sizes{32, 64, 128};
    std::vector<Workload> workloads;
    for (Pattern pattern : patterns) {
      for (std::size_t size : sizes) {
        workloads.push_back({pattern, size, size});
      }
    }
    for (Pattern pattern : {Pattern::uniform, Pattern::clustered,
                            Pattern::near_diagonal}) {
      workloads.push_back({pattern, 32, 128});
      workloads.push_back({pattern, 128, 32});
    }
    std::mt19937_64 data_generator(0xE1E5A55157ULL);
    std::mt19937_64 order_generator(0xE1E5B077EULL);
    double sink = 0.0;

    std::cout << "pattern,n,m,q,bottleneck_distance,wq,lower,upper,lower_ratio,upper_ratio,"
                 "candidate_count,candidate_reduction,matching_survival,cross_overlap_previous_q,"
                 "greedy_upper_ratio,greedy_5pct_upper_ratio,greedy_10pct_upper_ratio,"
                 "bottleneck_median_us,wq_reference_median_us,wq_reference_p95_us,"
                 "optimized_w1_median_us,prepass_floor_over_bottleneck,"
                 "w1_lower_assisted_median_us,assisted_speedup,hint_accepted,"
                 "assisted_decisions\n";
    for (const Workload& workload : workloads) {
      if ((std::max)(workload.first_size, workload.second_size) > options.max_points ||
          (options.pattern != "all" && options.pattern != pattern_name(workload.pattern))) {
        continue;
      }
      const PreparedDiagram first(
          generate_diagram(data_generator, workload.pattern, workload.first_size, false));
      const PreparedDiagram second(
          generate_diagram(data_generator, workload.pattern, workload.second_size, true));
      const double bottleneck_distance = bottleneck::bottleneck_distance(first, second);
      const std::vector<double> candidates = bottleneck_candidates(first, second);
      std::vector<MatchingResult> matching_results;
      std::vector<GreedyResult> greedy_results;
      matching_results.reserve(exponents.size());
      greedy_results.reserve(exponents.size());
      for (unsigned q : exponents) {
        matching_results.push_back(optimal_power_matching(first, second, q));
        greedy_results.push_back(greedy_power_matching(first, second, q));
      }
      bottleneck::WassersteinConfig w1_config;
      w1_config.metric = bottleneck::WassersteinMetric::w1_linf;
      const double optimized_w1 = bottleneck::wasserstein_distance(first, second, w1_config);
      if (!close(optimized_w1, matching_results.front().distance)) {
        throw std::runtime_error("optimized W1 disagrees with dense reference");
      }

      std::vector<std::vector<double>> q_samples(exponents.size());
      std::vector<double> bottleneck_samples;
      std::vector<double> optimized_w1_samples;
      std::vector<double> assisted_samples;
      bottleneck::SolverStats assisted_stats;
      std::vector<std::size_t> order(exponents.size());
      std::iota(order.begin(), order.end(), std::size_t{0});
      for (std::size_t round = 0; round < options.rounds; ++round) {
        std::shuffle(order.begin(), order.end(), order_generator);
        for (std::size_t q_index : order) {
          const auto start = Clock::now();
          MatchingResult measured;
          for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
            measured = optimal_power_matching(first, second, exponents[q_index]);
            sink += measured.distance;
          }
          const auto stop = Clock::now();
          if (!close(measured.distance, matching_results[q_index].distance)) {
            throw std::runtime_error("dense q-reference is not deterministic");
          }
          q_samples[q_index].push_back(
              std::chrono::duration<double, std::micro>(stop - start).count() /
              static_cast<double>(options.repetitions));
        }
        auto start = Clock::now();
        for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
          sink += bottleneck::bottleneck_distance(first, second);
        }
        auto stop = Clock::now();
        bottleneck_samples.push_back(
            std::chrono::duration<double, std::micro>(stop - start).count() /
            static_cast<double>(options.repetitions));
        start = Clock::now();
        for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
          sink += bottleneck::wasserstein_distance(first, second, w1_config);
        }
        stop = Clock::now();
        optimized_w1_samples.push_back(
            std::chrono::duration<double, std::micro>(stop - start).count() /
            static_cast<double>(options.repetitions));

        start = Clock::now();
        for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
          const double w1 = bottleneck::wasserstein_distance(first, second, w1_config);
          bottleneck::SolverConfig assisted_config;
          const std::size_t safe_summands = first.finite_points().size() +
                                            second.finite_points().size();
          assisted_config.lower_bound_hint =
              safe_summands == 0 ? 0.0 : w1 / static_cast<double>(safe_summands);
          const double value = bottleneck::bottleneck_distance(
              first, second, assisted_config, &assisted_stats);
          if (value != bottleneck_distance) {
            throw std::runtime_error("W1 lower-assisted Bottleneck mismatch");
          }
          sink += value;
        }
        stop = Clock::now();
        assisted_samples.push_back(
            std::chrono::duration<double, std::micro>(stop - start).count() /
            static_cast<double>(options.repetitions));
      }

      const double bottleneck_us = percentile(bottleneck_samples, 0.5);
      const double optimized_w1_us = percentile(optimized_w1_samples, 0.5);
      const double assisted_us = percentile(assisted_samples, 0.5);
      for (std::size_t q_index = 0; q_index < exponents.size(); ++q_index) {
        const unsigned q = exponents[q_index];
        const MatchingResult& matching = matching_results[q_index];
        const GreedyResult& greedy = greedy_results[q_index];
        const double lower = matching.summands == 0
                                 ? 0.0
                                 : matching.distance /
                                       std::pow(static_cast<double>(matching.summands),
                                                1.0 / static_cast<double>(q));
        const double reference_us = percentile(q_samples[q_index], 0.5);
        std::cout << pattern_name(workload.pattern) << ',' << workload.first_size << ','
                  << workload.second_size << ',' << q << ',' << std::fixed
                  << std::setprecision(9) << bottleneck_distance << ','
                  << matching.distance << ',' << lower << ',' << matching.upper << ','
                  << (bottleneck_distance == 0.0 ? 1.0 : lower / bottleneck_distance) << ','
                  << (bottleneck_distance == 0.0 ? 1.0
                                                  : matching.upper / bottleneck_distance)
                  << ',' << candidates.size() << ','
                  << candidate_reduction(candidates, lower, matching.upper) << ','
                  << matching_survival(matching, bottleneck_distance) << ','
                  << (q_index == 0 ? 1.0
                                   : cross_overlap(matching, matching_results[q_index - 1]))
                  << ','
                  << (bottleneck_distance == 0.0
                          ? 1.0
                          : greedy.matching.upper / bottleneck_distance)
                  << ','
                  << (bottleneck_distance == 0.0
                          ? 1.0
                          : greedy.upper_at_five_percent / bottleneck_distance)
                  << ','
                  << (bottleneck_distance == 0.0
                          ? 1.0
                          : greedy.upper_at_ten_percent / bottleneck_distance)
                  << ',' << std::setprecision(3) << bottleneck_us << ',' << reference_us
                  << ',' << percentile(q_samples[q_index], 0.95) << ','
                  << (q == 1 ? optimized_w1_us : 0.0) << ','
                  << (q == 1 ? optimized_w1_us : reference_us) / bottleneck_us << ','
                  << assisted_us << ',' << bottleneck_us / assisted_us << ','
                  << assisted_stats.lower_bound_hints_accepted << ','
                  << static_cast<double>(assisted_stats.threshold_decisions) /
                         static_cast<double>(options.repetitions * options.rounds)
                  << '\n';
      }
    }
    if (sink == -1.0) {
      std::cerr << "unreachable sink\n";
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
