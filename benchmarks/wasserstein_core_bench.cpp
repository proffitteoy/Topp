#include <bottleneck/wasserstein.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using bottleneck::WassersteinCandidateStrategy;
using bottleneck::WassersteinComponentStrategy;
using bottleneck::WassersteinConfig;
using bottleneck::WassersteinGraphStrategy;
using bottleneck::WassersteinMatcherStrategy;
using bottleneck::WassersteinMetric;
using bottleneck::WassersteinWarmStart;

struct Options {
  std::size_t repetitions = 10;
  std::size_t rounds = 5;
  std::size_t max_size = 512;
};

struct NamedConfig {
  std::string_view name;
  WassersteinConfig config;
};

struct PreparedPair {
  bottleneck::PreparedDiagram first;
  bottleneck::PreparedDiagram second;
};

bottleneck::Diagram random_diagram(std::mt19937_64& generator, std::size_t size,
                                   double birth_min, double birth_max,
                                   double persistence_min, double persistence_max) {
  std::uniform_real_distribution<double> birth(birth_min, birth_max);
  std::uniform_real_distribution<double> persistence(persistence_min, persistence_max);
  bottleneck::Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    const double value = birth(generator);
    result.push_back({value, value + persistence(generator)});
  }
  return result;
}

std::pair<bottleneck::Diagram, bottleneck::Diagram> random_pair(
    std::mt19937_64& generator, std::size_t rows, std::size_t columns,
    std::string_view pattern) {
  if (pattern == "near_diagonal") {
    return {random_diagram(generator, rows, -4.0, 4.0, 0.001, 0.05),
            random_diagram(generator, columns, -4.0, 4.0, 0.001, 0.05)};
  }
  if (pattern == "clustered") {
    return {random_diagram(generator, rows, -0.1, 0.1, 1.0, 3.0),
            random_diagram(generator, columns, -0.1, 0.1, 1.0, 3.0)};
  }
  if (pattern == "separated") {
    return {random_diagram(generator, rows, -8.0, -4.0, 0.1, 1.0),
            random_diagram(generator, columns, 4.0, 8.0, 0.1, 1.0)};
  }
  if (pattern == "duplicate_heavy") {
    const std::vector<bottleneck::Point> prototypes{{-2.0, -0.5}, {-0.5, 1.0},
                                                     {0.0, 2.0}, {1.0, 3.5}};
    bottleneck::Diagram first;
    bottleneck::Diagram second;
    first.reserve(rows);
    second.reserve(columns);
    for (std::size_t index = 0; index < rows; ++index) {
      first.push_back(prototypes[index % prototypes.size()]);
    }
    for (std::size_t index = 0; index < columns; ++index) {
      second.push_back(prototypes[(index + 1) % prototypes.size()]);
    }
    std::shuffle(first.begin(), first.end(), generator);
    std::shuffle(second.begin(), second.end(), generator);
    return {std::move(first), std::move(second)};
  }
  if (pattern == "adversarial_dense") {
    return {random_diagram(generator, rows, -0.01, 0.01, 1.99, 2.01),
            random_diagram(generator, columns, -0.01, 0.01, 1.99, 2.01)};
  }
  if (pattern == "adversarial_sparse") {
    bottleneck::Diagram first;
    bottleneck::Diagram second;
    first.reserve(rows);
    second.reserve(columns);
    for (std::size_t index = 0; index < rows; ++index) {
      const double birth = static_cast<double>(index) * 4.0;
      first.push_back({birth, birth + 1.0});
    }
    for (std::size_t index = 0; index < columns; ++index) {
      const double birth = static_cast<double>(index) * 4.0 + 0.05;
      second.push_back({birth, birth + 1.0});
    }
    return {std::move(first), std::move(second)};
  }
  return {random_diagram(generator, rows, -4.0, 4.0, 0.01, 3.0),
          random_diagram(generator, columns, -4.0, 4.0, 0.01, 3.0)};
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index + 1 < argc; index += 2) {
    const std::string_view name(argv[index]);
    const std::size_t value = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    if (name == "--repetitions") {
      options.repetitions = value;
    } else if (name == "--rounds") {
      options.rounds = value;
    } else if (name == "--max-size") {
      options.max_size = value;
    }
  }
  return options;
}

bool near(double first, double second) {
  const double scale = (std::max)({1.0, std::fabs(first), std::fabs(second)});
  return std::fabs(first - second) <= 2e-12 * scale;
}

std::vector<NamedConfig> experiment_configs(WassersteinMetric metric) {
  return {
      {"dense_hungarian",
       {metric, WassersteinCandidateStrategy::dense_scalar,
        WassersteinGraphStrategy::dense_matrix,
        WassersteinMatcherStrategy::dense_hungarian,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"dense_sap",
       {metric, WassersteinCandidateStrategy::dense_scalar,
        WassersteinGraphStrategy::dense_matrix,
        WassersteinMatcherStrategy::dense_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"blocked_dense",
       {metric, WassersteinCandidateStrategy::dense_blocked,
        WassersteinGraphStrategy::dense_matrix,
        WassersteinMatcherStrategy::dense_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"avx2_dense",
       {metric, WassersteinCandidateStrategy::dense_avx2,
        WassersteinGraphStrategy::dense_matrix,
        WassersteinMatcherStrategy::dense_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"sweep_dense",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::dense_matrix,
        WassersteinMatcherStrategy::dense_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"dsr_csr_sparse",
       {metric, WassersteinCandidateStrategy::dense_scalar,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"sweep_csr_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"two_pointer_rows_sparse",
       {metric, WassersteinCandidateStrategy::sweep_two_pointer,
        WassersteinGraphStrategy::row_vectors,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"component_dense",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::dense_sap,
        WassersteinComponentStrategy::dense, WassersteinWarmStart::none}},
      {"component_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::sparse, WassersteinWarmStart::none}},
      {"tiny_component_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::tiny_sparse, WassersteinWarmStart::none}},
      {"greedy_warm_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::tiny_sparse,
        WassersteinWarmStart::global_descending}},
      {"adaptive",
       {metric, WassersteinCandidateStrategy::adaptive,
        WassersteinGraphStrategy::adaptive, WassersteinMatcherStrategy::adaptive,
        WassersteinComponentStrategy::adaptive,
        WassersteinWarmStart::global_descending}},
  };
}

double percentile(std::vector<double> values, double quantile) {
  std::sort(values.begin(), values.end());
  const std::size_t index = static_cast<std::size_t>(
      std::ceil(quantile * static_cast<double>(values.size()))) - 1;
  return values[(std::min)(index, values.size() - 1)];
}

bool should_run(const NamedConfig& named, std::size_t rows, std::size_t columns,
                std::string_view pattern) {
  const std::size_t maximum = (std::max)(rows, columns);
  const bool explicitly_dense =
      named.config.matcher == WassersteinMatcherStrategy::dense_hungarian ||
      named.config.matcher == WassersteinMatcherStrategy::dense_sap ||
      named.config.components == WassersteinComponentStrategy::dense;
  if (explicitly_dense && maximum > 512) {
    return false;
  }
  if (maximum > 1024 && pattern != "separated" && pattern != "adversarial_sparse") {
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  const Options options = parse_options(argc, argv);
  std::mt19937_64 generator(0xD5A0A11ULL);
  const std::vector<std::size_t> all_sizes{8,   16,   32,   64,   128, 256,
                                           512, 1024, 2048, 4096, 8192};
  const std::vector<std::string_view> patterns{
      "uniform",          "near_diagonal",    "clustered",
      "separated",        "duplicate_heavy",  "imbalanced",
      "adversarial_dense", "adversarial_sparse"};
  const std::vector<WassersteinMetric> metrics{WassersteinMetric::w1_linf,
                                                WassersteinMetric::w2_l2};

  std::cout
      << "pattern,rows,columns,metric,experiment,repetitions,rounds,median_us,p95_us,"
         "prepare_us,candidate_us,graph_us,component_us,solver_us,candidate_density,"
         "edge_density,average_degree,max_component,components,augmentations,"
         "greedy_matches,warm_certificates,graph_bytes\n";
  double sink = 0.0;
  for (std::string_view pattern : patterns) {
    for (std::size_t size : all_sizes) {
      if (size > options.max_size) {
        continue;
      }
      const std::size_t rows = size;
      const std::size_t columns = pattern == "imbalanced" ? size * 4 : size;
      if (columns > options.max_size * 4) {
        continue;
      }
      const std::size_t maximum_dimension = (std::max)(rows, columns);
      if (maximum_dimension > 512 && pattern != "separated" &&
          pattern != "adversarial_sparse") {
        continue;
      }
      std::vector<std::pair<bottleneck::Diagram, bottleneck::Diagram>> raw_pairs;
      raw_pairs.reserve(options.repetitions);
      for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
        raw_pairs.push_back(random_pair(generator, rows, columns, pattern));
      }
      const auto prepare_start = Clock::now();
      std::vector<PreparedPair> pairs;
      pairs.reserve(raw_pairs.size());
      for (const auto& pair : raw_pairs) {
        pairs.push_back(
            {bottleneck::PreparedDiagram(pair.first),
             bottleneck::PreparedDiagram(pair.second)});
      }
      const auto prepare_stop = Clock::now();
      const double prepare_us =
          std::chrono::duration<double, std::micro>(prepare_stop - prepare_start).count() /
          static_cast<double>(pairs.size());

      for (WassersteinMetric metric : metrics) {
        const std::vector<NamedConfig> configs = experiment_configs(metric);
        const std::size_t reference_index = maximum_dimension <= 512 ? 0 : 6;
        std::vector<double> reference;
        reference.reserve(pairs.size());
        for (const PreparedPair& pair : pairs) {
          reference.push_back(bottleneck::wasserstein_distance(
              pair.first, pair.second, configs[reference_index].config));
        }

        std::vector<std::size_t> order(configs.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::vector<std::vector<double>> timings(configs.size());
        std::vector<bottleneck::WassersteinStats> totals(configs.size());
        std::vector<std::size_t> calls(configs.size(), 0);
        for (std::size_t round = 0; round < options.rounds; ++round) {
          std::shuffle(order.begin(), order.end(), generator);
          for (std::size_t config_index : order) {
            const NamedConfig& named = configs[config_index];
            if (!should_run(named, rows, columns, pattern)) {
              continue;
            }
            bottleneck::WassersteinStats round_stats;
            const auto start = Clock::now();
            for (std::size_t index = 0; index < pairs.size(); ++index) {
              const double result = bottleneck::wasserstein_distance(
                  pairs[index].first, pairs[index].second, named.config,
                  &round_stats);
              if (!near(result, reference[index])) {
                std::cerr << "result mismatch: pattern=" << pattern
                          << ", rows=" << rows << ", columns=" << columns
                          << ", metric=" << bottleneck::to_string(metric)
                          << ", experiment=" << named.name << '\n';
                return 2;
              }
              sink += result;
            }
            const auto stop = Clock::now();
            timings[config_index].push_back(
                std::chrono::duration<double, std::micro>(stop - start).count() /
                static_cast<double>(pairs.size()));
            auto& total = totals[config_index];
#define ADD_STAT(name) total.name += round_stats.name
            ADD_STAT(possible_pairs);
            ADD_STAT(candidate_pairs);
            ADD_STAT(positive_edges);
            ADD_STAT(pruned_pairs);
            ADD_STAT(active_rows);
            ADD_STAT(active_columns);
            ADD_STAT(augmentations);
            ADD_STAT(component_count);
            ADD_STAT(tiny_components);
            ADD_STAT(greedy_matches);
            ADD_STAT(warm_start_certificates);
            ADD_STAT(graph_bytes);
            ADD_STAT(candidate_time_ns);
            ADD_STAT(graph_time_ns);
            ADD_STAT(component_time_ns);
            ADD_STAT(solver_time_ns);
#undef ADD_STAT
            total.largest_component =
                (std::max)(total.largest_component, round_stats.largest_component);
            calls[config_index] += pairs.size();
          }
        }

        for (std::size_t config_index = 0; config_index < configs.size(); ++config_index) {
          if (timings[config_index].empty()) {
            continue;
          }
          const auto& stats = totals[config_index];
          const double call_count = static_cast<double>(calls[config_index]);
          const double possible = static_cast<double>(stats.possible_pairs);
          std::cout << pattern << ',' << rows << ',' << columns << ','
                    << bottleneck::to_string(metric) << ',' << configs[config_index].name
                    << ',' << options.repetitions << ',' << options.rounds << ','
                    << std::fixed << std::setprecision(3)
                    << percentile(timings[config_index], 0.5) << ','
                    << percentile(timings[config_index], 0.95) << ',' << prepare_us
                    << ',' << static_cast<double>(stats.candidate_time_ns) /
                                   (1000.0 * call_count)
                    << ',' << static_cast<double>(stats.graph_time_ns) /
                                   (1000.0 * call_count)
                    << ',' << static_cast<double>(stats.component_time_ns) /
                                   (1000.0 * call_count)
                    << ',' << static_cast<double>(stats.solver_time_ns) /
                                   (1000.0 * call_count)
                    << ',' << (possible == 0.0
                                   ? 0.0
                                   : static_cast<double>(stats.candidate_pairs) / possible)
                    << ',' << (possible == 0.0
                                   ? 0.0
                                   : static_cast<double>(stats.positive_edges) / possible)
                    << ',' << (stats.active_rows == 0
                                   ? 0.0
                                   : static_cast<double>(stats.positive_edges) /
                                         static_cast<double>(stats.active_rows))
                    << ',' << stats.largest_component << ',' << stats.component_count
                    << ',' << stats.augmentations << ',' << stats.greedy_matches << ','
                    << stats.warm_start_certificates << ','
                    << static_cast<double>(stats.graph_bytes) / call_count << '\n';
        }
      }
    }
  }
  if (sink < 0.0) {
    std::cerr << "unreachable sink value\n";
  }
  return 0;
}
