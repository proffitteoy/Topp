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
  std::size_t min_size = 8;
  std::size_t max_size = 512;
  std::string pattern;
  std::string metric;
  std::vector<std::string> experiments;
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
  if (pattern == "duplicate_light") {
    auto first = random_diagram(generator, rows, -4.0, 4.0, 0.01, 3.0);
    auto second = random_diagram(generator, columns, -4.0, 4.0, 0.01, 3.0);
    if (first.size() > 1) {
      first.back() = first.front();
    }
    if (second.size() > 1) {
      second.back() = second.front();
    }
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
  if (pattern == "multi_component") {
    constexpr std::size_t component_size = 32;
    bottleneck::Diagram first;
    bottleneck::Diagram second;
    first.reserve(rows);
    second.reserve(columns);
    for (std::size_t index = 0; index < rows; ++index) {
      const double group = static_cast<double>(index / component_size) * 16.0;
      const double offset = static_cast<double>(index % component_size) * 0.01;
      first.push_back({group + offset, group + offset + 2.0});
    }
    for (std::size_t index = 0; index < columns; ++index) {
      const double group = static_cast<double>(index / component_size) * 16.0;
      const double offset = static_cast<double>(index % component_size) * 0.01;
      second.push_back({group + offset + 0.003, group + offset + 2.003});
    }
    return {std::move(first), std::move(second)};
  }
  return {random_diagram(generator, rows, -4.0, 4.0, 0.01, 3.0),
          random_diagram(generator, columns, -4.0, 4.0, 0.01, 3.0)};
}

void append_experiments(std::vector<std::string>& experiments,
                        std::string_view value) {
  std::size_t begin = 0;
  while (begin <= value.size()) {
    const std::size_t end = value.find(',', begin);
    const std::string_view name = value.substr(begin, end - begin);
    if (!name.empty()) {
      experiments.emplace_back(name);
    }
    if (end == std::string_view::npos) {
      break;
    }
    begin = end + 1;
  }
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index + 1 < argc; index += 2) {
    const std::string_view name(argv[index]);
    if (name == "--repetitions") {
      options.repetitions = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--rounds") {
      options.rounds = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--min-size") {
      options.min_size = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--max-size") {
      options.max_size = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--pattern") {
      options.pattern = argv[index + 1];
    } else if (name == "--metric") {
      options.metric = argv[index + 1];
    } else if (name == "--experiment" || name == "--experiments") {
      append_experiments(options.experiments, argv[index + 1]);
    }
  }
  return options;
}

bool selected_experiment(const Options& options, std::string_view name) {
  return options.experiments.empty() ||
         std::find(options.experiments.begin(), options.experiments.end(), name) !=
             options.experiments.end();
}

bool near(double first, double second) {
  const double scale = (std::max)({1.0, std::fabs(first), std::fabs(second)});
  return std::fabs(first - second) <= 2e-12 * scale;
}

std::vector<NamedConfig> experiment_configs(WassersteinMetric metric) {
  std::vector<NamedConfig> configs{
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
      {"row_reduced_dense_sap",
       {metric, WassersteinCandidateStrategy::dense_scalar,
        WassersteinGraphStrategy::dense_matrix,
        WassersteinMatcherStrategy::dense_sap_row_reduction,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"jv_reduced_dense_sap",
       {metric, WassersteinCandidateStrategy::dense_scalar,
        WassersteinGraphStrategy::dense_matrix,
        WassersteinMatcherStrategy::dense_sap_jv_reduction,
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
      {"parallel_dense",
       {metric, WassersteinCandidateStrategy::dense_parallel,
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
      {"arena_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr,
        WassersteinMatcherStrategy::sparse_sap_arena,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"fixed4_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::fixed_degree_4,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"fixed8_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::fixed_degree_8,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"fixed16_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::fixed_degree_16,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"fixed32_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::fixed_degree_32,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"bitmask_lazy_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::bitmask_lazy,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"bitmask_lazy_tiny",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::bitmask_lazy,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::tiny_sparse, WassersteinWarmStart::none}},
      {"block16_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::block_sparse_16,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"block32_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::block_sparse_32,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none}},
      {"block16_tiny",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::block_sparse_16,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::tiny_sparse, WassersteinWarmStart::none}},
      {"block32_tiny",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::block_sparse_32,
        WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::tiny_sparse, WassersteinWarmStart::none}},
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
      {"parallel_component_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::parallel_sparse,
        WassersteinWarmStart::none}},
      {"parallel_component_dense",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::dense_sap,
        WassersteinComponentStrategy::parallel_dense,
        WassersteinWarmStart::none}},
      {"greedy_warm_sparse",
       {metric, WassersteinCandidateStrategy::sweep_binary,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::tiny_sparse,
        WassersteinWarmStart::global_descending}},
      {"priced_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_simd_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_simd,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_simd_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_simd,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_simd_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_simd,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_simd_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_simd,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_simd_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_simd,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_sweep_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_sweep_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_sweep_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_sweep_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_sweep_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_kdtree_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_kdtree_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_kdtree_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_kdtree_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_kdtree_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_kdtree_persistent_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_kdtree_persistent_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_kdtree_persistent_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_kdtree_persistent_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_kdtree_persistent_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_kdtree_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_adaptive_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_adaptive,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_adaptive_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_adaptive,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_adaptive_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_adaptive,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_adaptive_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_adaptive,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_adaptive_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_adaptive,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_incremental_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_incremental_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_incremental_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_incremental_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_incremental_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_persistent_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_persistent_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_persistent_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_persistent_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_persistent_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_persistent,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_dynamic_topk2",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_dynamic,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_dynamic_topk4",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_dynamic,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_dynamic_topk8",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_dynamic,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_dynamic_topk16",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_dynamic,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_dynamic_topk32",
       {metric, WassersteinCandidateStrategy::topk_pricing_sweep_dynamic,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"priced_dynamic_batched_topk2",
       {metric,
        WassersteinCandidateStrategy::topk_pricing_sweep_dynamic_batched,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 2}},
      {"priced_dynamic_batched_topk4",
       {metric,
        WassersteinCandidateStrategy::topk_pricing_sweep_dynamic_batched,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 4}},
      {"priced_dynamic_batched_topk8",
       {metric,
        WassersteinCandidateStrategy::topk_pricing_sweep_dynamic_batched,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 8}},
      {"priced_dynamic_batched_topk16",
       {metric,
        WassersteinCandidateStrategy::topk_pricing_sweep_dynamic_batched,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 16}},
      {"priced_dynamic_batched_topk32",
       {metric,
        WassersteinCandidateStrategy::topk_pricing_sweep_dynamic_batched,
        WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none,
        bottleneck::WassersteinDuplicateStrategy::none, 32}},
      {"adaptive_no_duplicates",
       {metric, WassersteinCandidateStrategy::adaptive,
        WassersteinGraphStrategy::adaptive, WassersteinMatcherStrategy::adaptive,
        WassersteinComponentStrategy::adaptive,
        WassersteinWarmStart::global_descending,
        bottleneck::WassersteinDuplicateStrategy::none}},
      {"duplicate_exact",
       {metric, WassersteinCandidateStrategy::adaptive,
        WassersteinGraphStrategy::adaptive, WassersteinMatcherStrategy::adaptive,
        WassersteinComponentStrategy::adaptive,
        WassersteinWarmStart::global_descending,
        bottleneck::WassersteinDuplicateStrategy::exact}},
      {"duplicate_adaptive",
       {metric, WassersteinCandidateStrategy::adaptive,
        WassersteinGraphStrategy::adaptive, WassersteinMatcherStrategy::adaptive,
        WassersteinComponentStrategy::adaptive,
        WassersteinWarmStart::global_descending,
        bottleneck::WassersteinDuplicateStrategy::adaptive}},
      {"adaptive",
       {metric, WassersteinCandidateStrategy::adaptive,
        WassersteinGraphStrategy::adaptive, WassersteinMatcherStrategy::adaptive,
        WassersteinComponentStrategy::adaptive,
        WassersteinWarmStart::global_descending}},
  };
  for (NamedConfig& named : configs) {
    if (named.name != "duplicate_exact" && named.name != "duplicate_adaptive" &&
        named.name != "adaptive") {
      named.config.duplicates =
          bottleneck::WassersteinDuplicateStrategy::none;
    }
  }
  return configs;
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
  if (named.config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_full_scan ||
      named.config.candidates == WassersteinCandidateStrategy::topk_pricing_simd ||
      named.config.candidates == WassersteinCandidateStrategy::topk_pricing_sweep ||
      named.config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_kdtree ||
      named.config.candidates == WassersteinCandidateStrategy::
                                         topk_pricing_kdtree_persistent ||
      named.config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_adaptive ||
      named.config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_sweep_incremental ||
      named.config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_sweep_persistent ||
      named.config.candidates ==
          WassersteinCandidateStrategy::topk_pricing_sweep_dynamic ||
      named.config.candidates == WassersteinCandidateStrategy::
                                     topk_pricing_sweep_dynamic_batched) {
    return true;
  }
  if (pattern == "multi_component" &&
      (named.name == "component_dense" ||
       named.name == "parallel_component_dense" ||
       named.name == "parallel_component_sparse" ||
       named.name == "adaptive")) {
    return true;
  }
  const bool explicitly_dense =
      named.config.matcher == WassersteinMatcherStrategy::dense_hungarian ||
      named.config.matcher == WassersteinMatcherStrategy::dense_sap ||
      named.config.matcher ==
          WassersteinMatcherStrategy::dense_sap_row_reduction ||
      named.config.matcher ==
          WassersteinMatcherStrategy::dense_sap_jv_reduction ||
      named.config.components == WassersteinComponentStrategy::dense;
  if (named.config.graph == WassersteinGraphStrategy::bitmask_lazy &&
      maximum > 256) {
    return false;
  }
  if ((named.config.graph == WassersteinGraphStrategy::block_sparse_16 ||
       named.config.graph == WassersteinGraphStrategy::block_sparse_32) &&
      maximum > 1024) {
    return false;
  }
  if (explicitly_dense && maximum > 512) {
    return false;
  }
  if (maximum > 1024 && pattern != "separated" && pattern != "adversarial_sparse") {
    return named.name == "adaptive";
  }
  if (maximum > 1024) {
    return named.name == "sweep_csr_sparse" || named.name == "arena_sparse" ||
           named.name == "fixed4_sparse" || named.name == "fixed8_sparse" ||
           named.name == "fixed16_sparse" || named.name == "fixed32_sparse" ||
           named.name == "component_sparse" ||
           named.name == "tiny_component_sparse" ||
           named.name == "greedy_warm_sparse" || named.name == "adaptive";
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  const Options options = parse_options(argc, argv);
  const bool pricing_requested =
      std::any_of(options.experiments.begin(), options.experiments.end(),
                  [](const std::string& name) {
                    return name.starts_with("priced_topk") ||
                           name.starts_with("priced_simd_topk") ||
                           name.starts_with("priced_sweep_topk") ||
                           name.starts_with("priced_kdtree_topk") ||
                           name.starts_with("priced_kdtree_persistent_topk") ||
                           name.starts_with("priced_adaptive_topk") ||
                           name.starts_with("priced_incremental_topk") ||
                           name.starts_with("priced_persistent_topk") ||
                           name.starts_with("priced_dynamic_topk") ||
                           name.starts_with("priced_dynamic_batched_topk");
                  });
  std::mt19937_64 generator(0xD5A0A11ULL);
  const std::vector<std::size_t> all_sizes{8,   16,   32,   64,   128, 256,
                                           512, 1024, 2048, 4096, 8192};
  const std::vector<std::string_view> patterns{
      "uniform",          "near_diagonal",    "clustered",
      "separated",        "duplicate_heavy",  "duplicate_light",
      "imbalanced",
      "adversarial_dense", "adversarial_sparse", "multi_component"};
  const std::vector<WassersteinMetric> metrics{WassersteinMetric::w1_linf,
                                                WassersteinMetric::w2_l2};

  std::cout
      << "pattern,rows,columns,metric,experiment,repetitions,rounds,median_us,p95_us,"
         "prepare_us,candidate_us,graph_us,component_us,solver_us,pricing_us,"
         "candidate_density,"
         "edge_density,average_degree,max_degree,max_component,components,augmentations,"
         "pricing_rounds,priced_edges,pricing_violations,peak_materialized_edges,"
         "pricing_full_rounds,pricing_simd_rounds,pricing_sweep_rounds,"
         "pricing_kdtree_rounds,pricing_kdtree_builds,pricing_kdtree_updates,"
         "sparse_fallbacks,sparse_arena_builds,sparse_scratch_reuses,"
         "sparse_heap_growths,peak_sparse_arena_bytes,jv_fallbacks,greedy_matches,"
         "dynamic_inserted_edges,dynamic_dijkstra_runs,dynamic_batch_groups,"
         "warm_certificates,graph_bytes,"
         "peak_graph_bytes,"
         "duplicate_groups,duplicate_points_removed\n";
  double sink = 0.0;
  for (std::string_view pattern : patterns) {
    if (!options.pattern.empty() && pattern != options.pattern) {
      continue;
    }
    for (std::size_t size : all_sizes) {
      if (size < options.min_size || size > options.max_size) {
        continue;
      }
      const std::size_t rows = size;
      const std::size_t columns = pattern == "imbalanced" ? size * 4 : size;
      if (columns > options.max_size * 4) {
        continue;
      }
      const std::size_t maximum_dimension = (std::max)(rows, columns);
      if (maximum_dimension > 512 && pattern != "separated" &&
          pattern != "adversarial_sparse" && pattern != "multi_component" &&
          !pricing_requested) {
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
        if (!options.metric.empty() && bottleneck::to_string(metric) != options.metric) {
          continue;
        }
        const std::vector<NamedConfig> configs = experiment_configs(metric);
        const auto adaptive_reference = std::find_if(
            configs.begin(), configs.end(),
            [](const NamedConfig& named) { return named.name == "adaptive"; });
        const auto priced_reference = std::find_if(
            configs.begin(), configs.end(),
            [](const NamedConfig& named) { return named.name == "priced_topk32"; });
        const bool ordinary_large = maximum_dimension > 512 &&
                                    pattern != "separated" &&
                                    pattern != "adversarial_sparse" &&
                                    pattern != "multi_component";
        const std::size_t reference_index = maximum_dimension <= 512
                                                ? 0
                                            : ordinary_large
                                                ? static_cast<std::size_t>(
                                                      priced_reference -
                                                      configs.begin())
                                                : static_cast<std::size_t>(
                                                      adaptive_reference -
                                                      configs.begin());
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
            if (!selected_experiment(options, named.name)) {
              continue;
            }
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
                          << ", experiment=" << named.name
                          << ", expected=" << std::setprecision(17) << reference[index]
                          << ", actual=" << result << '\n';
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
            ADD_STAT(pricing_rounds);
            ADD_STAT(pricing_full_scan_rounds);
            ADD_STAT(pricing_simd_rounds);
            ADD_STAT(pricing_sweep_rounds);
            ADD_STAT(pricing_kdtree_rounds);
            ADD_STAT(pricing_kdtree_builds);
            ADD_STAT(pricing_kdtree_updates);
            ADD_STAT(priced_edges);
            ADD_STAT(pricing_violations);
            ADD_STAT(sparse_fallbacks);
            ADD_STAT(sparse_arena_builds);
            ADD_STAT(sparse_scratch_reuses);
            ADD_STAT(sparse_heap_growths);
            ADD_STAT(dynamic_inserted_edges);
            ADD_STAT(dynamic_dijkstra_runs);
            ADD_STAT(dynamic_batch_groups);
            ADD_STAT(jv_fallbacks);
            ADD_STAT(component_count);
            ADD_STAT(tiny_components);
            ADD_STAT(greedy_matches);
            ADD_STAT(warm_start_certificates);
            ADD_STAT(duplicate_groups);
            ADD_STAT(duplicate_points_removed);
            ADD_STAT(graph_bytes);
            ADD_STAT(candidate_time_ns);
            ADD_STAT(graph_time_ns);
            ADD_STAT(component_time_ns);
            ADD_STAT(solver_time_ns);
            ADD_STAT(pricing_time_ns);
#undef ADD_STAT
            total.largest_component =
                (std::max)(total.largest_component, round_stats.largest_component);
            total.peak_materialized_edges =
                (std::max)(total.peak_materialized_edges,
                           round_stats.peak_materialized_edges);
            total.max_degree =
                (std::max)(total.max_degree, round_stats.max_degree);
            total.peak_graph_bytes =
                (std::max)(total.peak_graph_bytes, round_stats.peak_graph_bytes);
            total.peak_sparse_arena_bytes =
                (std::max)(total.peak_sparse_arena_bytes,
                           round_stats.peak_sparse_arena_bytes);
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
                    << ',' << static_cast<double>(stats.pricing_time_ns) /
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
                    << ',' << stats.max_degree << ',' << stats.largest_component << ','
                    << stats.component_count
                    << ',' << stats.augmentations << ',' << stats.pricing_rounds << ','
                    << stats.priced_edges << ',' << stats.pricing_violations << ','
                    << stats.peak_materialized_edges << ','
                    << stats.pricing_full_scan_rounds << ','
                    << stats.pricing_simd_rounds << ','
                    << stats.pricing_sweep_rounds << ','
                    << stats.pricing_kdtree_rounds << ','
                    << stats.pricing_kdtree_builds << ','
                    << stats.pricing_kdtree_updates << ','
                    << stats.sparse_fallbacks << ',' << stats.sparse_arena_builds << ','
                    << stats.sparse_scratch_reuses << ','
                    << stats.sparse_heap_growths << ','
                    << stats.peak_sparse_arena_bytes << ','
                    << stats.jv_fallbacks << ',' << stats.greedy_matches << ','
                    << stats.dynamic_inserted_edges << ','
                    << stats.dynamic_dijkstra_runs << ','
                    << stats.dynamic_batch_groups << ','
                    << stats.warm_start_certificates << ','
                    << static_cast<double>(stats.graph_bytes) / call_count << ','
                    << stats.peak_graph_bytes << ','
                    << static_cast<double>(stats.duplicate_groups) / call_count << ','
                    << static_cast<double>(stats.duplicate_points_removed) / call_count
                    << '\n';
        }
      }
    }
  }
  if (sink < 0.0) {
    std::cerr << "unreachable sink value\n";
  }
  return 0;
}
