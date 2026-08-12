#include <bottleneck/wasserstein.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

using bottleneck::Diagram;
using bottleneck::Point;
using bottleneck::PreparedDiagram;
using bottleneck::WassersteinCandidateStrategy;
using bottleneck::WassersteinComponentStrategy;
using bottleneck::WassersteinConfig;
using bottleneck::WassersteinDuplicateStrategy;
using bottleneck::WassersteinGraphStrategy;
using bottleneck::WassersteinMatcherStrategy;
using bottleneck::WassersteinMetric;
using bottleneck::WassersteinWarmStart;

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void expect_near(double actual, double expected, const std::string& label,
                 double tolerance = 2e-12) {
  if (actual == expected) {
    return;
  }
  const double scale = (std::max)({1.0, std::fabs(actual), std::fabs(expected)});
  if (std::fabs(actual - expected) > tolerance * scale) {
    fail(label + ": expected=" + std::to_string(expected) +
         ", actual=" + std::to_string(actual));
  }
}

long double diagonal_power(const Point& point, WassersteinMetric metric) {
  const long double persistence =
      (static_cast<long double>(point.death) - point.birth) / 2.0L;
  return metric == WassersteinMetric::w1_linf ? persistence
                                               : 2.0L * persistence * persistence;
}

long double cross_power(const Point& first, const Point& second,
                        WassersteinMetric metric) {
  const long double birth =
      std::fabs(static_cast<long double>(first.birth) - second.birth);
  const long double death =
      std::fabs(static_cast<long double>(first.death) - second.death);
  return metric == WassersteinMetric::w1_linf ? (std::max)(birth, death)
                                               : birth * birth + death * death;
}

double brute_force_distance(const Diagram& first, const Diagram& second,
                            WassersteinMetric metric) {
  long double baseline = 0;
  for (const Point& point : first) {
    baseline += diagonal_power(point, metric);
  }
  for (const Point& point : second) {
    baseline += diagonal_power(point, metric);
  }

  long double best_saving = 0;
  std::vector<bool> used(second.size(), false);
  const std::function<void(std::size_t, long double)> search =
      [&](std::size_t index, long double saving) {
        if (index == first.size()) {
          best_saving = (std::max)(best_saving, saving);
          return;
        }
        search(index + 1, saving);
        for (std::size_t candidate = 0; candidate < second.size(); ++candidate) {
          if (used[candidate]) {
            continue;
          }
          const long double edge = diagonal_power(first[index], metric) +
                                   diagonal_power(second[candidate], metric) -
                                   cross_power(first[index], second[candidate], metric);
          if (edge <= 0) {
            continue;
          }
          used[candidate] = true;
          search(index + 1, saving + edge);
          used[candidate] = false;
        }
      };
  search(0, 0);
  const long double powered = (std::max)(0.0L, baseline - best_saving);
  return metric == WassersteinMetric::w1_linf
             ? static_cast<double>(powered)
             : static_cast<double>(std::sqrt(powered));
}

void deterministic_cases() {
  const Diagram empty;
  const Diagram point{{0.0, 2.0}};
  const WassersteinConfig w1_dense{
      WassersteinMetric::w1_linf, WassersteinCandidateStrategy::dense_scalar,
      WassersteinGraphStrategy::dense_matrix,
      WassersteinMatcherStrategy::dense_hungarian,
      WassersteinComponentStrategy::none, WassersteinWarmStart::none};
  const WassersteinConfig w1_sweep{
      WassersteinMetric::w1_linf, WassersteinCandidateStrategy::sweep_binary,
      WassersteinGraphStrategy::dense_matrix, WassersteinMatcherStrategy::dense_sap,
      WassersteinComponentStrategy::none, WassersteinWarmStart::none};
  const WassersteinConfig w2_dense{
      WassersteinMetric::w2_l2, WassersteinCandidateStrategy::dense_scalar,
      WassersteinGraphStrategy::dense_matrix,
      WassersteinMatcherStrategy::dense_hungarian,
      WassersteinComponentStrategy::none, WassersteinWarmStart::none};
  const WassersteinConfig w2_sweep{
      WassersteinMetric::w2_l2, WassersteinCandidateStrategy::sweep_binary,
      WassersteinGraphStrategy::dense_matrix, WassersteinMatcherStrategy::dense_sap,
      WassersteinComponentStrategy::none, WassersteinWarmStart::none};

  expect_near(bottleneck::wasserstein_distance(empty, empty, w1_sweep), 0.0,
              "empty W1");
  expect_near(bottleneck::wasserstein_distance(empty, point, w1_sweep), 1.0,
              "empty-to-point W1");
  expect_near(bottleneck::wasserstein_distance(empty, point, w2_sweep), std::sqrt(2.0),
              "empty-to-point W2");
  expect_near(bottleneck::wasserstein_distance(point, point, w1_sweep), 0.0,
              "identity W1");
  expect_near(bottleneck::wasserstein_distance(point, point, w2_sweep), 0.0,
              "identity W2");

  const Diagram shifted{{1.0, 3.0}};
  expect_near(bottleneck::wasserstein_distance(point, shifted, w1_sweep), 1.0,
              "shifted W1");
  expect_near(bottleneck::wasserstein_distance(point, shifted, w2_sweep), std::sqrt(2.0),
              "shifted W2");

  const Diagram separated{{10.0, 12.0}};
  bottleneck::WassersteinStats stats;
  expect_near(bottleneck::wasserstein_distance(point, separated, w1_sweep, &stats), 2.0,
              "diagonal-separated W1");
  if (stats.positive_edges != 0 || stats.pruned_pairs != 1) {
    fail("W1 exact positive-saving pruning was not recorded");
  }
  expect_near(bottleneck::wasserstein_distance(point, separated, w2_sweep), 2.0,
              "diagonal-separated W2");

  const double infinity = std::numeric_limits<double>::infinity();
  const Diagram essential_first{{1.0, infinity}, {-infinity, 2.0}};
  const Diagram essential_second{{4.0, infinity}, {-infinity, 6.0}};
  expect_near(bottleneck::wasserstein_distance(essential_first, essential_second, w1_dense), 7.0,
              "essential W1");
  expect_near(bottleneck::wasserstein_distance(essential_first, essential_second, w2_dense), 5.0,
              "essential W2");
  if (!std::isinf(bottleneck::wasserstein_distance(
          Diagram{{1.0, infinity}}, empty, w1_sweep))) {
    fail("essential count mismatch must be infinite");
  }

  // A zero-saving edge is exactly equivalent to two diagonal assignments and
  // must be safely omitted by DSR.
  const Diagram boundary_first{{-1.0, 1.0}};
  const Diagram boundary_second{{1.0, 3.0}};
  expect_near(bottleneck::wasserstein_distance(boundary_first, boundary_second, w1_dense),
              bottleneck::wasserstein_distance(boundary_first, boundary_second, w1_sweep),
              "zero-saving W1 boundary");
  expect_near(bottleneck::wasserstein_distance(boundary_first, boundary_second, w2_dense),
              bottleneck::wasserstein_distance(boundary_first, boundary_second, w2_sweep),
              "zero-saving W2 boundary");

  const Diagram disconnected_first{{0.0, 2.0}, {20.0, 22.0}};
  const Diagram disconnected_second{{0.1, 2.1}, {20.1, 22.1}};
  const WassersteinConfig component_config{
      WassersteinMetric::w1_linf,
      WassersteinCandidateStrategy::sweep_binary,
      WassersteinGraphStrategy::csr,
      WassersteinMatcherStrategy::sparse_sap,
      WassersteinComponentStrategy::tiny_sparse,
      WassersteinWarmStart::global_descending,
  };
  bottleneck::WassersteinStats component_stats;
  expect_near(bottleneck::wasserstein_distance(
                  disconnected_first, disconnected_second, component_config,
                  &component_stats),
              bottleneck::wasserstein_distance(disconnected_first,
                                                disconnected_second, w1_dense),
              "disconnected tiny components");
  if (component_stats.component_count != 2 || component_stats.tiny_components != 2) {
    fail("component decomposition did not expose two tiny components");
  }

  const WassersteinConfig certified_warm{
      WassersteinMetric::w1_linf,
      WassersteinCandidateStrategy::dense_scalar,
      WassersteinGraphStrategy::row_vectors,
      WassersteinMatcherStrategy::sparse_sap,
      WassersteinComponentStrategy::none,
      WassersteinWarmStart::global_descending,
  };
  bottleneck::WassersteinStats warm_stats;
  expect_near(bottleneck::wasserstein_distance(disconnected_first,
                                                disconnected_second,
                                                certified_warm, &warm_stats),
              bottleneck::wasserstein_distance(disconnected_first,
                                                disconnected_second, w1_dense),
              "certified greedy warm start");
  if (warm_stats.warm_start_certificates == 0) {
    fail("greedy warm start failed to certify an independently optimal matching");
  }
}

Diagram random_diagram(std::mt19937_64& generator, std::size_t size) {
  std::uniform_real_distribution<double> birth(-4.0, 4.0);
  std::uniform_real_distribution<double> persistence(0.01, 3.0);
  Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    const double value = birth(generator);
    result.push_back({value, value + persistence(generator)});
  }
  return result;
}

std::vector<WassersteinConfig> experiment_configs(WassersteinMetric metric) {
  std::vector<WassersteinConfig> configs{
      {metric, WassersteinCandidateStrategy::dense_scalar,
       WassersteinGraphStrategy::dense_matrix,
       WassersteinMatcherStrategy::dense_hungarian,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::dense_scalar,
       WassersteinGraphStrategy::dense_matrix,
       WassersteinMatcherStrategy::dense_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::dense_scalar,
       WassersteinGraphStrategy::dense_matrix,
       WassersteinMatcherStrategy::dense_sap_row_reduction,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::dense_blocked,
       WassersteinGraphStrategy::dense_matrix,
       WassersteinMatcherStrategy::dense_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::dense_avx2,
       WassersteinGraphStrategy::dense_matrix,
       WassersteinMatcherStrategy::dense_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::dense_parallel,
       WassersteinGraphStrategy::dense_matrix,
       WassersteinMatcherStrategy::dense_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::dense_matrix,
       WassersteinMatcherStrategy::dense_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::dense_scalar,
       WassersteinGraphStrategy::csr,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::csr,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::fixed_degree_4,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::fixed_degree_8,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::fixed_degree_16,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::fixed_degree_32,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::bitmask_lazy,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::bitmask_lazy,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::tiny_sparse, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::block_sparse_16,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::block_sparse_32,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_two_pointer,
       WassersteinGraphStrategy::row_vectors,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::csr,
       WassersteinMatcherStrategy::dense_sap,
       WassersteinComponentStrategy::dense, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::csr,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::sparse, WassersteinWarmStart::none},
      {metric, WassersteinCandidateStrategy::sweep_binary,
       WassersteinGraphStrategy::csr,
       WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::tiny_sparse,
       WassersteinWarmStart::row_max},
      {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 2},
      {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 4},
      {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 8},
      {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 16},
      {metric, WassersteinCandidateStrategy::topk_pricing_full_scan,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 32},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 2},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 4},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 8},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 16},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 32},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 2},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 4},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 8},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 16},
      {metric, WassersteinCandidateStrategy::topk_pricing_sweep_incremental,
       WassersteinGraphStrategy::csr, WassersteinMatcherStrategy::sparse_sap,
       WassersteinComponentStrategy::none, WassersteinWarmStart::none,
       WassersteinDuplicateStrategy::none, 32},
      {metric, WassersteinCandidateStrategy::adaptive,
       WassersteinGraphStrategy::adaptive,
       WassersteinMatcherStrategy::adaptive,
       WassersteinComponentStrategy::adaptive,
       WassersteinWarmStart::global_descending},
      {metric, WassersteinCandidateStrategy::adaptive,
       WassersteinGraphStrategy::adaptive,
       WassersteinMatcherStrategy::adaptive,
       WassersteinComponentStrategy::adaptive,
       WassersteinWarmStart::global_descending,
       WassersteinDuplicateStrategy::exact},
      {metric, WassersteinCandidateStrategy::adaptive,
       WassersteinGraphStrategy::adaptive,
       WassersteinMatcherStrategy::adaptive,
       WassersteinComponentStrategy::adaptive,
       WassersteinWarmStart::global_descending,
       WassersteinDuplicateStrategy::adaptive},
  };
  for (std::size_t index = 0; index + 3 < configs.size(); ++index) {
    configs[index].duplicates = WassersteinDuplicateStrategy::none;
  }
  return configs;
}

void brute_force_differential() {
  std::mt19937_64 generator(0xD5A0A11ULL);
  std::uniform_int_distribution<int> size(0, 5);
  for (WassersteinMetric metric :
       {WassersteinMetric::w1_linf, WassersteinMetric::w2_l2}) {
    const auto configs = experiment_configs(metric);
    for (int trial = 0; trial < 150; ++trial) {
      const Diagram first = random_diagram(generator, static_cast<std::size_t>(size(generator)));
      const Diagram second = random_diagram(generator, static_cast<std::size_t>(size(generator)));
      const double expected = brute_force_distance(first, second, metric);
      for (const WassersteinConfig& config : configs) {
        expect_near(bottleneck::wasserstein_distance(first, second, config), expected,
                    "variant vs brute-force trial " + std::to_string(trial));
      }
    }
  }
}

void variant_differential() {
  std::mt19937_64 generator(0x5A33BULL);
  std::uniform_int_distribution<int> size(0, 24);
  constexpr int trial_count = 500;
  for (WassersteinMetric metric :
       {WassersteinMetric::w1_linf, WassersteinMetric::w2_l2}) {
    const auto configs = experiment_configs(metric);
    for (int trial = 0; trial < trial_count; ++trial) {
      Diagram first = random_diagram(generator, static_cast<std::size_t>(size(generator)));
      Diagram second = random_diagram(generator, static_cast<std::size_t>(size(generator)));
      if (trial % 13 == 0 && !first.empty()) {
        first.push_back(first.front());
      }
      if (trial % 17 == 0 && !second.empty()) {
        second.push_back(second.front());
      }
      const PreparedDiagram prepared_first(first);
      const PreparedDiagram prepared_second(second);
      const double expected =
          bottleneck::wasserstein_distance(prepared_first, prepared_second,
                                            configs.front());
      for (std::size_t config_index = 0; config_index < configs.size(); ++config_index) {
        const WassersteinConfig& config = configs[config_index];
        bottleneck::WassersteinStats stats;
        const double actual = bottleneck::wasserstein_distance(
            prepared_first, prepared_second, config, &stats);
        const double diagnostic_scale =
            (std::max)({1.0, std::fabs(actual), std::fabs(expected)});
        if (std::fabs(actual - expected) > 2e-12 * diagnostic_scale &&
            (config.candidates ==
                 WassersteinCandidateStrategy::topk_pricing_full_scan ||
             config.candidates ==
                 WassersteinCandidateStrategy::topk_pricing_sweep ||
             config.candidates == WassersteinCandidateStrategy::
                                      topk_pricing_sweep_incremental)) {
          std::cerr << "top-k diagnostic: rounds=" << stats.pricing_rounds
                    << ", priced=" << stats.priced_edges
                    << ", violations=" << stats.pricing_violations
                    << ", materialized=" << stats.peak_materialized_edges
                    << ", positive=" << stats.positive_edges << '\n';
        }
        expect_near(actual,
                    expected, "variant differential config " +
                                  std::to_string(config_index) + " trial " +
                                  std::to_string(trial));
      }
      expect_near(bottleneck::wasserstein_distance(prepared_second, prepared_first,
                                                    configs.back()),
                  expected, "adaptive symmetry trial " + std::to_string(trial));
    }
  }
}

void rectangular_dense_sap_differential() {
  std::mt19937_64 generator(0x5A9A11ULL);
  constexpr std::array<std::pair<std::size_t, std::size_t>, 6> shapes{{
      {1, 64}, {64, 1}, {7, 64}, {64, 7}, {31, 128}, {128, 31},
  }};
  for (WassersteinMetric metric :
       {WassersteinMetric::w1_linf, WassersteinMetric::w2_l2}) {
    const WassersteinConfig baseline{
        metric, WassersteinCandidateStrategy::dense_scalar,
        WassersteinGraphStrategy::dense_matrix,
        WassersteinMatcherStrategy::dense_hungarian,
        WassersteinComponentStrategy::none, WassersteinWarmStart::none};
    WassersteinConfig rectangular = baseline;
    rectangular.matcher = WassersteinMatcherStrategy::dense_sap;
    for (const auto [rows, columns] : shapes) {
      for (int trial = 0; trial < 12; ++trial) {
        Diagram first = random_diagram(generator, rows);
        Diagram second = random_diagram(generator, columns);
        if (trial % 4 == 0) {
          first.push_back(first.front());
          second.push_back(second.front());
        }
        const double expected = bottleneck::wasserstein_distance(first, second, baseline);
        expect_near(bottleneck::wasserstein_distance(first, second, rectangular),
                    expected, "rectangular dense SAP " + std::to_string(rows) + "x" +
                                  std::to_string(columns),
                    2e-11);
      }
    }
  }
}

void near_identical_dense_w2_regression() {
  Diagram first;
  Diagram second;
  first.reserve(128);
  second.reserve(128);
  long double expected_power = 0;
  for (std::size_t index = 0; index < 128; ++index) {
    const double birth = -0.01 + static_cast<double>(index) * 1e-4;
    const double death = birth + 2.0 + static_cast<double>(index % 7) * 1e-5;
    const double shift = static_cast<double>(1 + index % 3) * 1e-9;
    first.push_back({birth, death});
    second.push_back({birth + shift, death + shift});
    expected_power += cross_power(first.back(), second.back(),
                                  WassersteinMetric::w2_l2);
  }
  const double expected = static_cast<double>(std::sqrt(expected_power));
  const auto configs = experiment_configs(WassersteinMetric::w2_l2);
  for (std::size_t config_index = 0; config_index < configs.size(); ++config_index) {
    expect_near(bottleneck::wasserstein_distance(first, second,
                                                  configs[config_index]),
                expected,
                "near-identical dense W2 config " +
                    std::to_string(config_index),
                2e-14);
  }
}

void edge_case_differential() {
  const std::vector<std::pair<Diagram, Diagram>> cases{
      {Diagram{{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}},
       Diagram{{0.0, 1.0}, {0.0, 1.0}}},
      {Diagram{{0.0, 1.0}, {2.0, 3.0}},
       Diagram{{0.5, 1.5}, {0.5, 1.5}, {2.5, 3.5}, {8.0, 9.0}}},
      {Diagram{{-4.0, -3.999999999}, {0.0, 1e-9}, {4.0, 4.000000001}},
       Diagram{{-4.0, -3.999999998}, {4.0, 4.000000002}}},
      {Diagram{{1e100, 1e100 + 1e90}, {-1e100, -1e100 + 2e90}},
       Diagram{{1e100 + 2e89, 1e100 + 1.2e90}}},
      {Diagram{{-2.0, 0.0}},
       Diagram{{-2.0, 0.0}, {-1.0, 1.0}, {0.0, 2.0}, {1.0, 3.0}, {2.0, 4.0}}},
  };
  for (WassersteinMetric metric :
       {WassersteinMetric::w1_linf, WassersteinMetric::w2_l2}) {
    const auto configs = experiment_configs(metric);
    for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
      const double expected =
          brute_force_distance(cases[case_index].first, cases[case_index].second, metric);
      for (const WassersteinConfig& config : configs) {
        expect_near(bottleneck::wasserstein_distance(cases[case_index].first,
                                                      cases[case_index].second,
                                                      config),
                    expected,
                    "edge-case differential " + std::to_string(case_index), 2e-11);
      }
    }
  }

  const double infinity = std::numeric_limits<double>::infinity();
  const Diagram first{{0.0, 2.0}, {1.0, infinity}, {-infinity, 3.0}};
  const Diagram second{{0.2, 2.2}, {2.0, infinity}, {-infinity, 5.0}};
  for (WassersteinMetric metric :
       {WassersteinMetric::w1_linf, WassersteinMetric::w2_l2}) {
    const auto configs = experiment_configs(metric);
    const double expected = bottleneck::wasserstein_distance(first, second, configs.front());
    for (const WassersteinConfig& config : configs) {
      expect_near(bottleneck::wasserstein_distance(first, second, config), expected,
                  "essential plus finite differential");
    }
  }
}

void batch_cases() {
  const PreparedDiagram query(Diagram{{0.0, 2.0}, {4.0, 5.0}});
  const std::vector<PreparedDiagram> targets{
      PreparedDiagram(Diagram{}),
      PreparedDiagram(Diagram{{0.1, 2.1}}),
      PreparedDiagram(Diagram{{0.2, 2.2}, {4.1, 5.1}}),
  };
  const auto allocated = bottleneck::wasserstein_distances(query, targets);
  std::vector<double> output(targets.size(), -1.0);
  bottleneck::wasserstein_distances(query, targets, output);
  bottleneck::WassersteinWorkspace workspace;
  std::vector<double> workspace_output(targets.size(), -1.0);
  bottleneck::wasserstein_distances(query, targets, workspace_output, workspace);
  for (std::size_t index = 0; index < targets.size(); ++index) {
    const double expected = bottleneck::wasserstein_distance(query, targets[index]);
    expect_near(allocated[index], expected, "allocated Wasserstein batch");
    expect_near(output[index], expected, "caller-buffer Wasserstein batch");
    expect_near(workspace_output[index], expected,
                "workspace Wasserstein batch");
    expect_near(bottleneck::wasserstein_distance(query, targets[index], workspace),
                expected, "workspace Wasserstein distance");
  }
  bool rejected = false;
  try {
    std::vector<double> wrong_size(targets.size() - 1);
    bottleneck::wasserstein_distances(query, targets, wrong_size);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  if (!rejected) {
    fail("Wasserstein batch must reject a mismatched output span");
  }
}

}  // namespace

int main() {
  deterministic_cases();
  std::cout << "deterministic passed\n" << std::flush;
  brute_force_differential();
  std::cout << "brute-force differential passed\n" << std::flush;
  variant_differential();
  std::cout << "variant differential passed\n" << std::flush;
  rectangular_dense_sap_differential();
  std::cout << "rectangular dense SAP differential passed\n" << std::flush;
  near_identical_dense_w2_regression();
  std::cout << "near-identical dense W2 regression passed\n" << std::flush;
  edge_case_differential();
  batch_cases();
  std::cout << "wasserstein_core_tests: all checks passed\n";
  return 0;
}
