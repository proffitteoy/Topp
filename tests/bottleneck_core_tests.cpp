#include <bottleneck/core.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using bottleneck::AdjacencyStrategy;
using bottleneck::CandidateStrategy;
using bottleneck::Diagram;
using bottleneck::DistanceStrategy;
using bottleneck::MatcherStrategy;
using bottleneck::Point;
using bottleneck::PreparedDiagram;
using bottleneck::SolverConfig;
using bottleneck::SolverStats;
using bottleneck::ThresholdStrategy;
using bottleneck::VertexOrder;

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "FAILED: " << message << '\n';
  std::exit(1);
}

void expect_equal(double actual, double expected, const std::string& label) {
  if (actual != expected && !(std::isnan(actual) && std::isnan(expected))) {
    fail(label + ": expected " + std::to_string(expected) + ", got " + std::to_string(actual));
  }
}

std::vector<SolverConfig> all_configs() {
  std::vector<SolverConfig> result;
  for (CandidateStrategy candidates : {CandidateStrategy::sort_all, CandidateStrategy::sort_unique,
                                       CandidateStrategy::sort_unique_clipped,
                                       CandidateStrategy::sort_unique_greedy_clipped,
                                       CandidateStrategy::x_sweep_clipped}) {
    for (ThresholdStrategy threshold : {ThresholdStrategy::binary, ThresholdStrategy::gudhi_alpha,
                                        ThresholdStrategy::exponential,
                                        ThresholdStrategy::quickselect,
                                        ThresholdStrategy::incremental,
                                        ThresholdStrategy::incremental_blocked}) {
      for (DistanceStrategy distance : {DistanceStrategy::dense_aos, DistanceStrategy::dense_soa,
                                        DistanceStrategy::dense_soa_avx2,
                                        DistanceStrategy::recompute_aos,
                                        DistanceStrategy::recompute_soa}) {
        for (AdjacencyStrategy adjacency : {AdjacencyStrategy::on_demand,
                                            AdjacencyStrategy::dense_byte,
                                            AdjacencyStrategy::sparse_csr,
                                            AdjacencyStrategy::x_sweep_csr,
                                            AdjacencyStrategy::bitset64,
                                            AdjacencyStrategy::bitset128,
                                            AdjacencyStrategy::bitset256,
                                            AdjacencyStrategy::bitset_auto,
                                            AdjacencyStrategy::adaptive,
                                            AdjacencyStrategy::dynamic_bitset}) {
          for (MatcherStrategy matcher : {MatcherStrategy::kuhn, MatcherStrategy::greedy_kuhn,
                                          MatcherStrategy::reusable_greedy_kuhn,
                                          MatcherStrategy::fixed_greedy_kuhn,
                                          MatcherStrategy::constraint_kuhn,
                                          MatcherStrategy::component_kuhn,
                                          MatcherStrategy::mandatory_flow,
                                          MatcherStrategy::multiplicity_flow,
                                          MatcherStrategy::adaptive,
                                           MatcherStrategy::hopcroft_karp,
                                           MatcherStrategy::greedy_hopcroft_karp,
                                           MatcherStrategy::mandatory_sparse_flow}) {
            for (VertexOrder order : {VertexOrder::natural, VertexOrder::degree_ascending}) {
              result.push_back(SolverConfig{candidates, threshold, distance, adjacency, matcher, order});
            }
          }
        }
      }
    }
  }
  return result;
}

void deterministic_cases() {
  const double infinity = std::numeric_limits<double>::infinity();
  const SolverConfig reference{CandidateStrategy::sort_all, ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn,
                               VertexOrder::natural};
  expect_equal(bottleneck::bottleneck_distance(Diagram{}, Diagram{}, reference), 0.0, "empty");
  expect_equal(bottleneck::bottleneck_distance(Diagram{}, Diagram{{0.0, 1.0}}, reference), 0.5,
               "empty-to-point");
  expect_equal(bottleneck::bottleneck_distance(Diagram{{0.0, 1.0}}, Diagram{{0.0, 1.0}}, reference),
               0.0, "identity");
  expect_equal(bottleneck::bottleneck_distance(Diagram{{0.0, 1.0}}, Diagram{{0.25, 1.25}}, reference),
               0.25, "translated-point");
  expect_equal(bottleneck::bottleneck_distance(Diagram{{0.0, 0.0}}, Diagram{}, reference), 0.0,
               "diagonal-point-ignored");
  expect_equal(bottleneck::bottleneck_distance(Diagram{{1.0, infinity}},
                                                Diagram{{1.75, infinity}}, reference),
               0.75, "positive-essential");
  expect_equal(bottleneck::bottleneck_distance(Diagram{{-infinity, 1.0}},
                                                Diagram{{-infinity, 2.0}}, reference),
               1.0, "negative-essential");
  if (!std::isinf(bottleneck::bottleneck_distance(Diagram{{1.0, infinity}}, Diagram{}, reference))) {
    fail("essential-count-mismatch must be infinite");
  }
}

Diagram random_diagram(std::mt19937_64& generator, std::size_t size) {
  std::uniform_real_distribution<double> birth(-2.0, 2.0);
  std::uniform_real_distribution<double> persistence(0.0, 2.0);
  Diagram result;
  result.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    const double b = birth(generator);
    result.push_back(Point{b, b + persistence(generator)});
  }
  return result;
}

void randomized_differential() {
  std::mt19937_64 generator(0xB0771EULL);
  const auto configs = all_configs();
  const SolverConfig reference{CandidateStrategy::sort_all, ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn,
                               VertexOrder::natural};
  std::uniform_int_distribution<int> size_distribution(0, 12);

  for (int trial = 0; trial < 250; ++trial) {
    const Diagram first = random_diagram(generator, static_cast<std::size_t>(size_distribution(generator)));
    const Diagram second = random_diagram(generator, static_cast<std::size_t>(size_distribution(generator)));
    const double expected = bottleneck::bottleneck_distance(first, second, reference);
    const double reverse = bottleneck::bottleneck_distance(second, first, reference);
    expect_equal(reverse, expected, "symmetry trial " + std::to_string(trial));

    for (const SolverConfig& config : configs) {
      const double actual = bottleneck::bottleneck_distance(first, second, config);
      expect_equal(actual, expected, "variant differential trial " + std::to_string(trial));
      const PreparedDiagram prepared_first(first);
      const PreparedDiagram prepared_second(second);
      if (!bottleneck::bottleneck_within(prepared_first, prepared_second, expected, config)) {
        fail("within rejected exact radius at trial " + std::to_string(trial));
      }
      if (expected > 0.0) {
        const double previous = std::nextafter(expected, 0.0);
        if (bottleneck::bottleneck_within(prepared_first, prepared_second, previous, config)) {
          fail("within accepted radius below exact result at trial " + std::to_string(trial));
        }
      }
    }
  }
}

void geometric_differential() {
  std::mt19937_64 generator(0x6E0B0771EULL);
  const SolverConfig reference{CandidateStrategy::sort_all, ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn, VertexOrder::natural};
  const std::vector<SolverConfig> configs{
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::binary,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::geometric_hopcroft_karp, VertexOrder::natural},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::geometric_hopcroft_karp, VertexOrder::natural},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::incremental_blocked,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::geometric_hopcroft_karp, VertexOrder::natural},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::geometric_refinement,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::geometric_hopcroft_karp, VertexOrder::natural},
      {},
  };
  std::uniform_int_distribution<int> size_distribution(0, 40);

  for (int trial = 0; trial < 300; ++trial) {
    Diagram first = random_diagram(
        generator, static_cast<std::size_t>(size_distribution(generator)));
    Diagram second = random_diagram(
        generator, static_cast<std::size_t>(size_distribution(generator)));
    if (trial % 7 == 0 && !first.empty()) {
      first.push_back(first.front());
    }
    if (trial % 11 == 0 && !second.empty()) {
      second.push_back(second.front());
    }
    const double expected = bottleneck::bottleneck_distance(first, second, reference);
    const PreparedDiagram prepared_first(first);
    const PreparedDiagram prepared_second(second);
    for (const SolverConfig& config : configs) {
      const double actual = bottleneck::bottleneck_distance(
          prepared_first, prepared_second, config);
      expect_equal(actual, expected, "geometric differential trial " + std::to_string(trial));
      if (!bottleneck::bottleneck_within(
              prepared_first, prepared_second, expected, config)) {
        fail("geometric within rejected exact radius at trial " + std::to_string(trial));
      }
      if (expected > 0.0 && bottleneck::bottleneck_within(
                                prepared_first, prepared_second,
                                std::nextafter(expected, 0.0), config)) {
        fail("geometric within accepted radius below exact result at trial " +
             std::to_string(trial) + ", expected=" + std::to_string(expected) +
             ", n=" + std::to_string(first.size()) + ", m=" +
             std::to_string(second.size()) + ", threshold=" +
             bottleneck::to_string(config.threshold));
      }
    }
  }
}

void geometric_floating_boundaries() {
  const SolverConfig reference{CandidateStrategy::sort_all, ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn, VertexOrder::natural};
  const SolverConfig geometric{CandidateStrategy::sort_unique_clipped,
                               ThresholdStrategy::geometric_refinement,
                               DistanceStrategy::recompute_soa,
                               AdjacencyStrategy::on_demand,
                               MatcherStrategy::geometric_hopcroft_karp,
                               VertexOrder::natural};
  const double base = 1.0e12;
  const Diagram first{{base, base + 1.0}, {-base, -base + 0.5}};
  const Diagram second{{std::nextafter(base, std::numeric_limits<double>::infinity()),
                        std::nextafter(base + 1.0, std::numeric_limits<double>::infinity())},
                       {-base, -base + 0.5}};
  const double expected = bottleneck::bottleneck_distance(first, second, reference);
  expect_equal(bottleneck::bottleneck_distance(first, second, geometric), expected,
               "geometric large-coordinate ulp boundary");
  const PreparedDiagram prepared_first(first);
  const PreparedDiagram prepared_second(second);
  if (!bottleneck::bottleneck_within(prepared_first, prepared_second, expected, geometric)) {
    fail("geometric large-coordinate exact threshold rejected");
  }
  if (bottleneck::bottleneck_within(prepared_first, prepared_second,
                                    std::nextafter(expected, 0.0), geometric)) {
    fail("geometric large-coordinate previous threshold accepted");
  }
}

void multiplicity_flow_differential() {
  std::mt19937_64 generator(0xD0B1CA7EULL);
  const SolverConfig reference{CandidateStrategy::sort_all, ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn, VertexOrder::natural};
  const std::vector<SolverConfig> compressed{
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::binary,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::multiplicity_flow, VertexOrder::natural},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::multiplicity_flow, VertexOrder::natural},
  };
  std::uniform_int_distribution<int> unique_count(0, 7);
  std::uniform_int_distribution<int> multiplicity(1, 8);
  std::uniform_real_distribution<double> birth(-2.0, 2.0);
  std::uniform_real_distribution<double> persistence(0.0, 2.0);

  for (int trial = 0; trial < 300; ++trial) {
    const auto make_duplicate_heavy = [&]() {
      Diagram diagram;
      const int groups = unique_count(generator);
      for (int group = 0; group < groups; ++group) {
        const double point_birth = birth(generator);
        const Point point{point_birth, point_birth + persistence(generator)};
        diagram.insert(diagram.end(), static_cast<std::size_t>(multiplicity(generator)), point);
      }
      std::shuffle(diagram.begin(), diagram.end(), generator);
      return diagram;
    };
    const Diagram first = make_duplicate_heavy();
    const Diagram second = make_duplicate_heavy();
    const double expected = bottleneck::bottleneck_distance(first, second, reference);
    const PreparedDiagram prepared_first(first);
    const PreparedDiagram prepared_second(second);

    for (const SolverConfig& config : compressed) {
      SolverStats stats;
      const double actual = bottleneck::bottleneck_distance(
          prepared_first, prepared_second, config, &stats);
      expect_equal(actual, expected,
                   "multiplicity differential trial " + std::to_string(trial));
      if (!bottleneck::bottleneck_within(
              prepared_first, prepared_second, expected, config)) {
        fail("multiplicity within rejected exact radius at trial " +
             std::to_string(trial));
      }
      if (expected > 0.0 && bottleneck::bottleneck_within(
                                prepared_first, prepared_second,
                                std::nextafter(expected, 0.0), config)) {
        fail("multiplicity within accepted radius below exact result at trial " +
             std::to_string(trial));
      }
    }
  }

  const Diagram repeated_first(64, Point{0.0, 4.0});
  const Diagram repeated_second(64, Point{0.5, 4.5});
  const PreparedDiagram prepared_first(repeated_first);
  const PreparedDiagram prepared_second(repeated_second);
  SolverStats stats;
  if (!bottleneck::bottleneck_within(prepared_first, prepared_second, 0.5,
                                     compressed.front(), &stats)) {
    fail("multiplicity decision rejected a feasible repeated threshold");
  }
  if (stats.multiplicity_groups != 2 || stats.multiplicity_points_removed != 126 ||
      stats.capacity_edges != 1) {
    fail("multiplicity decision did not use the compressed capacity graph");
  }
}

void mandatory_partial_flow_differential() {
  std::mt19937_64 generator(0xE7A7DA7AULL);
  const SolverConfig reference{CandidateStrategy::sort_all, ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn, VertexOrder::natural};
  const std::vector<SolverConfig> variants{
      {CandidateStrategy::x_sweep_clipped, ThresholdStrategy::binary,
       DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
       MatcherStrategy::kuhn, VertexOrder::natural},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::mandatory_flow, VertexOrder::natural},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::mandatory_sparse_flow, VertexOrder::natural},
      {CandidateStrategy::x_sweep_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::mandatory_sparse_flow, VertexOrder::natural},
  };
  const std::vector<std::size_t> ratios{1, 2, 4, 8, 16};
  std::uniform_real_distribution<double> birth(-4.0, 4.0);
  std::uniform_real_distribution<double> near_persistence(0.001, 0.25);
  std::uniform_real_distribution<double> wide_persistence(0.5, 2.0);

  for (int trial = 0; trial < 150; ++trial) {
    const std::size_t first_size = 1 + static_cast<std::size_t>(trial % 8);
    const std::size_t second_size = first_size * ratios[static_cast<std::size_t>(trial) % ratios.size()];
    const auto make_diagram = [&](std::size_t size, int salt) {
      Diagram diagram;
      diagram.reserve(size);
      for (std::size_t index = 0; index < size; ++index) {
        const double point_birth = birth(generator);
        const bool wide = (index + static_cast<std::size_t>(salt)) % 5 == 0;
        const double persistence = wide ? wide_persistence(generator)
                                        : near_persistence(generator);
        diagram.push_back({point_birth, point_birth + persistence});
      }
      return diagram;
    };
    const Diagram first = make_diagram(first_size, trial);
    const Diagram second = make_diagram(second_size, trial + 2);
    const double expected = bottleneck::bottleneck_distance(first, second, reference);
    const PreparedDiagram prepared_first(first);
    const PreparedDiagram prepared_second(second);
    for (const SolverConfig& config : variants) {
      expect_equal(bottleneck::bottleneck_distance(prepared_first, prepared_second, config),
                   expected, "mandatory partial differential trial " +
                                 std::to_string(trial));
      expect_equal(bottleneck::bottleneck_distance(prepared_second, prepared_first, config),
                   expected, "mandatory partial symmetry trial " +
                                 std::to_string(trial));
      if (!bottleneck::bottleneck_within(prepared_first, prepared_second, expected, config)) {
        fail("mandatory partial flow rejected exact radius at trial " +
             std::to_string(trial));
      }
      if (expected > 0.0 && bottleneck::bottleneck_within(
                                prepared_first, prepared_second,
                                std::nextafter(expected, 0.0), config)) {
        fail("mandatory partial flow accepted previous ULP at trial " +
             std::to_string(trial));
      }
    }
  }

  const SolverConfig sparse{CandidateStrategy::x_sweep_clipped,
                            ThresholdStrategy::quickselect,
                            DistanceStrategy::recompute_soa,
                            AdjacencyStrategy::on_demand,
                            MatcherStrategy::mandatory_sparse_flow,
                            VertexOrder::natural};
  const PreparedDiagram optional_first(Diagram{{0.0, 0.1}, {2.0, 2.1}, {4.0, 4.1}});
  const PreparedDiagram optional_second(
      Diagram{{0.01, 0.11}, {1.0, 1.1}, {2.01, 2.11}, {3.0, 3.1}, {4.01, 4.11}});
  SolverStats optional_stats;
  if (!bottleneck::bottleneck_within(optional_first, optional_second, 0.051,
                                     sparse, &optional_stats)) {
    fail("all-optional graph must be feasible without cross edges");
  }
  if (optional_stats.mandatory_vertices != 0 ||
      optional_stats.optional_pairs_pruned != 15 ||
      optional_stats.adjacency_checks != 0 ||
      optional_stats.x_window_candidates != 0) {
    fail("all-optional graph did not bypass optional cross pairs");
  }

  const PreparedDiagram boundary_first(Diagram{{0.0, 0.41}, {10.0, 10.1}});
  const PreparedDiagram boundary_second(Diagram{{0.2, 0.4}, {20.0, 20.1}});
  SolverStats boundary_stats;
  if (!bottleneck::bottleneck_within(boundary_first, boundary_second, 0.2,
                                     sparse, &boundary_stats)) {
    fail("mandatory x-window rejected an exact birth-boundary edge");
  }
  if (boundary_stats.mandatory_vertices != 1 ||
      boundary_stats.optional_pairs_pruned != 2 ||
      boundary_stats.x_window_candidates != 1 ||
      boundary_stats.adjacency_checks != 1 ||
      boundary_stats.capacity_edges != 1) {
    fail("mandatory x-window boundary stats are inconsistent");
  }
  if (bottleneck::bottleneck_within(boundary_first, boundary_second,
                                    std::nextafter(0.2, 0.0), sparse)) {
    fail("mandatory x-window accepted a threshold below the boundary edge");
  }

  const PreparedDiagram death_outside(Diagram{{0.0, 0.41}});
  const PreparedDiagram same_birth_optional(Diagram{{0.1, 0.2}});
  SolverStats death_stats;
  if (bottleneck::bottleneck_within(death_outside, same_birth_optional, 0.2,
                                    sparse, &death_stats)) {
    fail("x-window candidate must still pass the exact death-coordinate check");
  }
  if (death_stats.x_window_candidates != 1 || death_stats.capacity_edges != 0) {
    fail("x-window exact-filter stats are inconsistent");
  }

  Diagram spaced_first;
  Diagram spaced_second;
  for (int index = 0; index < 16; ++index) {
    spaced_first.push_back({static_cast<double>(index) * 10.0,
                            static_cast<double>(index) * 10.0 + 0.1});
    spaced_second.push_back({static_cast<double>(index) * 10.0 + 0.01,
                             static_cast<double>(index) * 10.0 + 0.11});
  }
  const SolverConfig candidate_only{CandidateStrategy::x_sweep_clipped,
                                    ThresholdStrategy::binary,
                                    DistanceStrategy::dense_aos,
                                    AdjacencyStrategy::on_demand,
                                    MatcherStrategy::kuhn,
                                    VertexOrder::natural};
  const double spaced_expected = bottleneck::bottleneck_distance(
      PreparedDiagram(spaced_first), PreparedDiagram(spaced_second), reference);
  SolverStats candidate_stats;
  expect_equal(bottleneck::bottleneck_distance(PreparedDiagram(spaced_first),
                                               PreparedDiagram(spaced_second), candidate_only,
                                               &candidate_stats),
               spaced_expected, "x-window clipped candidate distance");
  const std::uint64_t full_candidate_count = 1 + 16 + 16 + 16 * 16;
  if (candidate_stats.raw_candidates != full_candidate_count ||
      candidate_stats.x_window_candidates >= 16 * 16 ||
      candidate_stats.clipped_candidates == 0) {
    fail("x-window candidate generation did not prune the full cross product");
  }
}

void component_structure_stats() {
  const double shift = std::ldexp(1.0, -10);
  const PreparedDiagram first(Diagram{{0.0, 2.0}, {8.0, 10.0}});
  const PreparedDiagram second(
      Diagram{{shift, 2.0 + shift}, {8.0 + shift, 10.0 + shift}});
  const SolverConfig component{CandidateStrategy::x_sweep_clipped,
                               ThresholdStrategy::quickselect,
                               DistanceStrategy::recompute_soa,
                               AdjacencyStrategy::x_sweep_csr,
                               MatcherStrategy::component_kuhn,
                               VertexOrder::degree_ascending};
  SolverStats stats;
  if (!bottleneck::bottleneck_within(first, second, shift, component, &stats)) {
    fail("component structure threshold must be feasible");
  }
  if (stats.component_count != 3 || stats.largest_component_vertices != 4) {
    fail("component structure stats did not expose the augmented dummy component");
  }
}

void router_x_window_stats() {
  Diagram first;
  Diagram second;
  first.reserve(64);
  second.reserve(64);
  for (int index = 0; index < 64; ++index) {
    const double first_birth = static_cast<double>(index) * 0.125;
    const double second_birth = static_cast<double>(index) * 0.125 + 0.03125;
    first.push_back({first_birth, first_birth + 0.5 +
                                  static_cast<double>(index % 5) * 0.0625});
    second.push_back({second_birth, second_birth + 0.375 +
                                    static_cast<double>(index % 7) * 0.0625});
  }
  const PreparedDiagram prepared_first(first);
  const PreparedDiagram prepared_second(second);
  const double upper = (std::max)(prepared_first.max_finite_diagonal_distance(),
                                  prepared_second.max_finite_diagonal_distance());
  std::uint64_t expected_pairs = 0;
  for (double left : prepared_first.finite_births()) {
    for (double right : prepared_second.finite_births()) {
      expected_pairs += static_cast<std::uint64_t>(std::fabs(left - right) <= upper);
    }
  }
  SolverStats stats;
  (void)bottleneck::bottleneck_distance(prepared_first, prepared_second, SolverConfig{}, &stats);
  if (stats.router_x_window_pairs != expected_pairs) {
    fail("linear router x-window count disagrees with brute force");
  }
}

void adaptive_dispatcher_differential() {
  std::mt19937_64 generator(0xADAF71EULL);
  const SolverConfig reference{CandidateStrategy::sort_all, ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn, VertexOrder::natural};
  const SolverConfig adaptive{};
  for (int trial = 0; trial < 40; ++trial) {
    const Diagram first = random_diagram(generator, 64);
    const Diagram second = random_diagram(generator, 64);
    expect_equal(bottleneck::bottleneck_distance(first, second, adaptive),
                 bottleneck::bottleneck_distance(first, second, reference),
                 "adaptive symmetric differential trial " + std::to_string(trial));
  }
  for (int trial = 0; trial < 20; ++trial) {
    const Diagram first = random_diagram(generator, 32);
    const Diagram second = random_diagram(generator, 256);
    expect_equal(bottleneck::bottleneck_distance(first, second, adaptive),
                 bottleneck::bottleneck_distance(first, second, reference),
                 "adaptive asymmetric differential trial " + std::to_string(trial));
  }
  const Diagram separated_first{{-10.0, -9.0}, {-8.0, -5.0}};
  Diagram separated_second;
  separated_second.reserve(128);
  for (int index = 0; index < 128; ++index) {
    const double birth = 10.0 + static_cast<double>(index) / 128.0;
    separated_second.push_back({birth, birth + 0.25});
  }
  SolverStats separated_stats;
  expect_equal(bottleneck::bottleneck_distance(PreparedDiagram(separated_first),
                                               PreparedDiagram(separated_second), adaptive,
                                               &separated_stats),
               bottleneck::bottleneck_distance(separated_first, separated_second, reference),
               "adaptive no-cross exact shortcut");
  if (separated_stats.router_no_cross_shortcuts != 1) {
    fail("adaptive dispatcher did not record the no-cross shortcut");
  }

  SolverStats small_stats;
  (void)bottleneck::bottleneck_distance(PreparedDiagram(random_diagram(generator, 32)),
                                        PreparedDiagram(random_diagram(generator, 32)),
                                        adaptive, &small_stats);
  if (small_stats.router_quickselect_routes != 1) {
    fail("adaptive dispatcher did not record the small-N quickselect route");
  }

  const Diagram repeated_first(512, Point{0.0, 4.0});
  const Diagram repeated_second(512, Point{0.25, 4.25});
  SolverStats stats;
  expect_equal(bottleneck::bottleneck_distance(
                   PreparedDiagram(repeated_first), PreparedDiagram(repeated_second),
                   adaptive, &stats),
               bottleneck::bottleneck_distance(repeated_first, repeated_second, reference),
               "adaptive multiplicity route");
  if (stats.multiplicity_points_removed == 0) {
    fail("adaptive dispatcher did not select multiplicity flow");
  }
  if (stats.router_multiplicity_routes != 1) {
    fail("adaptive dispatcher did not record the multiplicity route");
  }
  SolverStats within_stats;
  if (!bottleneck::bottleneck_within(
          PreparedDiagram(repeated_first), PreparedDiagram(repeated_second), 0.25,
          adaptive, &within_stats) || within_stats.multiplicity_points_removed == 0) {
    fail("adaptive within did not select multiplicity flow");
  }

  SolverStats identity_stats;
  expect_equal(bottleneck::bottleneck_distance(
                   PreparedDiagram(repeated_first), PreparedDiagram(repeated_first),
                   adaptive, &identity_stats),
               0.0, "adaptive repeated identity certificate");
  if (identity_stats.multiplicity_groups != 0 ||
      identity_stats.geometric_queries != 0 ||
      identity_stats.router_identity_shortcuts != 1) {
    fail("identity certificate should bypass matching kernels");
  }

  Diagram partial_first;
  Diagram partial_second;
  partial_first.reserve(256);
  partial_second.reserve(256);
  const double exact_shift = std::ldexp(1.0, -10);
  for (int index = 0; index < 256; ++index) {
    const double birth = static_cast<double>(index) * 0.125;
    const double persistence = index < 230 ? std::ldexp(1.0, -11) : 1.0;
    partial_first.push_back({birth, birth + persistence});
    partial_second.push_back({birth + exact_shift, birth + persistence + exact_shift});
  }
  SolverStats partial_stats;
  if (!bottleneck::bottleneck_within(PreparedDiagram(partial_first),
                                     PreparedDiagram(partial_second), exact_shift,
                                     adaptive, &partial_stats) ||
      partial_stats.mandatory_vertices == 0 ||
      partial_stats.optional_pairs_pruned == 0 ||
      partial_stats.x_window_candidates == 0) {
    fail("adaptive within did not select mandatory sparse flow");
  }
  SolverStats partial_distance_stats;
  expect_equal(bottleneck::bottleneck_distance(PreparedDiagram(partial_first),
                                               PreparedDiagram(partial_second), adaptive,
                                               &partial_distance_stats),
               exact_shift, "adaptive mandatory sparse distance route");
  if (partial_distance_stats.mandatory_vertices == 0 ||
      partial_distance_stats.optional_pairs_pruned == 0 ||
      partial_distance_stats.router_mandatory_routes != 1) {
    fail("adaptive distance did not select mandatory sparse flow");
  }

  Diagram medium_first;
  Diagram medium_second;
  medium_first.reserve(128);
  medium_second.reserve(128);
  for (std::size_t index = 0; index < 128; ++index) {
    const double birth = static_cast<double>(index) * 0.125;
    const double persistence = index < 115 ? std::ldexp(1.0, -11) : 1.0;
    medium_first.push_back({birth, birth + persistence});
    medium_second.push_back({birth + exact_shift, birth + persistence + exact_shift});
  }
  SolverStats medium_stats;
  expect_equal(bottleneck::bottleneck_distance(PreparedDiagram(medium_first),
                                               PreparedDiagram(medium_second), adaptive,
                                               &medium_stats),
               exact_shift, "adaptive medium mandatory sparse distance route");
  if (medium_stats.router_mandatory_routes != 1 ||
      medium_stats.mandatory_vertices == 0 ||
      medium_stats.optional_pairs_pruned == 0) {
    fail("adaptive dispatcher did not select the medium mandatory sparse route");
  }
}

void batch_cases() {
  const PreparedDiagram query(Diagram{{0.0, 1.0}, {2.0, 3.5}});
  const std::vector<PreparedDiagram> targets{
      PreparedDiagram(Diagram{}),
      PreparedDiagram(Diagram{{0.0, 1.0}}),
      PreparedDiagram(Diagram{{0.25, 1.25}, {2.1, 3.6}}),
  };
  const auto allocated = bottleneck::bottleneck_distances(query, targets);
  std::vector<double> caller_buffer(targets.size(), -1.0);
  bottleneck::bottleneck_distances(query, targets, caller_buffer);
  for (std::size_t index = 0; index < targets.size(); ++index) {
    const double expected = bottleneck::bottleneck_distance(query, targets[index]);
    expect_equal(allocated[index], expected, "allocated batch result");
    expect_equal(caller_buffer[index], expected, "caller-provided batch result");
  }
  bool rejected_bad_size = false;
  try {
    std::vector<double> bad_buffer(targets.size() - 1);
    bottleneck::bottleneck_distances(query, targets, bad_buffer);
  } catch (const std::invalid_argument&) {
    rejected_bad_size = true;
  }
  if (!rejected_bad_size) {
    fail("batch must reject an output span with the wrong size");
  }
}

void lower_bound_hint_cases() {
  const PreparedDiagram first(Diagram{{0.0, 2.0}, {4.0, 4.5}});
  const PreparedDiagram second(Diagram{{0.25, 2.25}, {4.0, 4.5}});
  const double expected = bottleneck::bottleneck_distance(first, second);

  SolverConfig exact_hint;
  exact_hint.lower_bound_hint = expected;
  SolverStats exact_stats;
  expect_equal(bottleneck::bottleneck_distance(first, second, exact_hint, &exact_stats),
               expected, "accepted exact lower-bound hint");
  if (exact_stats.lower_bound_hints_tested != 1 ||
      exact_stats.lower_bound_hints_accepted != 1) {
    fail("exact lower-bound hint was not verified and accepted");
  }

  SolverConfig invalid_hint;
  invalid_hint.lower_bound_hint = expected * 4.0;
  SolverStats invalid_stats;
  expect_equal(bottleneck::bottleneck_distance(first, second, invalid_hint, &invalid_stats),
               expected, "rejected invalid lower-bound hint");
  if (invalid_stats.lower_bound_hints_tested != 1 ||
      invalid_stats.lower_bound_hints_accepted != 0) {
    fail("invalid lower-bound hint was not safely rejected");
  }
}

}  // namespace

int main() {
  deterministic_cases();
  randomized_differential();
  geometric_differential();
  geometric_floating_boundaries();
  multiplicity_flow_differential();
  mandatory_partial_flow_differential();
  component_structure_stats();
  router_x_window_stats();
  adaptive_dispatcher_differential();
  batch_cases();
  lower_bound_hint_cases();
  std::cout << "bottleneck_core_tests: all checks passed\n";
  return 0;
}
