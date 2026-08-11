#include <bottleneck/core.hpp>

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
                                       CandidateStrategy::sort_unique_greedy_clipped}) {
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
                                          MatcherStrategy::adaptive,
                                          MatcherStrategy::hopcroft_karp,
                                          MatcherStrategy::greedy_hopcroft_karp}) {
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

}  // namespace

int main() {
  deterministic_cases();
  randomized_differential();
  batch_cases();
  std::cout << "bottleneck_core_tests: all checks passed\n";
  return 0;
}
