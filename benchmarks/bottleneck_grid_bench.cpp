#include <bottleneck/core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using bottleneck::AdjacencyStrategy;
using bottleneck::CandidateStrategy;
using bottleneck::Diagram;
using bottleneck::DistanceStrategy;
using bottleneck::MatcherStrategy;
using bottleneck::PreparedDiagram;
using bottleneck::SolverConfig;
using bottleneck::ThresholdStrategy;
using bottleneck::VertexOrder;

enum class Pattern { uniform, near_diagonal, clustered, repeated, separated };

struct Workload {
  Pattern pattern;
  std::size_t first_size;
  std::size_t second_size;
};

const char* pattern_name(Pattern pattern) {
  switch (pattern) {
    case Pattern::uniform:
      return "uniform";
    case Pattern::near_diagonal:
      return "near_diagonal";
    case Pattern::clustered:
      return "clustered";
    case Pattern::repeated:
      return "repeated";
    case Pattern::separated:
      return "separated";
  }
  return "unknown";
}

Diagram generate_diagram(std::mt19937_64& generator, Pattern pattern, std::size_t size,
                         bool second) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_real_distribution<double> signed_unit(-1.0, 1.0);
  Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    double birth = 0.0;
    double persistence = 0.0;
    switch (pattern) {
      case Pattern::uniform:
        birth = -2.0 + 4.0 * unit(generator);
        persistence = 0.01 + 1.99 * unit(generator);
        break;
      case Pattern::near_diagonal:
        birth = -2.0 + 4.0 * unit(generator);
        persistence = 0.0001 + 0.0499 * unit(generator);
        break;
      case Pattern::clustered: {
        const double center = static_cast<double>(static_cast<int>(index % 4) - 2) * 0.6;
        birth = center + 0.03 * signed_unit(generator) + (second ? 0.005 : 0.0);
        persistence = 0.45 + 0.08 * static_cast<double>(index % 3) +
                      0.015 * signed_unit(generator);
        break;
      }
      case Pattern::repeated:
        birth = static_cast<double>(static_cast<int>(index % 5) - 2) * 0.4;
        persistence = 0.25 + 0.15 * static_cast<double>(index % 4);
        break;
      case Pattern::separated:
        birth = (second ? 8.0 : -10.0) + 2.0 * unit(generator);
        persistence = 0.1 + 1.9 * unit(generator);
        break;
    }
    result.push_back({birth, birth + persistence});
  }
  return result;
}

double cross_density(const PreparedDiagram& first, const PreparedDiagram& second,
                     double threshold) {
  const auto& first_points = first.finite_points();
  const auto& second_points = second.finite_points();
  if (first_points.empty() || second_points.empty()) {
    return 0.0;
  }
  std::size_t edges = 0;
  for (const auto& left : first_points) {
    for (const auto& right : second_points) {
      const double distance = (std::max)(std::fabs(left.birth - right.birth),
                                         std::fabs(left.death - right.death));
      edges += static_cast<std::size_t>(distance <= threshold);
    }
  }
  return static_cast<double>(edges) /
         static_cast<double>(first_points.size() * second_points.size());
}

struct BenchmarkOptions {
  std::size_t repetitions = 20;
  std::size_t rounds = 1;
};

BenchmarkOptions parse_options(int argc, char** argv) {
  BenchmarkOptions options;
  for (int index = 1; index + 1 < argc; index += 2) {
    const std::string_view name(argv[index]);
    if (name == "--repetitions") {
      options.repetitions = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else if (name == "--rounds") {
      options.rounds = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else {
      throw std::invalid_argument("unknown benchmark option");
    }
  }
  return options;
}

double percentile(std::vector<double> values, double fraction) {
  std::sort(values.begin(), values.end());
  const std::size_t index = static_cast<std::size_t>(
      std::ceil(fraction * static_cast<double>(values.size()))) - 1;
  return values[(std::min)(index, values.size() - 1)];
}

}  // namespace

int main(int argc, char** argv) {
  const BenchmarkOptions options = parse_options(argc, argv);
  const SolverConfig reference{CandidateStrategy::sort_all, ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn, VertexOrder::natural};
  const std::vector<SolverConfig> configs{
      reference,
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::sparse_csr,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::bitset64,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::bitset128,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::bitset256,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::binary,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::hopcroft_karp, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::greedy_hopcroft_karp, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::constraint_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::component_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::incremental_blocked,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::bitset64,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::bitset128,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::bitset256,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::bitset_auto,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::bitset_auto,
       MatcherStrategy::reusable_greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::x_sweep_csr,
       MatcherStrategy::reusable_greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::adaptive,
       MatcherStrategy::reusable_greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_soa_avx2, AdjacencyStrategy::adaptive,
       MatcherStrategy::reusable_greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::on_demand,
       MatcherStrategy::mandatory_flow, VertexOrder::natural},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_aos, AdjacencyStrategy::adaptive,
       MatcherStrategy::fixed_greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
       DistanceStrategy::dense_soa_avx2, AdjacencyStrategy::adaptive,
       MatcherStrategy::adaptive, VertexOrder::degree_ascending},
      {},
      {CandidateStrategy::sort_unique_greedy_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::greedy_kuhn, VertexOrder::degree_ascending},
      {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::gudhi_alpha,
       DistanceStrategy::dense_aos, AdjacencyStrategy::dynamic_bitset,
       MatcherStrategy::greedy_kuhn, VertexOrder::natural},
  };
  const std::vector<Workload> workloads{
      {Pattern::uniform, 8, 8},       {Pattern::uniform, 16, 16},
      {Pattern::uniform, 32, 32},     {Pattern::uniform, 64, 64},
      {Pattern::uniform, 128, 128},   {Pattern::uniform, 8, 64},
      {Pattern::uniform, 256, 256},
      {Pattern::uniform, 64, 8},      {Pattern::uniform, 16, 128},
      {Pattern::uniform, 128, 16},    {Pattern::uniform, 32, 256},
      {Pattern::uniform, 256, 32},    {Pattern::near_diagonal, 32, 32},
      {Pattern::near_diagonal, 64, 64}, {Pattern::clustered, 32, 32},
      {Pattern::near_diagonal, 256, 256},
      {Pattern::clustered, 64, 64},   {Pattern::repeated, 32, 32},
      {Pattern::clustered, 256, 256}, {Pattern::repeated, 64, 64},
      {Pattern::repeated, 256, 256},  {Pattern::separated, 32, 32},
      {Pattern::separated, 64, 64},   {Pattern::separated, 256, 256},
  };

  std::cout << "pattern,n,m,candidates,threshold,distance,adjacency,matcher,order,repetitions,rounds,"
               "median_us,p95_us,relative_to_reference,mean_cross_density,decisions,clipped_candidates,"
               "forced_matches,components,component_rejects\n";
  std::mt19937_64 generator(0xD15A7C4ULL);
  std::mt19937_64 order_generator(0x0A11CEULL);
  double sink = 0.0;
  for (const Workload& workload : workloads) {
    std::vector<std::pair<PreparedDiagram, PreparedDiagram>> pairs;
    pairs.reserve(options.repetitions);
    for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
      pairs.emplace_back(generate_diagram(generator, workload.pattern, workload.first_size, false),
                         generate_diagram(generator, workload.pattern, workload.second_size, true));
    }

    std::vector<double> expected;
    expected.reserve(options.repetitions);
    double density_sum = 0.0;
    for (const auto& pair : pairs) {
      const double value = bottleneck::bottleneck_distance(pair.first, pair.second, reference);
      expected.push_back(value);
      density_sum += cross_density(pair.first, pair.second, value);
    }
    std::vector<std::vector<double>> timings(configs.size());
    std::vector<bottleneck::SolverStats> all_stats(configs.size());
    std::vector<std::size_t> order(configs.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    for (std::size_t round = 0; round < options.rounds; ++round) {
      std::shuffle(order.begin(), order.end(), order_generator);
      for (std::size_t config_index : order) {
        const SolverConfig& config = configs[config_index];
        const auto start = Clock::now();
        for (std::size_t index = 0; index < pairs.size(); ++index) {
          const double value = bottleneck::bottleneck_distance(
              pairs[index].first, pairs[index].second, config, &all_stats[config_index]);
          if (value != expected[index]) {
            std::cerr << "variant mismatch for " << pattern_name(workload.pattern) << ' '
                      << workload.first_size << 'x' << workload.second_size << '\n';
            return 2;
          }
          sink += value;
        }
        const auto stop = Clock::now();
        timings[config_index].push_back(
            std::chrono::duration<double, std::micro>(stop - start).count() /
            static_cast<double>(options.repetitions));
      }
    }
    const double reference_median = percentile(timings.front(), 0.5);
    for (std::size_t config_index = 0; config_index < configs.size(); ++config_index) {
      const SolverConfig& config = configs[config_index];
      const bottleneck::SolverStats& stats = all_stats[config_index];
      const double median_us = percentile(timings[config_index], 0.5);
      const double p95_us = percentile(timings[config_index], 0.95);
      std::cout << pattern_name(workload.pattern) << ',' << workload.first_size << ','
                << workload.second_size << ',' << bottleneck::to_string(config.candidates) << ','
                << bottleneck::to_string(config.threshold) << ','
                << bottleneck::to_string(config.distance) << ','
                << bottleneck::to_string(config.adjacency) << ','
                << bottleneck::to_string(config.matcher) << ','
                << bottleneck::to_string(config.vertex_order) << ',' << options.repetitions << ','
                << options.rounds << ',' << std::fixed << std::setprecision(3) << median_us << ','
                << p95_us << ',' << median_us / reference_median << ','
                << density_sum / static_cast<double>(options.repetitions) << ','
                << stats.threshold_decisions << ',' << stats.clipped_candidates << ','
                << stats.forced_matches << ',' << stats.component_count << ','
                << stats.component_rejects << '\n';
    }
  }
  if (sink == -1.0) {
    std::cerr << "unreachable sink value\n";
  }
  return 0;
}
