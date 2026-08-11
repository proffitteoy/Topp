#include <bottleneck/core.hpp>

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

bottleneck::Diagram random_diagram(std::mt19937_64& generator, std::size_t size) {
  std::uniform_real_distribution<double> birth(-2.0, 2.0);
  std::uniform_real_distribution<double> persistence(0.01, 2.0);
  bottleneck::Diagram result;
  result.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    const double b = birth(generator);
    result.push_back({b, b + persistence(generator)});
  }
  return result;
}

std::size_t parse_repetitions(int argc, char** argv) {
  if (argc == 3 && std::string_view(argv[1]) == "--repetitions") {
    return static_cast<std::size_t>(std::stoul(argv[2]));
  }
  return 20;
}

}  // namespace

int main(int argc, char** argv) {
  const std::size_t repetitions = parse_repetitions(argc, argv);
  std::mt19937_64 generator(0xB0771EULL);
  const std::vector<std::size_t> sizes{4, 8, 16, 32, 64};
  const std::vector<bottleneck::SolverConfig> configs{
      {bottleneck::CandidateStrategy::sort_all, bottleneck::ThresholdStrategy::binary,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::on_demand, bottleneck::MatcherStrategy::kuhn,
       bottleneck::VertexOrder::natural},
      {bottleneck::CandidateStrategy::sort_unique, bottleneck::ThresholdStrategy::binary,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::hopcroft_karp,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::binary,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::hopcroft_karp,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_greedy_clipped,
       bottleneck::ThresholdStrategy::gudhi_alpha, bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::constraint_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::component_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::sparse_csr, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::bitset64, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::bitset128, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::bitset256, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_soa,
       bottleneck::AdjacencyStrategy::bitset64, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::recompute_aos,
       bottleneck::AdjacencyStrategy::bitset64, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::recompute_soa,
       bottleneck::AdjacencyStrategy::bitset64, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::hopcroft_karp,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::gudhi_alpha,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset,
       bottleneck::MatcherStrategy::greedy_hopcroft_karp,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::exponential,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::hopcroft_karp,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped, bottleneck::ThresholdStrategy::incremental,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::incremental_blocked,
       bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect, bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::dynamic_bitset, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect, bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::bitset_auto, bottleneck::MatcherStrategy::greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect, bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::bitset_auto,
       bottleneck::MatcherStrategy::reusable_greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect, bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::x_sweep_csr,
       bottleneck::MatcherStrategy::reusable_greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect, bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::adaptive,
       bottleneck::MatcherStrategy::reusable_greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect,
       bottleneck::DistanceStrategy::dense_soa_avx2,
       bottleneck::AdjacencyStrategy::adaptive,
       bottleneck::MatcherStrategy::reusable_greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect, bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::on_demand,
       bottleneck::MatcherStrategy::mandatory_flow,
       bottleneck::VertexOrder::natural},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect, bottleneck::DistanceStrategy::dense_aos,
       bottleneck::AdjacencyStrategy::adaptive,
       bottleneck::MatcherStrategy::fixed_greedy_kuhn,
       bottleneck::VertexOrder::degree_ascending},
      {bottleneck::CandidateStrategy::sort_unique_clipped,
       bottleneck::ThresholdStrategy::quickselect,
       bottleneck::DistanceStrategy::dense_soa_avx2,
       bottleneck::AdjacencyStrategy::adaptive,
       bottleneck::MatcherStrategy::adaptive,
       bottleneck::VertexOrder::degree_ascending},
      {},
  };

  std::cout << "size,candidates,threshold,distance,adjacency,matcher,order,repetitions,total_us,mean_us,"
               "raw_candidates,retained_candidates,clipped_candidates,decisions,adjacency_checks,"
               "emitted_edges\n";
  double sink = 0.0;
  for (std::size_t size : sizes) {
    std::vector<std::pair<bottleneck::PreparedDiagram, bottleneck::PreparedDiagram>> pairs;
    pairs.reserve(repetitions);
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
      pairs.emplace_back(random_diagram(generator, size), random_diagram(generator, size));
    }

    double reference_sum = 0.0;
    for (const auto& pair : pairs) {
      reference_sum += bottleneck::bottleneck_distance(pair.first, pair.second, configs.front());
    }

    for (const auto& config : configs) {
      bottleneck::SolverStats stats;
      const auto start = Clock::now();
      double sum = 0.0;
      for (const auto& pair : pairs) {
        sum += bottleneck::bottleneck_distance(pair.first, pair.second, config, &stats);
      }
      const auto stop = Clock::now();
      if (sum != reference_sum) {
        std::cerr << "variant result mismatch at size " << size << '\n';
        return 2;
      }
      sink += sum;
      const double total_us =
          std::chrono::duration<double, std::micro>(stop - start).count();
      std::cout << size << ',' << bottleneck::to_string(config.candidates) << ','
                << bottleneck::to_string(config.threshold) << ','
                << bottleneck::to_string(config.distance) << ','
                << bottleneck::to_string(config.adjacency) << ','
                << bottleneck::to_string(config.matcher) << ','
                << bottleneck::to_string(config.vertex_order) << ',' << repetitions << ','
                << std::fixed << std::setprecision(3) << total_us << ','
                << total_us / static_cast<double>(repetitions) << ',' << stats.raw_candidates << ','
                << stats.retained_candidates << ',' << stats.clipped_candidates << ','
                << stats.threshold_decisions << ','
                << stats.adjacency_checks << ',' << stats.emitted_edges << '\n';
    }
  }
  if (sink == -1.0) {
    std::cerr << "unreachable sink value\n";
  }
  return 0;
}
