#include <bottleneck/core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string_view>
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
using bottleneck::SolverStats;
using bottleneck::ThresholdStrategy;
using bottleneck::VertexOrder;

struct Variant {
  const char* name;
  SolverConfig config;
};

struct Options {
  std::size_t repetitions = 3;
  std::size_t rounds = 7;
  std::size_t min_points = 128;
  std::size_t max_points = 1024;
  std::size_t block_size = 0;
};

Options parse_options(int argc, char** argv) {
  Options result;
  for (int index = 1; index + 1 < argc; index += 2) {
    const std::string_view option(argv[index]);
    if (option == "--repetitions") {
      result.repetitions = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else if (option == "--rounds") {
      result.rounds = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else if (option == "--min-points") {
      result.min_points = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else if (option == "--max-points") {
      result.max_points = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else if (option == "--block-size") {
      result.block_size = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else {
      throw std::invalid_argument("unknown benchmark option");
    }
  }
  if (result.repetitions == 0 || result.rounds == 0 || result.min_points == 0 ||
      result.max_points < result.min_points) {
    throw std::invalid_argument("invalid benchmark option value");
  }
  return result;
}

Diagram block_diagram(std::size_t size, std::size_t block_size, bool second) {
  const double step = std::ldexp(1.0, -12);
  const double shift = step / 2.0;
  Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    const std::size_t block = index / block_size;
    const std::size_t within = index % block_size;
    const double birth = static_cast<double>(block) * 4.0 +
                         static_cast<double>(within) * step +
                         (second ? shift : 0.0);
    result.push_back({birth, birth + 1.0});
  }
  return result;
}

double percentile(std::vector<double> values, double quantile) {
  std::sort(values.begin(), values.end());
  const std::size_t index = static_cast<std::size_t>(
      std::ceil(quantile * static_cast<double>(values.size())) - 1.0);
  return values[(std::min)(index, values.size() - 1)];
}

}  // namespace

int main(int argc, char** argv) {
  const Options options = parse_options(argc, argv);
  const SolverConfig reference{CandidateStrategy::sort_unique_clipped,
                               ThresholdStrategy::geometric_refinement,
                               DistanceStrategy::recompute_soa,
                               AdjacencyStrategy::on_demand,
                               MatcherStrategy::geometric_hopcroft_karp,
                               VertexOrder::natural};
  const SolverConfig component_config{CandidateStrategy::x_sweep_clipped,
                                      ThresholdStrategy::quickselect,
                                      DistanceStrategy::recompute_soa,
                                      AdjacencyStrategy::x_sweep_csr,
                                      MatcherStrategy::component_kuhn,
                                      VertexOrder::degree_ascending};
  const std::vector<Variant> variants{
      {"adaptive", {}},
      {"geometric_refined", reference},
      {"reusable_kuhn", {CandidateStrategy::x_sweep_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::recompute_soa,
       AdjacencyStrategy::x_sweep_csr, MatcherStrategy::reusable_greedy_kuhn,
       VertexOrder::degree_ascending}},
      {"augmented_component", component_config},
      {"mandatory_sparse", {CandidateStrategy::x_sweep_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::recompute_soa,
       AdjacencyStrategy::on_demand, MatcherStrategy::mandatory_sparse_flow,
       VertexOrder::natural}},
  };
  const std::vector<std::size_t> sizes{128, 256, 512, 1024, 2048};
  const std::vector<std::size_t> block_sizes{1, 4, 16, 64, 256, 2048};
  std::mt19937_64 generator(0xC09E07ULL);
  double sink = 0.0;

  std::cout << "n,block_size,expected_cross_components,exact_components,largest_component_vertices,"
               "s_max,variant,repetitions,rounds,median_us,p95_us,distance,decisions,"
               "adjacency_checks,emitted_edges,augment_searches\n";
  for (std::size_t size : sizes) {
    if (size < options.min_points || size > options.max_points) {
      continue;
    }
    for (std::size_t block_size : block_sizes) {
      if (block_size > size || size % block_size != 0 ||
          (options.block_size != 0 && block_size != options.block_size)) {
        continue;
      }
      const PreparedDiagram first(block_diagram(size, block_size, false));
      const PreparedDiagram second(block_diagram(size, block_size, true));
      const double expected = bottleneck::bottleneck_distance(first, second, reference);
      SolverStats structure_stats;
      if (!bottleneck::bottleneck_within(first, second, expected,
                                         component_config, &structure_stats)) {
        std::cerr << "component structure rejected exact threshold\n";
        return 2;
      }
      const double augmented_vertices = static_cast<double>(4 * size);
      const double s_max = static_cast<double>(structure_stats.largest_component_vertices) /
                           augmented_vertices;
      std::vector<std::vector<double>> samples(variants.size());
      std::vector<SolverStats> stats(variants.size());
      std::vector<std::size_t> order(variants.size());
      std::iota(order.begin(), order.end(), 0);
      for (const Variant& variant : variants) {
        const double value = bottleneck::bottleneck_distance(first, second, variant.config);
        if (value != expected) {
          std::cerr << "warmup mismatch: n=" << size << ", block=" << block_size
                    << ", variant=" << variant.name << '\n';
          return 2;
        }
        sink += value;
      }
      for (std::size_t round = 0; round < options.rounds; ++round) {
        std::shuffle(order.begin(), order.end(), generator);
        for (std::size_t variant_index : order) {
          const auto start = Clock::now();
          double value = 0.0;
          for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
            value = bottleneck::bottleneck_distance(first, second,
                                                    variants[variant_index].config,
                                                    &stats[variant_index]);
            sink += value;
          }
          const auto stop = Clock::now();
          if (value != expected) {
            std::cerr << "mismatch: n=" << size << ", block=" << block_size
                      << ", variant=" << variants[variant_index].name << '\n';
            return 2;
          }
          samples[variant_index].push_back(
              std::chrono::duration<double, std::micro>(stop - start).count() /
              static_cast<double>(options.repetitions));
        }
      }
      const double calls = static_cast<double>(options.repetitions * options.rounds);
      for (std::size_t index = 0; index < variants.size(); ++index) {
        std::cout << size << ',' << block_size << ',' << size / block_size << ','
                  << structure_stats.component_count << ','
                  << structure_stats.largest_component_vertices << ','
                  << std::fixed << std::setprecision(6) << s_max << ','
                  << variants[index].name << ',' << options.repetitions << ','
                  << options.rounds << ',' << std::setprecision(3)
                  << percentile(samples[index], 0.5) << ','
                  << percentile(samples[index], 0.95) << ',' << expected << ','
                  << static_cast<double>(stats[index].threshold_decisions) / calls << ','
                  << static_cast<double>(stats[index].adjacency_checks) / calls << ','
                  << static_cast<double>(stats[index].emitted_edges) / calls << ','
                  << static_cast<double>(stats[index].augment_searches) / calls << '\n';
      }
    }
  }
  if (sink == -1.0) {
    std::cerr << "unreachable sink\n";
  }
  return 0;
}
