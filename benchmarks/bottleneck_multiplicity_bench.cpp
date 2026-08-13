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

struct Variant { const char* name; SolverConfig config; };
struct Options {
  std::size_t repetitions = 3;
  std::size_t rounds = 7;
  std::size_t max_points = 2048;
  double min_duplicate_ratio = 0.0;
};

Options parse_options(int argc, char** argv) {
  Options result;
  for (int index = 1; index + 1 < argc; index += 2) {
    const std::string_view option(argv[index]);
    if (option == "--repetitions") result.repetitions = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    else if (option == "--rounds") result.rounds = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    else if (option == "--max-points") result.max_points = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    else if (option == "--min-duplicate-ratio") result.min_duplicate_ratio = std::stod(argv[index + 1]);
    else throw std::invalid_argument("unknown benchmark option");
  }
  return result;
}

Diagram duplicate_diagram(std::size_t raw_size, double duplicate_ratio, bool second) {
  const std::size_t unique_size = (std::max)(std::size_t{1}, static_cast<std::size_t>(
      std::llround(static_cast<double>(raw_size) * (1.0 - duplicate_ratio))));
  Diagram unique;
  unique.reserve(unique_size);
  for (std::size_t index = 0; index < unique_size; ++index) {
    const double position = static_cast<double>(index) /
        static_cast<double>((std::max)(std::size_t{1}, unique_size - 1));
    const double birth = -3.0 + 6.0 * position +
        0.07 * std::sin(static_cast<double>(index) * 1.7) + (second ? 0.004 : 0.0);
    const double persistence = 0.25 + 1.5 * (0.5 + 0.5 * std::cos(
        static_cast<double>(index) * 0.73 + (second ? 0.03 : 0.0)));
    unique.push_back({birth, birth + persistence});
  }
  Diagram result;
  result.reserve(raw_size);
  for (std::size_t index = 0; index < raw_size; ++index) result.push_back(unique[index % unique_size]);
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
      ThresholdStrategy::geometric_refinement, DistanceStrategy::recompute_soa,
      AdjacencyStrategy::on_demand, MatcherStrategy::geometric_hopcroft_karp,
      VertexOrder::natural};
  const std::vector<Variant> variants{
      {"adaptive", {}},
      {"geometric_default", reference},
      {"multiplicity_flow", {CandidateStrategy::sort_unique_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::recompute_soa,
       AdjacencyStrategy::on_demand, MatcherStrategy::multiplicity_flow,
       VertexOrder::natural}},
      {"legacy_quickselect", {CandidateStrategy::sort_unique_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::dense_aos,
       AdjacencyStrategy::adaptive, MatcherStrategy::adaptive,
       VertexOrder::degree_ascending}},
  };
  const std::vector<std::size_t> sizes{128, 512, 2048};
  const std::vector<double> duplicate_ratios{
      0.0, 0.1, 0.25, 0.5, 0.75, 0.8, 0.85, 0.875, 0.9, 0.95, 0.99};
  std::mt19937_64 generator(0xD0B1CA7EULL);
  double sink = 0.0;
  std::cout << "n,duplicate_ratio,unique_first,unique_second,variant,repetitions,rounds,median_us,p95_us,distance,decisions,groups,points_removed,capacity_edges\n";
  for (std::size_t size : sizes) {
    if (size > options.max_points) continue;
    for (double ratio : duplicate_ratios) {
      if (ratio < options.min_duplicate_ratio) continue;
      const PreparedDiagram first(duplicate_diagram(size, ratio, false));
      const PreparedDiagram second(duplicate_diagram(size, ratio, true));
      const double expected = bottleneck::bottleneck_distance(first, second, reference);
      std::vector<std::vector<double>> samples(variants.size());
      std::vector<SolverStats> stats(variants.size());
      std::vector<std::size_t> order(variants.size());
      std::iota(order.begin(), order.end(), 0);
      for (const Variant& variant : variants) sink += bottleneck::bottleneck_distance(first, second, variant.config);
      for (std::size_t round = 0; round < options.rounds; ++round) {
        std::shuffle(order.begin(), order.end(), generator);
        for (std::size_t variant_index : order) {
          const auto start = Clock::now();
          double value = 0.0;
          for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
            value = bottleneck::bottleneck_distance(first, second, variants[variant_index].config,
                                                    &stats[variant_index]);
            sink += value;
          }
          const auto stop = Clock::now();
          if (value != expected) {
            std::cerr << "mismatch at n=" << size << ", ratio=" << ratio
                      << ", variant=" << variants[variant_index].name << '\n';
            return 2;
          }
          samples[variant_index].push_back(std::chrono::duration<double, std::micro>(stop - start).count() /
                                           static_cast<double>(options.repetitions));
        }
      }
      for (std::size_t index = 0; index < variants.size(); ++index) {
        std::cout << size << ',' << ratio << ','
                  << first.finite_duplicate_representatives().size() << ','
                  << second.finite_duplicate_representatives().size() << ','
                  << variants[index].name << ',' << options.repetitions << ',' << options.rounds << ','
                  << std::fixed << std::setprecision(3) << percentile(samples[index], 0.5) << ','
                  << percentile(samples[index], 0.95) << ',' << expected << ','
                  << stats[index].threshold_decisions << ',' << stats[index].multiplicity_groups << ','
                  << stats[index].multiplicity_points_removed << ',' << stats[index].capacity_edges << '\n';
      }
    }
  }
  if (sink == -1.0) std::cerr << "unreachable sink\n";
  return 0;
}
