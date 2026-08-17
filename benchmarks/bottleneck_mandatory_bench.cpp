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
#include <string>
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
  std::size_t rounds = 5;
  std::size_t min_points = 128;
  std::size_t max_points = 2048;
  std::size_t ratio = 0;
  double optional_ratio = -1.0;
  std::string pattern = "all";
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
    } else if (option == "--ratio") {
      result.ratio = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else if (option == "--optional-ratio") {
      result.optional_ratio = std::stod(argv[index + 1]);
    } else if (option == "--pattern") {
      result.pattern = argv[index + 1];
    } else {
      throw std::invalid_argument("unknown benchmark option");
    }
  }
  if (result.repetitions == 0 || result.rounds == 0 || result.min_points == 0 ||
      result.max_points < result.min_points ||
      (result.ratio != 0 && result.ratio != 1 && result.ratio != 2 &&
       result.ratio != 4 && result.ratio != 8 && result.ratio != 16) ||
      (result.optional_ratio >= 0.0 && result.optional_ratio > 1.0)) {
    throw std::invalid_argument("invalid benchmark option value");
  }
  return result;
}

Diagram paired_near_diagonal(std::size_t small_size, std::size_t large_size,
                             double optional_ratio, bool second) {
  constexpr double shift = 0.01;
  Diagram result;
  result.reserve(second ? large_size : small_size);
  const std::size_t common_size = small_size;
  const std::size_t optional_count = static_cast<std::size_t>(
      std::llround(static_cast<double>(common_size) * optional_ratio));
  for (std::size_t index = 0; index < common_size; ++index) {
    const double position = static_cast<double>(index) /
        static_cast<double>((std::max)(std::size_t{1}, common_size - 1));
    const double birth = -5.0 + 10.0 * position +
                         0.025 * std::sin(static_cast<double>(index) * 1.7);
    const bool optional = index < optional_count;
    const double persistence = optional
                                   ? 0.012 + 0.004 * (0.5 + 0.5 * std::sin(
                                                            static_cast<double>(index) * 0.73))
                                   : 0.5 + 1.5 * (0.5 + 0.5 * std::cos(
                                                         static_cast<double>(index) * 0.37));
    const double offset = second ? shift : 0.0;
    result.push_back({birth + offset, birth + persistence + offset});
  }
  if (second) {
    for (std::size_t index = common_size; index < large_size; ++index) {
      const double position = static_cast<double>(index - common_size + 1) /
          static_cast<double>(large_size - common_size + 1);
      const double birth = -5.0 + 10.0 * position +
                           0.017 * std::cos(static_cast<double>(index) * 0.91);
      const double persistence = 0.008 +
          0.006 * (0.5 + 0.5 * std::sin(static_cast<double>(index) * 0.43));
      result.push_back({birth, birth + persistence});
    }
  }
  return result;
}

Diagram generated_diagram(std::string_view pattern, std::size_t size, bool second,
                          std::uint64_t seed) {
  std::mt19937_64 generator(seed + (second ? 0x9E3779B97F4A7C15ULL : 0));
  Diagram result;
  result.reserve(size);
  if (pattern == "uniform") {
    std::uniform_real_distribution<double> birth(-5.0, 5.0);
    std::uniform_real_distribution<double> persistence(0.05, 2.0);
    for (std::size_t index = 0; index < size; ++index) {
      const double point_birth = birth(generator);
      result.push_back({point_birth, point_birth + persistence(generator)});
    }
  } else if (pattern == "clustered") {
    std::normal_distribution<double> offset(0.0, 0.18);
    std::uniform_real_distribution<double> persistence(0.2, 1.8);
    const double centers[] = {-3.0, 0.0, 3.0};
    for (std::size_t index = 0; index < size; ++index) {
      const double point_birth = centers[index % 3] + offset(generator);
      result.push_back({point_birth, point_birth + persistence(generator)});
    }
  } else if (pattern == "separated") {
    std::uniform_real_distribution<double> offset(0.0, 2.0);
    std::uniform_real_distribution<double> persistence(0.05, 2.0);
    const double origin = second ? 10.0 : -12.0;
    for (std::size_t index = 0; index < size; ++index) {
      const double point_birth = origin + offset(generator);
      result.push_back({point_birth, point_birth + persistence(generator)});
    }
  } else {
    throw std::invalid_argument("unknown benchmark pattern");
  }
  return result;
}

double percentile(std::vector<double> values, double quantile) {
  std::sort(values.begin(), values.end());
  const std::size_t index = static_cast<std::size_t>(
      std::ceil(quantile * static_cast<double>(values.size())) - 1.0);
  return values[(std::min)(index, values.size() - 1)];
}

double optional_fraction(const PreparedDiagram& first, const PreparedDiagram& second,
                         double threshold) {
  std::size_t optional = 0;
  for (double diagonal : first.finite_diagonal_distances()) {
    optional += static_cast<std::size_t>(diagonal <= threshold);
  }
  for (double diagonal : second.finite_diagonal_distances()) {
    optional += static_cast<std::size_t>(diagonal <= threshold);
  }
  const std::size_t total = first.finite_points().size() + second.finite_points().size();
  return total == 0 ? 1.0 : static_cast<double>(optional) / static_cast<double>(total);
}

bool selected(std::string_view requested, std::string_view value) {
  return requested == "all" || requested == value;
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
  const std::vector<Variant> variants{
      {"adaptive", {}},
      {"legacy_full_scan", {CandidateStrategy::sort_unique_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::dense_aos,
       AdjacencyStrategy::adaptive, MatcherStrategy::reusable_greedy_kuhn,
       VertexOrder::degree_ascending}},
      {"x_sweep_only", {CandidateStrategy::x_sweep_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::dense_aos,
       AdjacencyStrategy::adaptive, MatcherStrategy::reusable_greedy_kuhn,
       VertexOrder::degree_ascending}},
      {"mandatory_dense", {CandidateStrategy::sort_unique_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::recompute_soa,
       AdjacencyStrategy::on_demand, MatcherStrategy::mandatory_flow,
       VertexOrder::natural}},
      {"mandatory_sparse", {CandidateStrategy::sort_unique_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::recompute_soa,
       AdjacencyStrategy::on_demand, MatcherStrategy::mandatory_sparse_flow,
       VertexOrder::natural}},
      {"e7_combined", {CandidateStrategy::x_sweep_clipped,
       ThresholdStrategy::quickselect, DistanceStrategy::recompute_soa,
       AdjacencyStrategy::on_demand, MatcherStrategy::mandatory_sparse_flow,
       VertexOrder::natural}},
  };
  const std::vector<std::size_t> sizes{128, 256, 512, 1024, 2048};
  const std::vector<std::size_t> ratios{1, 2, 4, 8, 16};
  const std::vector<std::string_view> patterns{
      "near_diagonal", "uniform", "clustered", "separated"};
  const std::vector<double> optional_ratios{0.5, 0.75, 0.9, 1.0};
  std::mt19937_64 order_generator(0xE7B0771EULL);
  double sink = 0.0;

  std::cout << "small_n,large_m,cardinality_ratio,pattern,target_optional_ratio,"
               "exact_optional_ratio,variant,repetitions,rounds,median_us,p95_us,distance,"
               "raw_candidates_per_call,x_window_candidates_per_call,mandatory_vertices_per_call,"
               "optional_pairs_pruned_per_call,capacity_edges_per_call,adjacency_checks_per_call,"
               "decisions_per_call\n";
  for (std::size_t small_size : sizes) {
    if (small_size < options.min_points || small_size > options.max_points) {
      continue;
    }
    for (std::size_t ratio : ratios) {
      if ((options.ratio != 0 && ratio != options.ratio) ||
          small_size > options.max_points / ratio) {
        continue;
      }
      const std::size_t large_size = small_size * ratio;
      for (std::string_view pattern : patterns) {
        if (!selected(options.pattern, pattern)) {
          continue;
        }
        const std::vector<double> case_optional_ratios =
            pattern == "near_diagonal" ? optional_ratios : std::vector<double>{-1.0};
        for (double target_optional_ratio : case_optional_ratios) {
          if (options.optional_ratio >= 0.0 &&
              (target_optional_ratio < 0.0 || target_optional_ratio != options.optional_ratio)) {
            continue;
          }
          Diagram first_diagram;
          Diagram second_diagram;
          if (pattern == "near_diagonal") {
            first_diagram = paired_near_diagonal(small_size, large_size,
                                                 target_optional_ratio, false);
            second_diagram = paired_near_diagonal(small_size, large_size,
                                                  target_optional_ratio, true);
          } else {
            const std::uint64_t seed = 0xB0771EULL + small_size * 131 + ratio * 17;
            first_diagram = generated_diagram(pattern, small_size, false, seed);
            second_diagram = generated_diagram(pattern, large_size, true, seed);
          }
          const PreparedDiagram first(first_diagram);
          const PreparedDiagram second(second_diagram);
          const double expected = bottleneck::bottleneck_distance(first, second, reference);
          const double exact_optional_ratio = optional_fraction(first, second, expected);
          std::vector<std::vector<double>> samples(variants.size());
          std::vector<SolverStats> stats(variants.size());
          std::vector<std::size_t> order(variants.size());
          std::iota(order.begin(), order.end(), 0);
          for (const Variant& variant : variants) {
            const double value = bottleneck::bottleneck_distance(first, second, variant.config);
            if (value != expected) {
              std::cerr << "warmup mismatch: pattern=" << pattern << ", n=" << small_size
                        << ", m=" << large_size << ", variant=" << variant.name << '\n';
              return 2;
            }
            sink += value;
          }
          for (std::size_t round = 0; round < options.rounds; ++round) {
            std::shuffle(order.begin(), order.end(), order_generator);
            for (std::size_t variant_index : order) {
              const auto start = Clock::now();
              double value = 0.0;
              for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
                value = bottleneck::bottleneck_distance(
                    first, second, variants[variant_index].config, &stats[variant_index]);
                sink += value;
              }
              const auto stop = Clock::now();
              if (value != expected) {
                std::cerr << "mismatch: pattern=" << pattern << ", n=" << small_size
                          << ", m=" << large_size << ", variant="
                          << variants[variant_index].name << '\n';
                return 2;
              }
              samples[variant_index].push_back(
                  std::chrono::duration<double, std::micro>(stop - start).count() /
                  static_cast<double>(options.repetitions));
            }
          }
          const double calls = static_cast<double>(options.repetitions * options.rounds);
          for (std::size_t index = 0; index < variants.size(); ++index) {
            const SolverStats& current = stats[index];
            std::cout << small_size << ',' << large_size << ',' << ratio << ',' << pattern << ','
                      << target_optional_ratio << ',' << exact_optional_ratio << ','
                      << variants[index].name << ',' << options.repetitions << ','
                      << options.rounds << ',' << std::fixed << std::setprecision(3)
                      << percentile(samples[index], 0.5) << ','
                      << percentile(samples[index], 0.95) << ',' << expected << ','
                      << static_cast<double>(current.raw_candidates) / calls << ','
                      << static_cast<double>(current.x_window_candidates) / calls << ','
                      << static_cast<double>(current.mandatory_vertices) / calls << ','
                      << static_cast<double>(current.optional_pairs_pruned) / calls << ','
                      << static_cast<double>(current.capacity_edges) / calls << ','
                      << static_cast<double>(current.adjacency_checks) / calls << ','
                      << static_cast<double>(current.threshold_decisions) / calls << '\n';
          }
        }
      }
    }
  }
  if (sink == -1.0) {
    std::cerr << "unreachable sink\n";
  }
  return 0;
}
