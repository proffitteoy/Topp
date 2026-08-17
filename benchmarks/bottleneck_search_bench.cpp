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
using bottleneck::Point;
using bottleneck::PreparedDiagram;
using bottleneck::SolverConfig;
using bottleneck::SolverStats;
using bottleneck::ThresholdStrategy;
using bottleneck::VertexOrder;

enum class Pattern { uniform, clustered, near_diagonal, threshold_shell };

struct Variant {
  const char* name;
  SolverConfig config;
};

struct Options {
  std::size_t repetitions = 3;
  std::size_t rounds = 7;
  std::size_t min_points = 64;
  std::size_t max_points = 256;
  std::string_view pattern = "all";
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
    } else if (option == "--pattern") {
      result.pattern = argv[index + 1];
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

const char* pattern_name(Pattern pattern) {
  switch (pattern) {
    case Pattern::uniform:
      return "uniform";
    case Pattern::clustered:
      return "clustered";
    case Pattern::near_diagonal:
      return "near_diagonal";
    case Pattern::threshold_shell:
      return "threshold_shell";
  }
  return "unknown";
}

Diagram generate_diagram(std::mt19937_64& generator, Pattern pattern,
                         std::size_t size, bool second) {
  Diagram result;
  result.reserve(size);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 0.15);
  const double shell_radius = std::ldexp(1.0, -10);
  for (std::size_t index = 0; index < size; ++index) {
    double birth = 0.0;
    double death = 0.0;
    if (pattern == Pattern::uniform) {
      birth = -4.0 + 8.0 * unit(generator);
      death = birth + 0.05 + 1.95 * unit(generator);
    } else if (pattern == Pattern::clustered) {
      const double center = static_cast<double>(static_cast<int>(index % 4) - 2) * 1.5;
      birth = center + normal(generator) + (second ? 0.01 : 0.0);
      death = birth + 0.4 + 0.8 * unit(generator);
    } else if (pattern == Pattern::near_diagonal) {
      birth = -4.0 + 8.0 * unit(generator);
      death = birth + 0.0001 + 0.04 * unit(generator);
    } else {
      const double offset = static_cast<double>(index) * shell_radius /
                            static_cast<double>((std::max)(std::size_t{1}, size));
      birth = offset;
      death = 2.0 + (second ? shell_radius : 0.0);
    }
    result.push_back({birth, death});
  }
  return result;
}

struct CandidateShape {
  std::size_t raw = 0;
  std::size_t clipped = 0;
  std::size_t unique = 0;
};

CandidateShape candidate_shape(const PreparedDiagram& first,
                               const PreparedDiagram& second) {
  const double upper = (std::max)(first.max_finite_diagonal_distance(),
                                  second.max_finite_diagonal_distance());
  std::vector<double> values;
  values.reserve(1 + first.finite_points().size() + second.finite_points().size() +
                 first.finite_points().size() * second.finite_points().size());
  values.push_back(0.0);
  values.insert(values.end(), first.finite_diagonal_distances().begin(),
                first.finite_diagonal_distances().end());
  values.insert(values.end(), second.finite_diagonal_distances().begin(),
                second.finite_diagonal_distances().end());
  std::size_t clipped = 0;
  for (const Point& left : first.finite_points()) {
    for (const Point& right : second.finite_points()) {
      const double value = (std::max)(std::fabs(left.birth - right.birth),
                                      std::fabs(left.death - right.death));
      if (value <= upper) {
        values.push_back(value);
      } else {
        ++clipped;
      }
    }
  }
  const std::size_t raw = values.size() + clipped;
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
  return {raw, clipped, values.size()};
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
  const SolverConfig reference{CandidateStrategy::sort_all,
                               ThresholdStrategy::binary,
                               DistanceStrategy::dense_aos,
                               AdjacencyStrategy::on_demand,
                               MatcherStrategy::kuhn,
                               VertexOrder::natural};
  const auto rebuild = [](ThresholdStrategy threshold) {
    return SolverConfig{CandidateStrategy::x_sweep_clipped, threshold,
                        DistanceStrategy::recompute_soa,
                        AdjacencyStrategy::x_sweep_csr,
                        MatcherStrategy::reusable_greedy_kuhn,
                        VertexOrder::degree_ascending};
  };
  const std::vector<Variant> variants{
      {"adaptive", {}},
      {"binary_rebuild", rebuild(ThresholdStrategy::binary)},
      {"gudhi_alpha_rebuild", rebuild(ThresholdStrategy::gudhi_alpha)},
      {"galloping_rebuild", rebuild(ThresholdStrategy::exponential)},
      {"quickselect_rebuild", rebuild(ThresholdStrategy::quickselect)},
      {"incremental_monotone", rebuild(ThresholdStrategy::incremental)},
      {"incremental_blocked", rebuild(ThresholdStrategy::incremental_blocked)},
      {"geometric_refinement", {CandidateStrategy::x_sweep_clipped,
       ThresholdStrategy::geometric_refinement, DistanceStrategy::recompute_soa,
       AdjacencyStrategy::on_demand, MatcherStrategy::geometric_hopcroft_karp,
       VertexOrder::natural}},
  };
  const std::vector<std::size_t> sizes{64, 128, 256, 512};
  const std::vector<Pattern> patterns{Pattern::uniform, Pattern::clustered,
                                      Pattern::near_diagonal, Pattern::threshold_shell};
  std::mt19937_64 data_generator(0xE10B0771EULL);
  std::mt19937_64 order_generator(0x5EA2C4ULL);
  double sink = 0.0;

  std::cout << "pattern,n,raw_candidates,clipped_candidates,unique_candidates,duplicate_ratio,"
               "variant,repetitions,rounds,median_us,p95_us,distance,decisions,refinement_rounds,"
               "incremental_groups_added,incremental_rollbacks,adjacency_checks,emitted_edges,"
               "augment_searches\n";
  for (Pattern pattern : patterns) {
    if (options.pattern != "all" && options.pattern != pattern_name(pattern)) {
      continue;
    }
    for (std::size_t size : sizes) {
      if (size < options.min_points || size > options.max_points) {
        continue;
      }
      const PreparedDiagram first(generate_diagram(data_generator, pattern, size, false));
      const PreparedDiagram second(generate_diagram(data_generator, pattern, size, true));
      const double expected = bottleneck::bottleneck_distance(first, second, reference);
      const CandidateShape shape = candidate_shape(first, second);
      const double duplicate_ratio = shape.raw == shape.clipped
                                         ? 0.0
                                         : 1.0 - static_cast<double>(shape.unique) /
                                                     static_cast<double>(shape.raw - shape.clipped);
      std::vector<std::vector<double>> samples(variants.size());
      std::vector<SolverStats> stats(variants.size());
      std::vector<std::size_t> order(variants.size());
      std::iota(order.begin(), order.end(), 0);
      for (const Variant& variant : variants) {
        const double value = bottleneck::bottleneck_distance(first, second, variant.config);
        if (value != expected) {
          std::cerr << "warmup mismatch: pattern=" << pattern_name(pattern)
                    << ", n=" << size << ", variant=" << variant.name << '\n';
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
            value = bottleneck::bottleneck_distance(first, second,
                                                    variants[variant_index].config,
                                                    &stats[variant_index]);
            sink += value;
          }
          const auto stop = Clock::now();
          if (value != expected) {
            std::cerr << "mismatch: pattern=" << pattern_name(pattern)
                      << ", n=" << size << ", variant="
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
        std::cout << pattern_name(pattern) << ',' << size << ',' << shape.raw << ','
                  << shape.clipped << ',' << shape.unique << ',' << std::fixed
                  << std::setprecision(6) << duplicate_ratio << ','
                  << variants[index].name << ',' << options.repetitions << ','
                  << options.rounds << ',' << std::setprecision(3)
                  << percentile(samples[index], 0.5) << ','
                  << percentile(samples[index], 0.95) << ',' << expected << ','
                  << static_cast<double>(current.threshold_decisions) / calls << ','
                  << static_cast<double>(current.refinement_rounds) / calls << ','
                  << static_cast<double>(current.incremental_groups_added) / calls << ','
                  << static_cast<double>(current.incremental_rollbacks) / calls << ','
                  << static_cast<double>(current.adjacency_checks) / calls << ','
                  << static_cast<double>(current.emitted_edges) / calls << ','
                  << static_cast<double>(current.augment_searches) / calls << '\n';
      }
    }
  }
  if (sink == -1.0) {
    std::cerr << "unreachable sink\n";
  }
  return 0;
}
