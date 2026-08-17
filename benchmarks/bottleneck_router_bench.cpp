#include <bottleneck/core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
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

enum class Pattern {
  uniform,
  clustered,
  near_diagonal,
  duplicate_heavy,
  separated,
  threshold_shell
};

struct Options {
  std::size_t repetitions = 3;
  std::size_t rounds = 7;
  std::size_t max_points = 512;
  std::string_view pattern = "all";
};

struct Workload {
  Pattern pattern;
  std::size_t first_size;
  std::size_t second_size;
};

struct Variant {
  const char* name;
  SolverConfig config;
  std::size_t max_points;
};

struct Result {
  std::vector<double> samples;
  SolverStats stats;
};

const char* pattern_name(Pattern pattern) {
  switch (pattern) {
    case Pattern::uniform:
      return "uniform";
    case Pattern::clustered:
      return "clustered";
    case Pattern::near_diagonal:
      return "near_diagonal";
    case Pattern::duplicate_heavy:
      return "duplicate_heavy";
    case Pattern::separated:
      return "separated";
    case Pattern::threshold_shell:
      return "threshold_shell";
  }
  return "unknown";
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; index += 2) {
    if (index + 1 >= argc) {
      throw std::invalid_argument("missing benchmark option value");
    }
    const std::string_view name(argv[index]);
    if (name == "--repetitions") {
      options.repetitions = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--rounds") {
      options.rounds = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--max-points") {
      options.max_points = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--pattern") {
      options.pattern = argv[index + 1];
    } else {
      throw std::invalid_argument("unknown benchmark option");
    }
  }
  if (options.repetitions == 0 || options.rounds == 0 || options.max_points == 0) {
    throw std::invalid_argument("benchmark sizes must be positive");
  }
  return options;
}

Diagram generate_diagram(std::mt19937_64& generator, Pattern pattern, std::size_t size,
                         bool second) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 0.15);
  const double shell_radius = std::ldexp(1.0, -10);
  Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    double birth = 0.0;
    double persistence = 0.0;
    switch (pattern) {
      case Pattern::uniform:
        birth = -4.0 + 8.0 * unit(generator);
        persistence = 0.05 + 1.95 * unit(generator);
        break;
      case Pattern::clustered: {
        const double center = static_cast<double>(static_cast<int>(index % 4) - 2) * 1.5;
        birth = center + normal(generator) + (second ? 0.01 : 0.0);
        persistence = 0.4 + 0.8 * unit(generator);
        break;
      }
      case Pattern::near_diagonal:
        birth = -4.0 + 8.0 * unit(generator);
        persistence = 0.0001 + 0.04 * unit(generator);
        break;
      case Pattern::duplicate_heavy:
        birth = static_cast<double>(static_cast<int>(index % 8) - 4) * 0.5 +
                (second ? 0.01 : 0.0);
        persistence = 0.25 + 0.125 * static_cast<double>(index % 4);
        break;
      case Pattern::separated:
        birth = (second ? 10.0 : -12.0) + 2.0 * unit(generator);
        persistence = 0.05 + 1.95 * unit(generator);
        break;
      case Pattern::threshold_shell:
        birth = static_cast<double>(index) * shell_radius /
                static_cast<double>((std::max)(std::size_t{1}, size));
        persistence = 2.0 + (second ? shell_radius : 0.0) - birth;
        break;
    }
    result.push_back({birth, birth + persistence});
  }
  return result;
}

double percentile(std::vector<double> values, double quantile) {
  std::sort(values.begin(), values.end());
  const std::size_t index = static_cast<std::size_t>(
      std::ceil(quantile * static_cast<double>(values.size())) - 1.0);
  return values[(std::min)(index, values.size() - 1)];
}

double exact_cross_density(const PreparedDiagram& first, const PreparedDiagram& second,
                           double threshold) {
  if (first.finite_points().empty() || second.finite_points().empty()) {
    return 0.0;
  }
  std::size_t edges = 0;
  for (const auto& left : first.finite_points()) {
    for (const auto& right : second.finite_points()) {
      edges += static_cast<std::size_t>(
          (std::max)(std::fabs(left.birth - right.birth),
                     std::fabs(left.death - right.death)) <= threshold);
    }
  }
  return static_cast<double>(edges) /
         static_cast<double>(first.finite_points().size() * second.finite_points().size());
}

double duplicate_fraction(const PreparedDiagram& first, const PreparedDiagram& second) {
  const std::size_t raw = first.finite_points().size() + second.finite_points().size();
  if (raw == 0) {
    return 0.0;
  }
  const std::size_t unique = first.finite_duplicate_representatives().size() +
                             second.finite_duplicate_representatives().size();
  return 1.0 - static_cast<double>(unique) / static_cast<double>(raw);
}

double near_diagonal_fraction(const PreparedDiagram& first, const PreparedDiagram& second) {
  const std::size_t total = first.finite_points().size() + second.finite_points().size();
  const double upper = (std::max)(first.max_finite_diagonal_distance(),
                                  second.max_finite_diagonal_distance());
  if (total == 0 || upper <= 0.0) {
    return 0.0;
  }
  const double cutoff = upper * 0.05;
  std::size_t count = 0;
  for (double value : first.finite_diagonal_distances()) {
    count += static_cast<std::size_t>(value <= cutoff);
  }
  for (double value : second.finite_diagonal_distances()) {
    count += static_cast<std::size_t>(value <= cutoff);
  }
  return static_cast<double>(count) / static_cast<double>(total);
}

const char* adaptive_route(const SolverStats& stats) {
  if (stats.router_identity_shortcuts != 0) {
    return "identity";
  }
  if (stats.router_no_cross_shortcuts != 0) {
    return "no_cross";
  }
  if (stats.router_multiplicity_routes != 0) {
    return "multiplicity";
  }
  if (stats.router_mandatory_routes != 0) {
    return "mandatory_sparse";
  }
  if (stats.router_refinement_routes != 0) {
    return "refinement";
  }
  if (stats.router_quickselect_routes != 0) {
    return "quickselect";
  }
  return "explicit";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    const SolverConfig refinement{
        CandidateStrategy::x_sweep_clipped, ThresholdStrategy::geometric_refinement,
        DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
        MatcherStrategy::geometric_hopcroft_karp, VertexOrder::natural};
    const SolverConfig quickselect{
        CandidateStrategy::x_sweep_clipped, ThresholdStrategy::quickselect,
        DistanceStrategy::recompute_soa, AdjacencyStrategy::x_sweep_csr,
        MatcherStrategy::reusable_greedy_kuhn, VertexOrder::degree_ascending};
    const std::vector<Variant> variants{
        {"adaptive", {}, 512},
        {"refinement", refinement, 512},
        {"quickselect", quickselect, 512},
        {"geometric_candidates",
         {CandidateStrategy::x_sweep_clipped, ThresholdStrategy::quickselect,
          DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
          MatcherStrategy::geometric_hopcroft_karp, VertexOrder::natural},
         512},
        {"multiplicity",
         {CandidateStrategy::x_sweep_clipped, ThresholdStrategy::quickselect,
          DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
          MatcherStrategy::multiplicity_flow, VertexOrder::natural},
         256},
        {"mandatory_sparse",
         {CandidateStrategy::x_sweep_clipped, ThresholdStrategy::quickselect,
          DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
          MatcherStrategy::mandatory_sparse_flow, VertexOrder::natural},
         512}};
    const std::vector<Pattern> patterns{Pattern::uniform, Pattern::clustered,
                                        Pattern::near_diagonal, Pattern::duplicate_heavy,
                                        Pattern::separated, Pattern::threshold_shell};
    const std::vector<std::size_t> sizes{64, 128, 256, 512};
    std::vector<Workload> workloads;
    for (Pattern pattern : patterns) {
      for (std::size_t size : sizes) {
        workloads.push_back({pattern, size, size});
      }
    }
    for (Pattern pattern : {Pattern::uniform, Pattern::clustered, Pattern::near_diagonal}) {
      workloads.push_back({pattern, 32, 256});
      workloads.push_back({pattern, 256, 32});
      workloads.push_back({pattern, 64, 512});
      workloads.push_back({pattern, 512, 64});
    }

    std::mt19937_64 data_generator(0xE12ADA771ULL);
    std::mt19937_64 order_generator(0xE12B077EULL);
    double sink = 0.0;
    std::cout << "pattern,n,m,duplicate_fraction,near_diagonal_fraction,exact_density,"
                 "variant,route,repetitions,rounds,median_us,p95_us,regret_to_best,"
                 "router_x_window_pairs,decisions,refinement_rounds,emitted_edges,"
                 "augment_searches\n";
    for (const Workload& workload : workloads) {
      if ((std::max)(workload.first_size, workload.second_size) > options.max_points ||
          (options.pattern != "all" && options.pattern != pattern_name(workload.pattern))) {
        continue;
      }
      const PreparedDiagram first(
          generate_diagram(data_generator, workload.pattern, workload.first_size, false));
      const PreparedDiagram second(
          generate_diagram(data_generator, workload.pattern, workload.second_size, true));
      const double expected = bottleneck::bottleneck_distance(first, second, refinement);
      const double density = exact_cross_density(first, second, expected);
      std::vector<Result> results(variants.size());
      std::vector<std::size_t> enabled;
      for (std::size_t index = 0; index < variants.size(); ++index) {
        if ((std::max)(workload.first_size, workload.second_size) <= variants[index].max_points) {
          const double value =
              bottleneck::bottleneck_distance(first, second, variants[index].config);
          if (value != expected) {
            throw std::runtime_error("warmup mismatch: " +
                                     std::string(pattern_name(workload.pattern)) + "/" +
                                     variants[index].name);
          }
          enabled.push_back(index);
        }
      }
      for (std::size_t round = 0; round < options.rounds; ++round) {
        std::shuffle(enabled.begin(), enabled.end(), order_generator);
        for (std::size_t index : enabled) {
          const auto start = Clock::now();
          double value = 0.0;
          for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
            value = bottleneck::bottleneck_distance(first, second, variants[index].config,
                                                    &results[index].stats);
            sink += value;
          }
          const auto stop = Clock::now();
          if (value != expected) {
            throw std::runtime_error("distance mismatch: " +
                                     std::string(pattern_name(workload.pattern)) + "/" +
                                     variants[index].name);
          }
          results[index].samples.push_back(
              std::chrono::duration<double, std::micro>(stop - start).count() /
              static_cast<double>(options.repetitions));
        }
      }
      double best = std::numeric_limits<double>::infinity();
      for (std::size_t index : enabled) {
        best = (std::min)(best, percentile(results[index].samples, 0.5));
      }
      const double calls = static_cast<double>(options.repetitions * options.rounds);
      for (std::size_t index : enabled) {
        const double median = percentile(results[index].samples, 0.5);
        const SolverStats& stats = results[index].stats;
        std::cout << pattern_name(workload.pattern) << ',' << workload.first_size << ','
                  << workload.second_size << ',' << std::fixed << std::setprecision(6)
                  << duplicate_fraction(first, second) << ','
                  << near_diagonal_fraction(first, second) << ',' << density << ','
                  << variants[index].name << ','
                  << (index == 0 ? adaptive_route(stats) : "explicit") << ','
                  << options.repetitions << ',' << options.rounds << ','
                  << std::setprecision(3) << median << ','
                  << percentile(results[index].samples, 0.95) << ',' << median / best << ','
                  << static_cast<double>(stats.router_x_window_pairs) / calls << ','
                  << static_cast<double>(stats.threshold_decisions) / calls << ','
                  << static_cast<double>(stats.refinement_rounds) / calls << ','
                  << static_cast<double>(stats.emitted_edges) / calls << ','
                  << static_cast<double>(stats.augment_searches) / calls << '\n';
      }
    }
    if (sink == -1.0) {
      std::cerr << "unreachable sink\n";
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
