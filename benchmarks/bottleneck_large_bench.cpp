#include <bottleneck/core.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
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

enum class Pattern { uniform, near_diagonal, clustered, repeated, separated };

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
        birth = -4.0 + 8.0 * unit(generator) + (second ? 0.002 : 0.0);
        persistence = 0.01 + 2.99 * unit(generator);
        break;
      case Pattern::near_diagonal:
        birth = -4.0 + 8.0 * unit(generator);
        persistence = 0.00001 + 0.01999 * unit(generator);
        break;
      case Pattern::clustered: {
        const double center = static_cast<double>(static_cast<int>(index % 8) - 4) * 0.5;
        birth = center + 0.02 * signed_unit(generator) + (second ? 0.003 : 0.0);
        persistence = 0.4 + 0.05 * static_cast<double>(index % 5) +
                      0.01 * signed_unit(generator);
        break;
      }
      case Pattern::repeated:
        birth = static_cast<double>(static_cast<int>(index % 7) - 3) * 0.35;
        persistence = 0.2 + 0.1 * static_cast<double>(index % 6);
        break;
      case Pattern::separated:
        birth = (second ? 10.0 : -12.0) + 2.0 * unit(generator);
        persistence = 0.05 + 2.95 * unit(generator);
        break;
    }
    result.push_back({birth, birth + persistence});
  }
  return result;
}

double cross_density(const PreparedDiagram& first, const PreparedDiagram& second,
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

struct Options {
  std::size_t repetitions = 3;
  std::size_t max_points = 4096;
};

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index + 1 < argc; index += 2) {
    const std::string_view name(argv[index]);
    if (name == "--repetitions") {
      options.repetitions = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else if (name == "--max-points") {
      options.max_points = static_cast<std::size_t>(std::stoul(argv[index + 1]));
    } else {
      throw std::invalid_argument("unknown benchmark option");
    }
  }
  return options;
}

}  // namespace

int main(int argc, char** argv) {
  const Options options = parse_options(argc, argv);
  const std::size_t repeat = options.repetitions;
  const SolverConfig refined{
      CandidateStrategy::sort_unique_clipped, ThresholdStrategy::geometric_refinement,
      DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
      MatcherStrategy::geometric_hopcroft_karp, VertexOrder::natural};
  const std::vector<Variant> variants{
      {"geometric_refinement", refined, 4096},
      {"geometric_candidates",
       {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
        DistanceStrategy::recompute_soa, AdjacencyStrategy::on_demand,
       MatcherStrategy::geometric_hopcroft_karp, VertexOrder::natural},
       2048},
      {"adaptive_dispatcher", {}, 4096},
      {"legacy_adaptive",
       {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
        DistanceStrategy::dense_aos, AdjacencyStrategy::adaptive,
        MatcherStrategy::adaptive, VertexOrder::degree_ascending},
       2048},
      {"legacy_mandatory_flow",
       {CandidateStrategy::sort_unique_clipped, ThresholdStrategy::quickselect,
        DistanceStrategy::dense_soa_avx2, AdjacencyStrategy::adaptive,
        MatcherStrategy::mandatory_flow, VertexOrder::natural},
       1024},
  };
  const std::vector<Workload> workloads{
      {Pattern::uniform, 32, 32},            {Pattern::uniform, 64, 64},
      {Pattern::uniform, 96, 96},            {Pattern::uniform, 128, 128},
      {Pattern::uniform, 192, 192},
      {Pattern::uniform, 256, 256},          {Pattern::uniform, 384, 384},
      {Pattern::uniform, 512, 512},          {Pattern::uniform, 1024, 1024},
      {Pattern::uniform, 2048, 2048},        {Pattern::uniform, 4096, 4096},
      {Pattern::near_diagonal, 64, 64},      {Pattern::near_diagonal, 128, 128},
      {Pattern::near_diagonal, 256, 256},    {Pattern::near_diagonal, 384, 384},
      {Pattern::near_diagonal, 512, 512},    {Pattern::near_diagonal, 2048, 2048},
      {Pattern::clustered, 64, 64},          {Pattern::clustered, 128, 128},
      {Pattern::clustered, 256, 256},        {Pattern::clustered, 384, 384},
      {Pattern::clustered, 512, 512},        {Pattern::clustered, 2048, 2048},
      {Pattern::repeated, 512, 512},         {Pattern::repeated, 2048, 2048},
      {Pattern::separated, 64, 64},          {Pattern::separated, 128, 128},
      {Pattern::separated, 256, 256},        {Pattern::separated, 384, 384},
      {Pattern::separated, 512, 512},
      {Pattern::separated, 2048, 2048},
      {Pattern::uniform, 512, 2048},         {Pattern::uniform, 2048, 512},
      {Pattern::uniform, 8, 128},            {Pattern::uniform, 128, 8},
      {Pattern::uniform, 16, 128},           {Pattern::uniform, 128, 16},
      {Pattern::uniform, 32, 256},           {Pattern::uniform, 256, 32},
      {Pattern::clustered, 16, 128},         {Pattern::clustered, 128, 16},
      {Pattern::clustered, 32, 256},         {Pattern::clustered, 256, 32},
      {Pattern::clustered, 512, 4096},       {Pattern::clustered, 4096, 512},
  };

  std::cout << "pattern,n,m,variant,repetitions,mean_us,relative_to_refined,distance,density,"
               "decisions,retained_candidates,clipped_candidates,geometric_queries,kd_nodes,"
               "bfs_phases,augment_searches\n";
  std::mt19937_64 generator(0x1A26EULL);
  double sink = 0.0;
  for (const Workload& workload : workloads) {
    if ((std::max)(workload.first_size, workload.second_size) > options.max_points) {
      continue;
    }
    std::vector<std::pair<PreparedDiagram, PreparedDiagram>> pairs;
    pairs.reserve(repeat);
    for (std::size_t index = 0; index < repeat; ++index) {
      pairs.emplace_back(generate_diagram(generator, workload.pattern, workload.first_size, false),
                         generate_diagram(generator, workload.pattern, workload.second_size, true));
    }
    std::vector<double> expected;
    expected.reserve(repeat);
    for (const auto& pair : pairs) {
      expected.push_back(bottleneck::bottleneck_distance(pair.first, pair.second, refined));
    }
    double density = 0.0;
    for (std::size_t index = 0; index < repeat; ++index) {
      density += cross_density(pairs[index].first, pairs[index].second, expected[index]);
    }
    density /= static_cast<double>(repeat);

    double refined_us = 0.0;
    for (std::size_t variant_index = 0; variant_index < variants.size(); ++variant_index) {
      const Variant& variant = variants[variant_index];
      if ((std::max)(workload.first_size, workload.second_size) > variant.max_points) {
        continue;
      }
      SolverStats stats;
      const auto start = Clock::now();
      for (std::size_t index = 0; index < repeat; ++index) {
        const double value = bottleneck::bottleneck_distance(
            pairs[index].first, pairs[index].second, variant.config, &stats);
        if (value != expected[index]) {
          std::cerr << "mismatch: " << pattern_name(workload.pattern) << ' '
                    << workload.first_size << 'x' << workload.second_size << ' '
                    << variant.name << '\n';
          return 2;
        }
        sink += value;
      }
      const auto stop = Clock::now();
      const double mean_us =
          std::chrono::duration<double, std::micro>(stop - start).count() /
          static_cast<double>(repeat);
      if (variant_index == 0) {
        refined_us = mean_us;
      }
      std::cout << pattern_name(workload.pattern) << ',' << workload.first_size << ','
                << workload.second_size << ',' << variant.name << ',' << repeat << ','
                << std::fixed << std::setprecision(3) << mean_us << ','
                << mean_us / refined_us << ',' << expected.front() << ',' << density << ','
                << stats.threshold_decisions << ',' << stats.retained_candidates << ','
                << stats.clipped_candidates << ',' << stats.geometric_queries << ','
                << stats.kd_nodes_visited << ',' << stats.bfs_phases << ','
                << stats.augment_searches << '\n';
    }
  }
  if (sink == -1.0) {
    std::cerr << "unreachable sink\n";
  }
  return 0;
}
