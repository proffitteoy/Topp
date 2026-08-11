#include <bottleneck/wasserstein.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
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

bottleneck::Diagram random_diagram(std::mt19937_64& generator, std::size_t size) {
  std::uniform_real_distribution<double> birth(-4.0, 4.0);
  std::uniform_real_distribution<double> persistence(0.01, 3.0);
  bottleneck::Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    const double value = birth(generator);
    result.push_back({value, value + persistence(generator)});
  }
  return result;
}

struct Options {
  std::size_t repetitions = 10;
  std::size_t rounds = 7;
};

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index + 1 < argc; index += 2) {
    const std::string_view name(argv[index]);
    if (name == "--repetitions") {
      options.repetitions = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else if (name == "--rounds") {
      options.rounds = static_cast<std::size_t>(std::stoull(argv[index + 1]));
    } else {
      throw std::invalid_argument("unknown benchmark option");
    }
  }
  if (argc % 2 == 0 || options.repetitions == 0 || options.rounds == 0) {
    throw std::invalid_argument("benchmark options require positive values");
  }
  return options;
}

double percentile(std::vector<double> samples, double quantile) {
  std::sort(samples.begin(), samples.end());
  const std::size_t index = static_cast<std::size_t>(
      std::ceil(quantile * static_cast<double>(samples.size())) - 1.0);
  return samples[(std::min)(index, samples.size() - 1)];
}

}  // namespace

int main(int argc, char** argv) {
  const Options options = parse_options(argc, argv);
  std::mt19937_64 generator(0xBA7C4ULL);
  constexpr std::size_t diagram_size = 64;
  constexpr std::size_t target_count = 128;
  const bottleneck::Diagram raw_query = random_diagram(generator, diagram_size);
  std::vector<bottleneck::Diagram> raw_targets;
  std::vector<bottleneck::PreparedDiagram> targets;
  raw_targets.reserve(target_count);
  targets.reserve(target_count);
  for (std::size_t index = 0; index < target_count; ++index) {
    raw_targets.push_back(random_diagram(generator, diagram_size));
    targets.emplace_back(raw_targets.back());
  }
  const bottleneck::PreparedDiagram query(raw_query);
  const std::vector<double> reference = bottleneck::wasserstein_distances(query, targets);
  std::vector<double> output(target_count);
  bottleneck::WassersteinWorkspace workspace;

  struct Mode {
    std::string_view name;
    std::function<void()> run;
  };
  const std::array<Mode, 5> modes{{
      {"one_shot",
       [&] {
         for (std::size_t index = 0; index < target_count; ++index) {
           output[index] =
               bottleneck::wasserstein_distance(raw_query, raw_targets[index]);
         }
       }},
      {"prepared_loop",
       [&] {
         for (std::size_t index = 0; index < target_count; ++index) {
           output[index] = bottleneck::wasserstein_distance(query, targets[index]);
         }
       }},
      {"native_caller_buffer",
       [&] { bottleneck::wasserstein_distances(query, targets, output); }},
      {"native_workspace",
       [&] {
         bottleneck::wasserstein_distances(query, targets, output, workspace);
       }},
      {"native_allocated",
       [&] { output = bottleneck::wasserstein_distances(query, targets); }},
  }};

  std::array<std::vector<double>, 5> timings;
  std::array<std::size_t, 5> order{};
  std::iota(order.begin(), order.end(), std::size_t{0});
  double sink = 0.0;
  for (const Mode& mode : modes) {
    mode.run();
  }
  for (std::size_t round = 0; round < options.rounds; ++round) {
    std::shuffle(order.begin(), order.end(), generator);
    for (std::size_t mode_index : order) {
      const auto start = Clock::now();
      for (std::size_t repetition = 0; repetition < options.repetitions;
           ++repetition) {
        modes[mode_index].run();
        sink += output[repetition % output.size()];
      }
      const auto stop = Clock::now();
      timings[mode_index].push_back(
          std::chrono::duration<double, std::milli>(stop - start).count());
      for (std::size_t index = 0; index < target_count; ++index) {
        if (std::fabs(output[index] - reference[index]) > 2e-12) {
          std::cerr << modes[mode_index].name << " result mismatch at " << index
                    << '\n';
          return 2;
        }
      }
    }
  }

  std::cout << "mode,repetitions,rounds,target_count,median_ms,p95_ms,"
               "median_per_distance_us\n";
  for (std::size_t mode_index = 0; mode_index < modes.size(); ++mode_index) {
    const double median_ms = percentile(timings[mode_index], 0.5);
    std::cout << modes[mode_index].name << ',' << options.repetitions << ','
              << options.rounds << ',' << target_count << ',' << std::fixed
              << std::setprecision(3) << median_ms << ','
              << percentile(timings[mode_index], 0.95) << ','
              << median_ms * 1000.0 /
                     static_cast<double>(options.repetitions * target_count)
              << '\n';
  }
  if (sink < 0.0) {
    std::cerr << "unreachable sink value\n";
  }
  return 0;
}
