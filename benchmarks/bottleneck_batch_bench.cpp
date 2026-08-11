#include <bottleneck/core.hpp>

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string_view>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

bottleneck::Diagram random_diagram(std::mt19937_64& generator, std::size_t size) {
  std::uniform_real_distribution<double> birth(-2.0, 2.0);
  std::uniform_real_distribution<double> persistence(0.01, 2.0);
  bottleneck::Diagram result;
  result.reserve(size);
  for (std::size_t index = 0; index < size; ++index) {
    const double value = birth(generator);
    result.push_back({value, value + persistence(generator)});
  }
  return result;
}

std::size_t repetitions(int argc, char** argv) {
  if (argc == 3 && std::string_view(argv[1]) == "--repetitions") {
    return static_cast<std::size_t>(std::stoul(argv[2]));
  }
  return 5;
}

}  // namespace

int main(int argc, char** argv) {
  const std::size_t repeat = repetitions(argc, argv);
  std::mt19937_64 generator(0xBA7C4ULL);
  std::cout << "points,batch_size,repetitions,raw_loop_us,prepared_loop_us,native_allocated_us,"
               "native_caller_buffer_us,raw_over_caller,prepared_over_caller\n";
  double sink = 0.0;
  for (std::size_t points : {std::size_t{16}, std::size_t{64}, std::size_t{256}}) {
    for (std::size_t batch_size : {std::size_t{20}, std::size_t{100}}) {
      const bottleneck::Diagram query_diagram = random_diagram(generator, points);
      std::vector<bottleneck::Diagram> raw_targets;
      std::vector<bottleneck::PreparedDiagram> prepared_targets;
      raw_targets.reserve(batch_size);
      prepared_targets.reserve(batch_size);
      for (std::size_t index = 0; index < batch_size; ++index) {
        raw_targets.push_back(random_diagram(generator, points));
        prepared_targets.emplace_back(raw_targets.back());
      }
      const bottleneck::PreparedDiagram query(query_diagram);
      std::vector<double> output(batch_size);

      const auto raw_start = Clock::now();
      for (std::size_t round = 0; round < repeat; ++round) {
        for (const auto& target : raw_targets) {
          sink += bottleneck::bottleneck_distance(query_diagram, target);
        }
      }
      const auto raw_stop = Clock::now();

      const auto prepared_start = Clock::now();
      for (std::size_t round = 0; round < repeat; ++round) {
        for (const auto& target : prepared_targets) {
          sink += bottleneck::bottleneck_distance(query, target);
        }
      }
      const auto prepared_stop = Clock::now();

      const auto allocated_start = Clock::now();
      for (std::size_t round = 0; round < repeat; ++round) {
        const auto values = bottleneck::bottleneck_distances(query, prepared_targets);
        sink += values.front();
      }
      const auto allocated_stop = Clock::now();

      const auto caller_start = Clock::now();
      for (std::size_t round = 0; round < repeat; ++round) {
        bottleneck::bottleneck_distances(query, prepared_targets, output);
        sink += output.front();
      }
      const auto caller_stop = Clock::now();

      const auto micros = [&](auto start, auto stop) {
        return std::chrono::duration<double, std::micro>(stop - start).count() /
               static_cast<double>(repeat);
      };
      const double raw_us = micros(raw_start, raw_stop);
      const double prepared_us = micros(prepared_start, prepared_stop);
      const double allocated_us = micros(allocated_start, allocated_stop);
      const double caller_us = micros(caller_start, caller_stop);
      std::cout << points << ',' << batch_size << ',' << repeat << ',' << std::fixed
                << std::setprecision(3) << raw_us << ',' << prepared_us << ',' << allocated_us
                << ',' << caller_us << ',' << raw_us / caller_us << ','
                << prepared_us / caller_us << '\n';
    }
  }
  if (sink == -1.0) {
    std::cerr << "unreachable sink value\n";
  }
  return 0;
}
