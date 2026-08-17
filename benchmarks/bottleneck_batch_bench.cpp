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
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
  std::size_t repetitions = 3;
  std::size_t rounds = 7;
  std::size_t min_points = 16;
  std::size_t max_points = 256;
  std::size_t min_batch = 1;
  std::size_t max_batch = 1024;
};

struct Variant {
  std::string_view name;
  std::vector<double> samples;
};

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

std::size_t parse_size(const char* value, std::string_view option) {
  const std::size_t parsed = static_cast<std::size_t>(std::stoull(value));
  if (parsed == 0) {
    throw std::invalid_argument(std::string(option) + " must be positive");
  }
  return parsed;
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string_view option(argv[index]);
    if (index + 1 >= argc) {
      throw std::invalid_argument("missing value for " + std::string(option));
    }
    const char* value = argv[++index];
    if (option == "--repetitions") {
      options.repetitions = parse_size(value, option);
    } else if (option == "--rounds") {
      options.rounds = parse_size(value, option);
    } else if (option == "--min-points") {
      options.min_points = parse_size(value, option);
    } else if (option == "--max-points") {
      options.max_points = parse_size(value, option);
    } else if (option == "--min-batch") {
      options.min_batch = parse_size(value, option);
    } else if (option == "--max-batch") {
      options.max_batch = parse_size(value, option);
    } else {
      throw std::invalid_argument("unknown option: " + std::string(option));
    }
  }
  if (options.min_points > options.max_points) {
    throw std::invalid_argument("--min-points must not exceed --max-points");
  }
  if (options.min_batch > options.max_batch) {
    throw std::invalid_argument("--min-batch must not exceed --max-batch");
  }
  return options;
}

double percentile(std::vector<double> values, double quantile) {
  std::sort(values.begin(), values.end());
  const std::size_t index = static_cast<std::size_t>(
      std::ceil(quantile * static_cast<double>(values.size())) - 1.0);
  return values[(std::min)(index, values.size() - 1)];
}

void require_equal(const std::vector<double>& actual, const std::vector<double>& expected,
                   std::string_view variant, std::size_t points, std::size_t batch_size) {
  if (actual.size() != expected.size()) {
    throw std::runtime_error("batch output size mismatch");
  }
  for (std::size_t index = 0; index < actual.size(); ++index) {
    if (actual[index] != expected[index]) {
      throw std::runtime_error("distance mismatch: variant=" + std::string(variant) +
                               ", points=" + std::to_string(points) +
                               ", batch=" + std::to_string(batch_size));
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    const std::vector<std::size_t> point_counts{16, 64, 256};
    const std::vector<std::size_t> batch_sizes{1, 4, 16, 64, 256, 1024};
    std::mt19937_64 data_generator(0xBA7C4ULL);
    std::mt19937_64 order_generator(0xE11BA7C4ULL);
    double sink = 0.0;

    std::cout << "points,batch_size,repetitions,rounds,query_prepare_median_us,"
                 "targets_prepare_median_us,variant,median_us,p95_us,per_pair_median_us,"
                 "speedup_vs_raw\n";
    for (std::size_t points : point_counts) {
      if (points < options.min_points || points > options.max_points) {
        continue;
      }
      for (std::size_t batch_size : batch_sizes) {
        if (batch_size < options.min_batch || batch_size > options.max_batch) {
          continue;
        }
        const bottleneck::Diagram query_diagram = random_diagram(data_generator, points);
        std::vector<bottleneck::Diagram> raw_targets;
        raw_targets.reserve(batch_size);
        for (std::size_t index = 0; index < batch_size; ++index) {
          raw_targets.push_back(random_diagram(data_generator, points));
        }
        const bottleneck::PreparedDiagram query(query_diagram);
        std::vector<bottleneck::PreparedDiagram> prepared_targets;
        prepared_targets.reserve(batch_size);
        for (const auto& target : raw_targets) {
          prepared_targets.emplace_back(target);
        }

        std::vector<double> expected(batch_size);
        bottleneck::bottleneck_distances(query, prepared_targets, expected);
        std::vector<double> warmup(batch_size);
        for (std::size_t index = 0; index < batch_size; ++index) {
          warmup[index] = bottleneck::bottleneck_distance(query_diagram, raw_targets[index]);
        }
        require_equal(warmup, expected, "raw_loop", points, batch_size);
        require_equal(bottleneck::bottleneck_distances(query, prepared_targets), expected,
                      "batch_allocated", points, batch_size);

        std::vector<double> query_prepare_samples;
        std::vector<double> targets_prepare_samples;
        std::vector<Variant> variants{{"raw_loop", {}},
                                      {"prepared_query", {}},
                                      {"prepared_loop", {}},
                                      {"batch_allocated", {}},
                                      {"batch_caller_buffer", {}}};
        std::vector<std::size_t> order(variants.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::vector<double> output(batch_size);

        for (std::size_t round = 0; round < options.rounds; ++round) {
          auto start = Clock::now();
          for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
            const bottleneck::PreparedDiagram prepared(query_diagram);
            sink += static_cast<double>(prepared.finite_points().size());
          }
          auto stop = Clock::now();
          query_prepare_samples.push_back(
              std::chrono::duration<double, std::micro>(stop - start).count() /
              static_cast<double>(options.repetitions));

          start = Clock::now();
          for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
            std::vector<bottleneck::PreparedDiagram> prepared;
            prepared.reserve(batch_size);
            for (const auto& target : raw_targets) {
              prepared.emplace_back(target);
            }
            sink += static_cast<double>(prepared.front().finite_points().size());
          }
          stop = Clock::now();
          targets_prepare_samples.push_back(
              std::chrono::duration<double, std::micro>(stop - start).count() /
              static_cast<double>(options.repetitions));

          std::shuffle(order.begin(), order.end(), order_generator);
          for (std::size_t variant_index : order) {
            start = Clock::now();
            for (std::size_t repetition = 0; repetition < options.repetitions; ++repetition) {
              switch (variant_index) {
                case 0:
                  for (std::size_t index = 0; index < batch_size; ++index) {
                    output[index] =
                        bottleneck::bottleneck_distance(query_diagram, raw_targets[index]);
                  }
                  break;
                case 1:
                  for (std::size_t index = 0; index < batch_size; ++index) {
                    output[index] = bottleneck::bottleneck_distance(
                        query, bottleneck::PreparedDiagram(raw_targets[index]));
                  }
                  break;
                case 2:
                  for (std::size_t index = 0; index < batch_size; ++index) {
                    output[index] = bottleneck::bottleneck_distance(query,
                                                                    prepared_targets[index]);
                  }
                  break;
                case 3:
                  output = bottleneck::bottleneck_distances(query, prepared_targets);
                  break;
                case 4:
                  bottleneck::bottleneck_distances(query, prepared_targets, output);
                  break;
                default:
                  throw std::runtime_error("unknown benchmark variant");
              }
              sink += output.front();
            }
            stop = Clock::now();
            require_equal(output, expected, variants[variant_index].name, points, batch_size);
            variants[variant_index].samples.push_back(
                std::chrono::duration<double, std::micro>(stop - start).count() /
                static_cast<double>(options.repetitions));
          }
        }

        const double query_prepare_median = percentile(query_prepare_samples, 0.5);
        const double targets_prepare_median = percentile(targets_prepare_samples, 0.5);
        const double raw_median = percentile(variants.front().samples, 0.5);
        for (const Variant& variant : variants) {
          const double median = percentile(variant.samples, 0.5);
          std::cout << points << ',' << batch_size << ',' << options.repetitions << ','
                    << options.rounds << ',' << std::fixed << std::setprecision(3)
                    << query_prepare_median << ',' << targets_prepare_median << ','
                    << variant.name << ',' << median << ','
                    << percentile(variant.samples, 0.95) << ','
                    << median / static_cast<double>(batch_size) << ','
                    << raw_median / median << '\n';
        }
      }
    }
    if (sink == -1.0) {
      std::cerr << "unreachable sink value\n";
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
