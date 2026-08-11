#include <algorithm>
#include <cmath>
#include <cstddef>
#include <immintrin.h>

namespace bottleneck::detail {

void fill_w1_savings_avx2(double first_u, double first_v, const double* second_u,
                          const double* second_v, std::size_t size,
                          double* output) noexcept {
  const __m256d u = _mm256_set1_pd(first_u);
  const __m256d v = _mm256_set1_pd(first_v);
  const __m256d two = _mm256_set1_pd(2.0);
  const __m256d sign = _mm256_set1_pd(-0.0);
  std::size_t index = 0;
  for (; index + 4 <= size; index += 4) {
    const __m256d other_u = _mm256_loadu_pd(second_u + index);
    const __m256d other_v = _mm256_loadu_pd(second_v + index);
    const __m256d delta = _mm256_andnot_pd(sign, _mm256_sub_pd(u, other_u));
    const __m256d minimum = _mm256_min_pd(v, other_v);
    _mm256_storeu_pd(output + index,
                     _mm256_sub_pd(_mm256_mul_pd(two, minimum), delta));
  }
  for (; index < size; ++index) {
    output[index] = 2.0 * (std::min)(first_v, second_v[index]) -
                    std::fabs(first_u - second_u[index]);
  }
}

void fill_w2_savings_avx2(double first_u, double first_v, const double* second_u,
                          const double* second_v, std::size_t size,
                          double* output) noexcept {
  const __m256d u = _mm256_set1_pd(first_u);
  const __m256d four_v = _mm256_set1_pd(4.0 * first_v);
  const __m256d two = _mm256_set1_pd(2.0);
  std::size_t index = 0;
  for (; index + 4 <= size; index += 4) {
    const __m256d other_u = _mm256_loadu_pd(second_u + index);
    const __m256d other_v = _mm256_loadu_pd(second_v + index);
    const __m256d delta = _mm256_sub_pd(u, other_u);
    const __m256d persistence = _mm256_mul_pd(four_v, other_v);
    const __m256d squared = _mm256_mul_pd(two, _mm256_mul_pd(delta, delta));
    _mm256_storeu_pd(output + index, _mm256_sub_pd(persistence, squared));
  }
  for (; index < size; ++index) {
    const double delta = first_u - second_u[index];
    output[index] = 4.0 * first_v * second_v[index] - 2.0 * delta * delta;
  }
}

}  // namespace bottleneck::detail
