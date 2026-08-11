#include <algorithm>
#include <cmath>
#include <cstddef>
#include <immintrin.h>

namespace bottleneck::detail {

void fill_cross_row_avx2(double first_birth, double first_death, const double* second_births,
                         const double* second_deaths, std::size_t size, double* output) noexcept {
  const __m256d birth = _mm256_set1_pd(first_birth);
  const __m256d death = _mm256_set1_pd(first_death);
  const __m256d sign = _mm256_set1_pd(-0.0);
  std::size_t index = 0;
  for (; index + 4 <= size; index += 4) {
    const __m256d birth_delta =
        _mm256_andnot_pd(sign, _mm256_sub_pd(birth, _mm256_loadu_pd(second_births + index)));
    const __m256d death_delta =
        _mm256_andnot_pd(sign, _mm256_sub_pd(death, _mm256_loadu_pd(second_deaths + index)));
    _mm256_storeu_pd(output + index, _mm256_max_pd(birth_delta, death_delta));
  }
  for (; index < size; ++index) {
    output[index] = (std::max)(std::fabs(first_birth - second_births[index]),
                               std::fabs(first_death - second_deaths[index]));
  }
}

}  // namespace bottleneck::detail
