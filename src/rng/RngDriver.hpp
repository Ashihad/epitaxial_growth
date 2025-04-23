#pragma once

#include <cstddef>
#include <random>

class RNGDriver {
 public:
  RNGDriver(const std::size_t& grid_x, const std::size_t& diffusion_range);
  double gen_uniform();
  std::size_t gen_discrete_1_grid_x();
  int gen_discrete_plus_minus_diffusion_range();

  const std::size_t& m_grid_x_ref;
  const std::size_t& m_diffusion_range_ref;
  std::random_device rd;
  // Mersenne Twister generator, 19937 variant
  // https://pl.wikipedia.org/wiki/Mersenne_Twister
  std::mt19937 mt;
};

namespace LegacyRNG {
double gen_uniform();
int gen_discrete_1_K(const std::size_t& K);
int gen_sign();
int gen_discrete_1_K_multiply_sign(const std::size_t& K);
void test_gen_discrete_sign();
}  // namespace LegacyRNG