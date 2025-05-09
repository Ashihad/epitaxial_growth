#include "RngDriver.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#define DEBUG

/***************************************************************************
 *   generatory liczb pseudolosowych
 *
 *    UWAGA: generator rand() z biblioteki kompilatora C jest restartowany
 *    gdy wywolywana jest procedura PARDISO do rozwiazania ukladu
 *    (teraz uzywamy metody iteracyjnej CG)
 *
 ***************************************************************************/

RNGDriver::RNGDriver(const std::size_t& grid_x,
                     const std::size_t& diffusion_range)
    : m_grid_x_ref{grid_x},
      m_diffusion_range_ref{diffusion_range},
      rd{},
      mt{rd()} {}

double RNGDriver::gen_uniform() {
  static std::uniform_real_distribution<> uniform(0.0, 1.0);
  return uniform(mt);
}

std::size_t RNGDriver::gen_discrete_1_grid_x() {
  static std::uniform_int_distribution<std::size_t> one_to_grid_x(1,
                                                                  m_grid_x_ref);
  return one_to_grid_x(mt);
}

int RNGDriver::gen_discrete_plus_minus_diffusion_range() {
  static std::uniform_int_distribution<int> plus_minus_diffusion_range(
      -static_cast<int>(m_diffusion_range_ref),
      static_cast<int>(m_diffusion_range_ref));
  return plus_minus_diffusion_range(mt);
}

namespace LegacyRNG {

double gen_uniform() {
#ifdef DEBUG
  static unsigned fake_rand{};
  fake_rand = (1664525U * fake_rand + 1013904223U) % (RAND_MAX * 2U + 1U);
  return (static_cast<double>(fake_rand) /
          static_cast<double>(RAND_MAX * 2U + 1U));
#else
  return (static_cast<double>(rand()) / RAND_MAX);
#endif
}

// generate int from 1 to K INCLUSIVE
int gen_discrete_1_K(const std::size_t& K) {
  double uniform{gen_uniform()};
  return static_cast<int>(lround(floor(uniform * static_cast<double>(K)) + 1));
}

int gen_sign() {
  if (gen_uniform() < 0.5)
    return -1;
  else
    return 1;
}

int gen_discrete_1_K_multiply_sign(const std::size_t& K) {
  int discrete{gen_discrete_1_K(K)};
  int sign{gen_sign()};
  return static_cast<int>(discrete) * sign;
}

}  // namespace LegacyRNG
