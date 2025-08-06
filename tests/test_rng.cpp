#include <gtest/gtest.h>
#include "RngDriver.hpp"
#include "gtest_common_define.hpp"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

class RngTest : public testing::Test {
 protected:
  RngTest() { std::srand(static_cast<unsigned int>(std::time(0))); }
  double uniform_eps{1e-1};
};

TEST_F(RngTest, uniform_generator_basic) {
  // Primitive generator test, TODO: make something real
  std::size_t grid_x{0};
  std::size_t diff_range{0};
  RNGDriver driver(grid_x, diff_range);

  double sum{};
  std::size_t max_iter{100000};
  for (std::size_t iter = 0; iter < max_iter; iter++) {
    double res{driver.gen_uniform()};
    EXPECT_LE(res, 1);
    EXPECT_GE(res, 0);
    sum += res;
  }
  bool is_ok = (sum / static_cast<double>(max_iter) <= 0.5 + uniform_eps &&
                sum / static_cast<double>(max_iter) >= 0.5 - uniform_eps);
  std::cout << COUT_GTEST << ANSI_TXT_DFT
            << "sum/max_iter=" << sum / static_cast<double>(max_iter)
            << ", expected: ~" << 1. / 2 << '\n';
  ASSERT_TRUE(is_ok);
}

TEST_F(RngTest, gen_discrete_1_grid_x_basic) {
  // Primitive generator test, TODO: make something real
  std::size_t grid_x{9};
  std::size_t diff_range{0};
  RNGDriver driver(grid_x, diff_range);

  std::size_t sum{};
  std::size_t max_iter{100000};
  for (std::size_t iter = 0; iter < max_iter; iter++) {
    std::size_t res{driver.gen_discrete_1_grid_x()};
    EXPECT_LE(res, grid_x);
    EXPECT_GE(res, 1);
    sum += res;
  }
  bool is_ok =
      (static_cast<double>(sum) / static_cast<double>(max_iter) <=
           std::floor(static_cast<double>(grid_x) / 2.) + 1 + uniform_eps &&
       static_cast<double>(sum) / static_cast<double>(max_iter) >=
           std::floor(static_cast<double>(grid_x) / 2.) + 1 - uniform_eps);
  std::cout << COUT_GTEST << ANSI_TXT_DFT << "sum/max_iter="
            << static_cast<double>(sum) / static_cast<double>(max_iter)
            << ", expected: ~" << grid_x / 2 + 1 << '\n';
  ASSERT_TRUE(is_ok);
}

TEST_F(RngTest, gen_discrete_plus_minus_diffusion_range_basic) {
  // Primitive generator test, TODO: make something real
  std::size_t grid_x{0};
  std::size_t diff_range{9};
  RNGDriver driver(grid_x, diff_range);

  int sum{};
  std::size_t max_iter{50000};
  for (std::size_t iter = 0; iter < max_iter; iter++) {
    int res{driver.gen_discrete_plus_minus_diffusion_range()};
    EXPECT_LE(res, static_cast<int>(diff_range));
    EXPECT_GE(res, -static_cast<int>(diff_range));
    sum += res;
  }
  bool is_ok = (sum / static_cast<double>(max_iter) <= 0.0 + uniform_eps &&
                sum / static_cast<double>(max_iter) >= 0.0 - uniform_eps);
  std::cout << COUT_GTEST << ANSI_TXT_DFT
            << "sum/max_iter=" << sum / static_cast<double>(max_iter)
            << ", expected: ~0" << '\n';
  ASSERT_TRUE(is_ok);
}
