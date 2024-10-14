#include <gtest/gtest.h>
#include "gtest_common_define.hpp"
#include "rng.hpp"

#include <cstdlib>
#include <ctime>
#include <iostream>

class RngTest : public testing::Test {
 protected:
  RngTest() {
     std::srand(static_cast<unsigned int>(std::time(0)));
  }
  double uniform_eps {1e-1};
};

TEST_F(RngTest, uniform_generator_basic) {
  // Primitive generator test, TODO: make something real
  double sum {};
  std::size_t max_iter {10000};
  for(std::size_t iter = 0; iter < max_iter; iter++)
  {
    sum += gen_uniform();
  }
  bool is_ok =  (sum/static_cast<double>(max_iter) <= 0.5+uniform_eps && 
                 sum/static_cast<double>(max_iter) >= 0.5-uniform_eps);
  std::cout << COUT_GTEST << ANSI_TXT_DFT <<  "sum/max_iter=" << sum/static_cast<double>(max_iter) << '\n';
  EXPECT_EQ(is_ok, true);
}

TEST_F(RngTest, gen_discrete_1_K_basic) {
  // Primitive generator test, TODO: make something real
  int k {9};

  int sum {};
  std::size_t max_iter {10000};
  for(std::size_t iter = 0; iter < max_iter; iter++)
  {
    sum += gen_discrete_1_K(k);
  }
  bool is_ok =  (sum/static_cast<double>(max_iter) <= k/2+1+uniform_eps && 
                 sum/static_cast<double>(max_iter) >= k/2+1-uniform_eps);
  std::cout << COUT_GTEST << ANSI_TXT_DFT <<  "sum/max_iter=" << sum/static_cast<double>(max_iter) << '\n';
  EXPECT_EQ(is_ok, true);
}

TEST_F(RngTest, gen_sign) {
  // Primitive generator test, TODO: make something real
  int sum {};
  std::size_t max_iter {10000};
  for(std::size_t iter = 0; iter < max_iter; iter++)
  {
    sum += gen_sign();
  }
  bool is_ok =  (sum/static_cast<double>(max_iter) <= 0.0+uniform_eps && 
                 sum/static_cast<double>(max_iter) >= 0.0-uniform_eps);
  std::cout << COUT_GTEST << ANSI_TXT_DFT <<  "sum/max_iter=" << sum/static_cast<double>(max_iter) << '\n';
  EXPECT_EQ(is_ok, true);
}
