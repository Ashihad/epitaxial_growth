#include <gtest/gtest.h>
#include "gtest_common_define.hpp"
#include "linalg.hpp"

#include <vector>

class LinalgTest : public testing::Test {};

TEST_F(LinalgTest, scalar_product_basic) {
  // Primitive generator test, TODO: make something real
  double example_x_vector[] = {1, -2, 3, -4, 5, -6};
  double example_y_vector[] = {7, 1/8., 9, 1/10., 11, 1/12.};
  std::size_t vec_size = 6ul;

  double result = scalar_product_x_y(vec_size, example_x_vector, example_y_vector);
  EXPECT_DOUBLE_EQ(result, 87.85);
}
