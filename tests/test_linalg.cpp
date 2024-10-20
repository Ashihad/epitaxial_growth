#include <gtest/gtest.h>
#include "gtest_common_define.hpp"
#include "linalg.hpp"

#include <vector>

class LinalgTest : public testing::Test {};

TEST_F(LinalgTest, scalar_product_basic) {
  // Primitive generator test, TODO: make something real
  double example_x_vector[] = {1, -2, 3, -4, 5, -6};
  double example_y_vector[] = {7, 1 / 8., 9, 1 / 10., 11, 1 / 12.};
  std::size_t vec_size = 6ul;

  double result =
      scalar_product_x_y(vec_size, example_x_vector, example_y_vector);
  EXPECT_DOUBLE_EQ(result, 87.85);
}

TEST_F(LinalgTest, compute_sparse_Ax_y_basic) {
  std::size_t n = 3;
  // CSR format matrix:
  // [1 0 0]
  // [0 2 3]
  // [4 0 0]
  double csr_val[] = {1, 2, 3, 4};
  int csr_row[] = {0, 1, 3, 4};
  int csr_column[] = {0, 1, 2, 0};

  double x[] = {1, 2, 3};
  double y[] = {0, 0, 0};

  double expected_y[] = {1, 13, 4};

  // Call the function
  compute_sparse_Ax_y(n, csr_val, csr_row, csr_column, x, y);

  // Check if the result matches the expected result
  for (std::size_t i = 0; i < n; i++) {
    EXPECT_DOUBLE_EQ(y[i], expected_y[i]);
  }
}

TEST(LinalgTest, compute_sparse_Ax_y_empty_matrix) {
  // There should be no access when matrix is empty (nullptr access -> segfault)
  std::size_t n = 0;

  // Empty matrix (all arrays empty)
  double* csr_val{nullptr};
  int* csr_row{nullptr};
  int* csr_column{nullptr};

  double* x{nullptr};
  double* y{nullptr};

  // Expected output: y = empty
  double* expected_y{nullptr};

  // Call the function
  compute_sparse_Ax_y(n, csr_val, csr_row, csr_column, x, y);

  // There should be no result since matrix is empty
  EXPECT_EQ(y, expected_y);
}

TEST(LinalgTest, LargeSparseMatrixTest) {
  std::size_t n = 10000;  // 10,000 rows

  // Allocate arrays for CSR representation
  double* csr_val = static_cast<double*>(
      malloc(sizeof(double) * (n * 2)));  // Overestimate number of non-zeros
  int* csr_row = static_cast<int*>(malloc(sizeof(int) * (n + 1)));
  int* csr_column = static_cast<int*>(
      malloc(sizeof(int) * (n * 2)));  // Overestimate number of non-zeros

  // Input vector x
  double* x = static_cast<double*>(malloc(sizeof(double) * n));
  double* y =
      static_cast<double*>(malloc(sizeof(double) * n));  // Output vector y

  // Initialize x to all ones
  for (std::size_t i = 0; i < n; i++) {
    x[i] = 1.0;
    y[i] = 0.0;  // Initialize y to 0
  }

  // Fill CSR arrays with a sparsity pattern (diagonal + off-diagonal elements)
  int nnz = 0;     // Non-zero count
  csr_row[0] = 0;  // First row start index
  for (std::size_t i = 0; i < n; i++) {
    if (i % 2 == 0) {
      // Diagonal element
      csr_val[nnz] = 2.0;                     // Non-zero value
      csr_column[nnz] = static_cast<int>(i);  // Column index
      nnz++;

      // Off-diagonal element
      if (i + 1 < n) {
        csr_val[nnz] = 1.0;                         // Non-zero value
        csr_column[nnz] = static_cast<int>(i + 1);  // Column index
        nnz++;
      }
    } else {
      // Single non-zero element in odd rows
      csr_val[nnz] = 3.0;
      csr_column[nnz] = static_cast<int>(i);
      nnz++;
    }
    csr_row[i + 1] = nnz;  // Update row start/end index
  }

  // Expected output: since x is all 1's, y should reflect the sum of non-zero
  // elements in each row
  double* expected_y = static_cast<double*>(malloc(sizeof(double) * n));
  for (std::size_t i = 0; i < n; i++) {
    if (i % 2 == 0) {
      expected_y[i] = 2.0 + (i + 1 < n ? 1.0 : 0.0);  // Diagonal + off-diagonal
    } else {
      expected_y[i] = 3.0;
    }
  }

  // Call the function
  compute_sparse_Ax_y(n, csr_val, csr_row, csr_column, x, y);

  // Check if the result matches the expected result
  for (std::size_t i = 0; i < n; i++) {
    EXPECT_DOUBLE_EQ(y[i], expected_y[i]);
  }

  // Free all allocated memory
  free(csr_val);
  free(csr_row);
  free(csr_column);
  free(x);
  free(y);
  free(expected_y);
}
