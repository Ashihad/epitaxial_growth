#include "GSDriver.hpp"

#include <cmath>
#include <cstring>
#include <iostream>

GSDriver::GSDriver(const ConfigMathDriver& conf_md,
                   const ConfigPhysics& conf_ph)
    : MathDriver{conf_md, conf_ph} {}

void GSDriver::solve_linear_system(const CSRMatrix& A,
                                   const FastVector<double>& b,
                                   FastVector<double>& x) {
  std::size_t n{A.n_rows};

  // if b is zero vector, return trivial solution x={0}
  double basically_zero{1.0E-10};
  double b_dot_b = dot(b, b);
  if (b_dot_b < std::abs(basically_zero)) {
    std::fill_n(x.get(), n, 0.);
    m_last_iterations_no = 0;
    m_last_tolerance = 0.;
    return;
  }

  FastVector<double> x_old(n);
  std::fill_n(x_old.get(), n, 0.);

  double error;
  for (std::size_t iter{}; iter < m_max_iterations; iter++) {
    std::memcpy(x_old.get(), x.get(), x.size * sizeof(double));
    for (std::size_t i{}; i < n; i++) {
      double sum{};
      double diag{};
      for (std::size_t j{A.row_index[i]}; j < A.row_index[i + 1]; j++) {
        std::size_t col{A.col_index[j]};
        double val{A.value[j]};
        if (col == i)
          diag = val;
        else
          sum += val * x[col];
      }
      if (std::abs(diag) < 1e-12) {
        std::cerr << "Zero in diag at row " << i << ", cannot divide by zero\n";
        std::abort();
      }
      x[i] = (b[i] - sum) / diag;
    }
    error = 0.;
    for (std::size_t i = 0; i < n; i++)
      error += std::pow(x[i] - x_old[i], 2);

    if (std::sqrt(error) < m_tolerance && iter > 0) {
      m_last_iterations_no = iter;
      break;
    }
  }
  m_last_tolerance = error;
}