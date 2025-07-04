#include "CGDriver.hpp"

#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

#include "Simulator.hpp"

CGDriver::CGDriver(const ConfigMathDriver& conf_md,
                   const ConfigPhysics& conf_ph)
    : MathDriver{conf_md, conf_ph} {}

// Solve Ax = b using Conjugate gradient algorithm
// https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm
void CGDriver::solve_linear_system(CSRMatrix& A,
                                   FastVector<double>& b,
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

  // residual vector
  FastVector<double> r(n);
  // vector for storing A*pj (temporary storage)
  FastVector<double> A_times_pj(n);
  // direction vector
  FastVector<double> p(n);

  // r = b - A * x
  matrix_times_vector(A, x, r);
  for (std::size_t i = 0; i < A.n_rows; i++) {
    r[i] = b[i] - r[i];
  }

  // p = r
  std::memcpy(p.get(), r.get(), r.size * sizeof(double));

  double approximation_error{};
  double r_dot_r{dot(r, r)};
  // iterate for maximum of itmax iterations
  for (std::size_t iter = 0; iter < m_max_iterations; iter++) {
    // compute A*pj
    matrix_times_vector(A, p, A_times_pj);

    // compute step size alpha = r_j * r_j / (p_j * A * p_j)
    double alfa = r_dot_r / dot(A_times_pj, p);
    if (std::fabs(alfa) < 1.0E-5) {
      // step too small, abort
      std::cerr << "Conjugate Gradient method error, step too small:  alfa="
                << alfa << std::endl;
      std::abort();
    }

    // update solution vector
    for (std::size_t i = 0; i < n; i++) {
      x[i] += alfa * p[i];
    }

    // update residual vector
    for (std::size_t i = 0; i < n; i++) {
      r[i] -= alfa * A_times_pj[i];
    }

    // compute the update factor beta = r_{j+1} * r_{j+1} / (r_j * r_j)
    double r_dot_r_new = dot(r, r);
    double beta = r_dot_r_new / r_dot_r;
    for (std::size_t i = 0; i < n; i++)
      p[i] = r[i] + beta * p[i];

    r_dot_r = r_dot_r_new;

    // compute solution error
    approximation_error = std::sqrt(r_dot_r) / std::sqrt(b_dot_b);
    // if error is satisfying, end procedure
    if (approximation_error < m_tolerance && iter > 0) {
      m_last_iterations_no = iter;
      break;
    }
  }
  m_last_tolerance = approximation_error;
}  // CG-standard
