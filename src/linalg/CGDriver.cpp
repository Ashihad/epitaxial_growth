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
  // if b is zero vector, return trivial solution x={0}
  double basically_zero{1.0E-10};
  double b_dot_b = dot(b, b);
  if (b_dot_b < std::abs(basically_zero)) {
    std::fill_n(x.get(), A.n_rows, 0.);
    m_last_iterations_no = 0;
    m_last_tolerance = 0.;
    return;
  }

  // residual vectors
  FastVector<double> rj(A.n_rows);
  FastVector<double> rj_proposed(A.n_rows);
  // approximate solution vectors
  FastVector<double> xj(A.n_rows);
  FastVector<double> xj_proposed(A.n_rows);
  // vectos for storing A*x0 and A*pj (temporary storage)
  FastVector<double> A_times_x0(A.n_rows);
  FastVector<double> A_times_pj(A.n_rows);
  // direction vectors
  FastVector<double> pj(A.n_rows);
  FastVector<double> pj_proposed(A.n_rows);

  // compute initial guess for A*x0, store in A_times_x0
  matrix_times_vector(A, x, A_times_x0);

  // initial guess
  std::memcpy(xj.get(), x.get(), x.size * sizeof(double));

  for (std::size_t i = 0; i < A.n_rows; i++) {
    rj[i] = b[i] - A_times_x0[i];
    pj[i] = rj[i];
  }

  double approximation_error{};
  // iterate for maximum of itmax iterations
  for (std::size_t j = 0; j < m_max_iterations; j++) {
    // compute A*pj
    matrix_times_vector(A, pj, A_times_pj);

    // compute step size alpha = r_j * r_j / (p_j * A * p_j)
    double rj_dot_rj = dot(rj, rj);
    double apj_dot_apj = dot(A_times_pj, pj);
    double alfa = rj_dot_rj / apj_dot_apj;
    if (std::fabs(alfa) < 1.0E-5)
      // step too small
      std::cerr << "Conjugate Gradient method error, step too small:  alfa="
                << alfa << std::endl;

    // update approximate solution vector
    for (std::size_t i = 0; i < A.n_rows; i++) {
      xj_proposed[i] = xj[i] + alfa * pj[i];
    }

    // update residual vector
    for (std::size_t i = 0; i < A.n_rows; i++) {
      rj_proposed[i] = rj[i] - alfa * A_times_pj[i];
    }

    // compute the update factor beta = r_{j+1} * r_{j+1} / (r_j * r_j)
    double rj_prop_dot_rj_prop = dot(rj_proposed, rj_proposed);
    double beta = rj_prop_dot_rj_prop / rj_dot_rj;
    for (std::size_t i = 0; i < A.n_rows; i++)
      pj_proposed[i] = rj_proposed[i] + beta * pj[i];

    // insert computed tmp solutions to corresponding vectors
    std::memcpy(xj.get(), xj_proposed.get(), xj_proposed.size * sizeof(double));
    std::memcpy(rj.get(), rj_proposed.get(), rj_proposed.size * sizeof(double));
    std::memcpy(pj.get(), pj_proposed.get(), pj_proposed.size * sizeof(double));

    // compute solution error
    approximation_error = std::sqrt(rj_dot_rj) / std::sqrt(b_dot_b);
    // if error is satisfying, end procedure
    if (approximation_error < m_tolerance && j > 0) {
      m_last_iterations_no = j;
      break;
    }
  }
  m_last_tolerance = approximation_error;

  // save computed solution
  std::memcpy(x.get(), xj.get(), xj.size * sizeof(double));
}  // CG-standard
