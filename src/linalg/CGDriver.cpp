#include "CGDriver.hpp"

#include <iostream>
#include <limits>
#include <vector>

#include "Simulator.hpp"

CGDriver::CGDriver(const ConfigMathDriver& conf_md,
                   const ConfigPhysics& conf_ph)
    : MathDriver{conf_md, conf_ph} {}

void CGDriver::solve_linear_system(const std::size_t row_count,
                                   double* csr_val,
                                   int* csr_row,
                                   int* csr_col,
                                   double* b,
                                   double* x) {
  // if b is zero vector, return trivial solution x={0,...}
  double f_2 = scalar_product(row_count, b, b);
  if (f_2 < 1.0E-10) {
    for (std::size_t i = 0; i < row_count; i++) {
      x[i] = 0.;
    }
    m_last_iterations_no = 0;
    m_last_tolerance = 0.;
    return;
  }

  // residual vectors
  std::unique_ptr<double[]> rj{new double[row_count]};
  std::unique_ptr<double[]> rj_proposed{new double[row_count]};
  // approximate solution vectors
  std::unique_ptr<double[]> xj{new double[row_count]};
  std::unique_ptr<double[]> xj_proposed{new double[row_count]};
  // vector for storing A*x0
  std::unique_ptr<double[]> A_times_x0{new double[row_count]};
  // matrix-vector product result (A*pj)
  std::unique_ptr<double[]> A_times_pj{new double[row_count]};
  // direction vectors
  std::unique_ptr<double[]> pj{new double[row_count]};
  std::unique_ptr<double[]> pj_proposed{new double[row_count]};

  // compute initial guess for A*x0, store in tmp1
  compute_sparse_Ax_y(row_count, csr_val, csr_row, csr_col, x,
                      A_times_x0.get());

  for (std::size_t i = 0; i < row_count; i++) {
    xj[i] = x[i];                  // initial guess
    rj[i] = b[i] - A_times_x0[i];  // b-A*x0
    pj[i] = rj[i];
  }

  double Apj_2;
  double rj_2;
  double rjp1_2;
  double approximation_error;
  double alfa;

  // iterate for maximum of itmax iterations
  for (std::size_t j = 0; j < m_max_iterations; j++) {
    // compute A*pj
    compute_sparse_Ax_y(row_count, csr_val, csr_row, csr_col, pj.get(),
                        A_times_pj.get());

    // compute step size alpha = r_j * r_j / (p_j * A * p_j)
    rj_2 = scalar_product(row_count, rj.get(), rj.get());
    Apj_2 = scalar_product(row_count, A_times_pj.get(), pj.get());
    alfa = rj_2 / Apj_2;
    if (std::fabs(alfa) < 1.0E-5)
      // step too small
      std::cerr << "Conjugate Gradient method error, step too small:  alfa="
                << alfa << std::endl;

    // update approximate solution vector
    for (std::size_t i = 0; i < row_count; i++) {
      xj_proposed[i] = xj[i] + alfa * pj[i];
    }

    // update residual vector
    for (std::size_t i = 0; i < row_count; i++) {
      rj_proposed[i] = rj[i] - alfa * A_times_pj[i];
    }

    // compute the update factor beta = r_{j+1} * r_{j+1} / (r_j * r_j)
    rjp1_2 = scalar_product(row_count, rj_proposed.get(), rj_proposed.get());
    double beta = rjp1_2 / rj_2;
    for (std::size_t i = 0; i < row_count; i++)
      pj_proposed[i] = rj_proposed[i] + beta * pj[i];

    // insert computed tmp solutions to corresponding vectors
    for (std::size_t i = 0; i < row_count; i++) {
      xj[i] = xj_proposed[i];
      rj[i] = rj_proposed[i];
      pj[i] = pj_proposed[i];
    }

    // compute solution error
    approximation_error = sqrt(rj_2) / sqrt(f_2);
    // if error is satisfying, end procedure
    if (approximation_error < m_tolerance && j > 0) {
      m_last_iterations_no = j;
      break;
    }
  }
  m_last_tolerance = approximation_error;

  // save computed solution
  for (std::size_t i = 0; i < row_count; i++)
    x[i] = xj[i];
}  // CG-standard
