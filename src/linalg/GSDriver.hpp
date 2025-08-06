#pragma once

#include "MathDriver.hpp"

class GSDriver : public MathDriver {
 public:
  GSDriver(const ConfigMathDriver&, const ConfigPhysics&);
  virtual ~GSDriver() = default;
  GSDriver(const GSDriver&) = delete;
  GSDriver(GSDriver&&) = delete;
  GSDriver& operator=(const GSDriver&) = delete;
  GSDriver&& operator=(GSDriver&&) = delete;

  /**
   * Solve Ax = b in CSR format, Gauss-Seidel algorithm
   *
   * @param A - matrix
   * @param b - free vector
   * @param x - solution vector
   */
  void solve_linear_system(const CSRMatrix& A,
                           const FastVector<double>& b,
                           FastVector<double>& x) override;
};