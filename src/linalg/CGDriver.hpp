#pragma once

#include "MathDriver.hpp"

class CGDriver : public MathDriver {
 public:
  CGDriver(const ConfigMathDriver&, const ConfigPhysics&);
  virtual ~CGDriver() = default;
  CGDriver(const CGDriver&) = delete;
  CGDriver(CGDriver&&) = delete;
  CGDriver& operator=(const CGDriver&) = delete;
  CGDriver&& operator=(CGDriver&&) = delete;

  /**
   * solve Ax = b, Conjugate gradient algorithm
   * https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm
   *
   * @param A - matrix
   * @param b - free vector
   * @param x - solution vector
   */
  void solve_linear_system(const CSRMatrix& A,
                           const FastVector<double>& b,
                           FastVector<double>& x) override;
};