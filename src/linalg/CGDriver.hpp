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
   *	solve linear equations system: CG - standard algorithm - Saad
   *
   * @param n
   * @param acsr
   * @param icsr
   * @param jcsr
   * @param b
   * @param x
   * @param itmax
   * @param tol
   */
  void solve_linear_system(const std::size_t row_count,
                           double* csr_val,
                           int* csr_row,
                           int* csr_col,
                           double* b,
                           double* x) override;
};