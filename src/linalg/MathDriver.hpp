#pragma once

#include <array>
#include <cstdlib>

class Simulator;

// maksymalna liczba elementow w wierszu - liczba sasiadow * liczba kierunkow
constexpr std::size_t column_count = 9 * 2;

class MathDriver {
 public:
  MathDriver(Simulator*);
  virtual ~MathDriver() = default;
  MathDriver(const MathDriver&) = delete;
  MathDriver(MathDriver&&) = delete;
  MathDriver& operator=(const MathDriver&) = delete;
  MathDriver&& operator=(MathDriver&&) = delete;
  // must be defined in child
  virtual void solve_linear_system(const std::size_t,
                                   const std::size_t,
                                   const std::size_t,
                                   const std::size_t,
                                   const int,
                                   double*,
                                   bool = false) = 0;

  Simulator* sim;
};
