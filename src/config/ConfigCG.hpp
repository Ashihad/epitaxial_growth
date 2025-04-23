#pragma once

struct ConfigCG {
  /* Conjugate Gradient method specific parameters */
  const long unsigned max_iterations{1000};
  const double tolerance{1.0E-3};
};