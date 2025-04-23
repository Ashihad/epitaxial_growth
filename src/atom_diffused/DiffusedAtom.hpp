#pragma once

#include "AtomTypes.hpp"

#include <cstdlib>

struct DiffusedAtom {
  enum ATOM_TYPE type;  // param 0

  // distances? positions?
  std::size_t idx;         // param 1
  std::size_t idy;         // param 2
  unsigned int neighbors;  // param 3
  double el_energy;        // param 4
  // we put rate of diffusion here
  double r_diff;  // param 5
  // we put copy of rate of diffusion here
  double r_diff_copy;           // param 6
  unsigned int neighbors_copy;  // param 7 ??????
  double delta_w;               // param 8
};
