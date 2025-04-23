#pragma once

#include "AtomTypes.hpp"

// main data container structure
enum STRUCTURE { ATOM_TYPE = 0, U = 1, V = 2, ELASTIC_ENERGY = 8 };

struct Atom {
  Atom()
      : type{ATOM_TYPE::NO_ATOM},
        u{},
        v{},
        boundary1{},
        boundary2{},
        param_5{},
        param_6{},
        grad_x{},
        grad_y{},
        el_energy{} {};
  Atom(enum ATOM_TYPE atom_type)
      : type{atom_type},
        u{},
        v{},
        boundary1{},
        boundary2{},
        param_5{},
        param_6{},
        grad_x{},
        grad_y{},
        el_energy{} {}
  enum ATOM_TYPE type;  // param 0
  double u;             // param 1
  double v;             // param 2
  int boundary1;        // param 3
  int boundary2;        // param 4
  double param_5;
  double param_6;
  double grad_x;
  double grad_y;
  double el_energy;
};
