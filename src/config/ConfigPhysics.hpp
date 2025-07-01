#pragma once

#include <cstdint>

// parameters taken from
// https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.89.216102

/* Physical properties of substrate and adatom materials */

struct ConfigPhysics {
  virtual ~ConfigPhysics() = default;
  // a_{Si}
  const double substrate_lattice_constant{5.43 / 2};  // [angstrem
  // a_{Ge}
  const double adatom_lattice_constant{5.66 / 2};  // [angstrem]
  // mu
  const double misfit_coeff{
      (adatom_lattice_constant - substrate_lattice_constant) /
      substrate_lattice_constant};
  // k_{N}
  const double spring_const_neighbors{
      15.85 / substrate_lattice_constant /
      substrate_lattice_constant};  //[eV/angstrem^2] - stala sprzezystosci
  // k_{NN}
  const double spring_const_next_neighbors{spring_const_neighbors /
                                           2.};  //[Ha/ab^2]

  // vertical lattice spacing
  const double vertical_lat_spacing{
      adatom_lattice_constant +
      substrate_lattice_constant * misfit_coeff * spring_const_next_neighbors /
          (spring_const_neighbors +
           spring_const_next_neighbors)};  // vertical lattice spacing

  // fitting factor, from experimental data
  const double D0{3.83E+13};  //[angstrem^2/second] -
  const double D{D0 / substrate_lattice_constant /
                 substrate_lattice_constant};  //[1/second]

  // gamma
  const double bond_energy{0.4};  // [eV]
  // fitting factor, from experimental data
  const double E{0.53};  // [eV]

  const double deposition_speed{0.08};  // [layers/second]
  const double temperature{600};        // [K]
};
