#pragma once

#include <chrono>
#include <cstddef>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Atom.hpp"
#include "AtomContainers.hpp"
#include "ConfigCG.hpp"
#include "ConfigPhysics.hpp"
#include "ConfigSimulation.hpp"
#include "FileHandler.hpp"
#include "MathDriver.hpp"
#include "RngDriver.hpp"

class Simulator {
 public:
  Simulator(const ConfigPhysics&, const ConfigSimulation&, const ConfigCG&);
  void init_grid();
  void print_header();
  void print_iter_header();
  void perform_periodic_actions();

  double get_duv_max();

  void run_loop();

  double get_elastic_energy(const std::size_t, const std::size_t, bool = false);
  void conduct_local_relaxation(const std::size_t,
                                const std::size_t,
                                bool = false);
  void conduct_global_relaxation(bool = false);
  const double m_substrate_lattice_constant;
  const double m_adatom_lattice_constant;
  const double m_vertical_lat_spacing;
  const double m_spring_const_neighbors;
  const double m_spring_const_next_neighbors;
  const double m_misfit_coeff;
  const double m_D;
  const double m_E;
  const double m_bond_energy;

  const std::size_t m_grid_x;
  const std::size_t m_grid_y;
  const std::size_t m_substrate_height;

  const std::size_t m_local_relaxation_range_min;
  const std::size_t m_local_relaxation_range_max;
  const double m_local_relaxation_tolerance;

  const double m_time_fluence;
  const double m_simulation_time;

  const std::size_t m_diffusion_range;
  const double m_deposition_speed;
  const double m_temperature;

  const int m_diffusion_mode;
  const long unsigned m_dump_data_freq;
  const long unsigned m_global_relaxation_freq;

  // physics constants
  const double m_kbt;  // k_b * T / 1eV

  double m_duv_max;

  Grid m_grid;
  Grid m_grid_copy;
  DiffStructure m_atoms_diff;

  // drivers
  std::unique_ptr<MathDriver> mathDriver;
  RNGDriver rngDriver;
  FileHandler fHandler;
};