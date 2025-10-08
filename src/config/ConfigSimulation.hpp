#pragma once

#include <cstdlib>

struct ConfigSimulation {
  virtual ~ConfigSimulation() = default;
  const std::size_t grid_x{250};
  const std::size_t grid_y{80};
  const std::size_t substrate_height{16};
  const std::size_t initial_adatom_height{0};

  // relaxation parameters
  const std::size_t diffusion_range{8};
  const std::size_t local_relaxation_range_min{10};
  const std::size_t local_relaxation_range_max{10};
  const double local_relaxation_tolerance{1.0E-2};
  // global relaxation will be performed every n_global_relaxation iterations
  const long unsigned global_relaxation_freq{1000};

  // Legacy
  // 1 indicates diffusion of both substrate and adatoms
  // 2 indicates diffusion of only adatoms
  // TODO: does 1 work?
  const int diffusion_mode{2};
  const double time_fluence{8.0};
  const double simulation_time{2.};

  // every dump_data_freq data is saved to files
  const long unsigned dump_data_freq{1000};

  // whether to add island which should modify growth region shape
  const bool add_island{false};
};
