#include "Simulator.hpp"

#include "CGDriver.hpp"
#include "GSDriver.hpp"
#include "MatrixTypes.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

Simulator::Simulator(const ConfigPhysics& conf_ph,
                     const ConfigSimulation& conf_sim,
                     const ConfigMathDriver& conf_cg)
    : m_substrate_lattice_constant{conf_ph.substrate_lattice_constant},
      m_adatom_lattice_constant{conf_ph.adatom_lattice_constant},
      m_vertical_lat_spacing{conf_ph.vertical_lat_spacing},
      m_spring_const_neighbors{conf_ph.spring_const_neighbors},
      m_spring_const_next_neighbors{conf_ph.spring_const_next_neighbors},
      m_misfit_coeff{conf_ph.misfit_coeff},
      m_D{conf_ph.D},
      m_E{conf_ph.E},
      m_bond_energy{conf_ph.bond_energy},
      m_grid_x{conf_sim.grid_x},
      m_grid_y{conf_sim.grid_y},
      m_substrate_height{conf_sim.substrate_height},
      m_initial_adatom_height{conf_sim.initial_adatom_height},
      m_local_relaxation_range_min{conf_sim.local_relaxation_range_min},
      m_local_relaxation_range_max{conf_sim.local_relaxation_range_max},
      m_local_relaxation_tolerance{conf_sim.local_relaxation_tolerance},
      m_time_fluence{conf_sim.time_fluence},
      m_simulation_time{conf_sim.simulation_time},
      m_diffusion_range{conf_sim.diffusion_range},
      m_deposition_speed{conf_ph.deposition_speed},
      m_temperature{conf_ph.temperature},
      m_diffusion_mode{conf_sim.diffusion_mode},
      m_dump_data_freq{conf_sim.dump_data_freq},
      m_global_relaxation_freq{conf_sim.global_relaxation_freq},
      m_island{conf_sim.add_island},
      m_kbt{1.38E-23 * m_temperature / 1.602E-19},
      m_duv_max{},
      m_grid{},
      m_grid_copy{},
      m_atoms_diffused{},
      mathDriver{new CGDriver(conf_cg, conf_ph)},
      rngDriver{m_grid_x, m_diffusion_range},
      fHandler{} {
  init_grid();
}

void Simulator::init_grid() {
  // dimensionality: (m_grid_x, m_grid_y, n_values)
  m_grid.resize(m_grid_x, std::vector<Atom>(m_grid_y, Atom()));
  m_grid_copy.resize(m_grid_x, std::vector<Atom>(m_grid_y, Atom()));
  // deposition probability + no of atoms*diffusion probability
  m_atoms_diffused.resize(m_grid_x + 1);

  // setting up substrate atoms up to m_substrate_height
  for (auto& col : m_grid) {
    for (std::size_t j = 0; j < m_substrate_height; ++j) {
      col[j].type = ATOM_TYPE::SUBSTRATE;  // 0-empty, 1-Si, 2-Ge
    }
  }

  // setting up five monolayers of adatoms
  for (auto& col : m_grid) {
    for (std::size_t j = m_substrate_height;
         j < m_substrate_height + m_initial_adatom_height; ++j) {
      col[j].type = ATOM_TYPE::ADATOM;  // 0-empty, 1-Si, 2-Ge
    }
  }
  if (m_island)
    add_island();
}

void Simulator::add_island() {
  const std::size_t x_min{m_grid_x / 2 - 15 - 40};
  const std::size_t x_max{m_grid_x / 2 + 15 + 1 - 40};
  const std::size_t y_max{m_substrate_height - 2 + 1};
  const std::size_t y_min{m_substrate_height - 12};

  for (std::size_t i{x_min}; i < x_max; i++) {
    for (std::size_t j{y_min}; j < y_max; j++) {
      m_grid[i][j].type = ATOM_TYPE::ADATOM;
    }
  }
}

void Simulator::print_header() {
  std::cout << "KMC 1+1 simulation started\n\n";
  std::cout << "Simulation parameters:\n";
  std::cout
      << "=====================================================\n"
      << std::fixed << std::setprecision(5) << std::setw(40)
      << "substrate_lattice_constant [A] = " << m_substrate_lattice_constant
      << "\n"
      << std::setw(40)
      << "adatom_lattice_constant [A] = " << m_adatom_lattice_constant << "\n"
      << std::setw(40)
      << "vertical_lat_spacing [A] = " << m_vertical_lat_spacing << "\n"
      << std::setw(40)
      << "spring_const_neighbors [eV/A/A] = " << m_spring_const_neighbors
      << "\n"
      << std::setw(40) << "spring_const_next_neighbors [eV/A/A] = "
      << m_spring_const_next_neighbors << "\n"
      << std::setw(40) << "misfit_coeff [A] = " << m_misfit_coeff << "\n"
      << std::setw(40) << "D [1/s] = " << std::scientific << m_D << "\n"
      << std::setw(40) << "E [eV] = " << std::fixed << m_E << "\n"
      << std::setw(40) << "bond_energy [eV] = " << m_bond_energy << "\n"
      << std::setw(40) << "grid_x = " << m_grid_x << "\n"
      << std::setw(40) << "grid_y = " << m_grid_y << "\n"
      << std::setw(40) << "substrate_height = " << m_substrate_height << "\n"
      << std::setw(40)
      << "local_relaxation_range_min = " << m_local_relaxation_range_min << "\n"
      << std::setw(40)
      << "local_relaxation_range_max = " << m_local_relaxation_range_max << "\n"
      << std::setw(40) << "time_fluence [s] = " << std::setprecision(3)
      << m_time_fluence << "\n"
      << std::setw(40) << "simulation_time [s] = " << m_simulation_time << "\n"
      << std::setw(40) << "diffusion_range = " << m_diffusion_range << "\n"
      << std::setw(40) << "deposition_speed [ML/s] = " << m_deposition_speed
      << "\n"
      << std::setw(40) << "temperature [K] = " << std::setprecision(0)
      << m_temperature << "\n"
      << "=====================================================\n\n";
}

void Simulator::print_iter_header() {
  // clang-format off
  std::cout << std::setw(15) << "iter"
            << std::setw(15) << "sim_time"
            << std::setw(15) << "new_atoms"
            << std::setw(15) << "pr_acc"
            << std::setw(15) << "real_time[s]"
            << std::setw(15) << "duv_max" << std::endl;
  // clang-format on
}

double Simulator::get_duv_max() {
  // calculate duv_max for all atoms in the grid
  double duv_max{};
  for (const auto& col : m_grid) {
    for (const auto& atom : col) {
      double duv = std::sqrt(std::pow(atom.u, 2) + std::pow(atom.v, 2));
      if (duv > duv_max && atom.type != ATOM_TYPE::NO_ATOM)
        duv_max = duv;
    }
  }
  return duv_max;
}

void Simulator::perform_periodic_actions() {
  // calculate elastic energy for each atom
  std::size_t x_pos{};
  for (auto& col : m_grid) {
    std::size_t y_pos{};
    for (auto& atom : col) {
      atom.el_energy = get_elastic_energy(x_pos, y_pos);
      y_pos++;
    }
    x_pos++;
  }

  // get dE/duv
  calculate_dEduv();

  // dump data
  fHandler.save_grid(m_grid);
  fHandler.save_elastic_energy(m_grid);
  fHandler.save_tmp(m_grid);

  m_duv_max = get_duv_max();
}

void Simulator::run_loop() {
  srand(0);

  print_header();
  print_iter_header();

  // counters
  long unsigned added_atoms_count = 0;
  long unsigned proposed_diffusions = 1;
  long unsigned accepted_diffusions = 1;

  // timers
  double sim_time_passed{};
  long unsigned sim_iter{};
  auto start_time = std::chrono::high_resolution_clock::now();

  while (sim_time_passed < m_simulation_time) {
    sim_iter++;
    if (sim_iter % m_dump_data_freq == 0) {
      auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now() - start_time);
      perform_periodic_actions();

      // print iteration sanity-check info
      // clang-format off
      std::cout << std::setw(15) << sim_iter
                << std::fixed << std::setprecision(6) << std::setw(15) << sim_time_passed
                << std::setw(15) << added_atoms_count
                << std::setprecision(4) << std::setw(15) << (static_cast<double>(accepted_diffusions) / static_cast<double>(proposed_diffusions))
                << std::setprecision(1) << std::setw(15) << (static_cast<double>(time_elapsed.count()) * 1.0E-3)
                << std::setprecision(3) << std::setw(15) << m_duv_max << '\n';
      // clang-format on
    }

    // perform global relaxation every global_relaxation_freq iterations
    if (sim_iter % m_global_relaxation_freq == 0 || sim_iter == 1) {
      conduct_global_relaxation();
    }

    // okreslamy tempo depozycji atomow: deposition rate
    // rate= fluence*nx -> tempo (prawdopodobienstwo/sekunde) depozycji
    // jednego atomu w calym ukladzie (gdziekolwiek)
    double fluence = sim_time_passed < m_time_fluence
                         ? m_deposition_speed * static_cast<double>(m_grid_x)
                         : 0;

    double r_deposition = (static_cast<double>(m_diffusion_range) + 1.0) *
                          (2 * static_cast<double>(m_diffusion_range) + 1.0) *
                          fluence / 6.;

    // m_atoms_diffused[0] holds information about deposition probability, next
    // elements holds probabilities for diffusion of specific atoms
    m_atoms_diffused[0].r_diff = r_deposition;
    m_atoms_diffused[0].r_diff_copy = r_deposition;

    // calculate how many atoms can be diffused
    std::size_t n_at_diff = 0;
    for (std::size_t i = 0; i < m_grid_x; i++) {
      // search for first existing atom, top to bottom
      for (std::size_t j = m_grid_y - 2; j >= 1; j--) {
        enum ATOM_TYPE curr_atom_type = m_grid[i][j].type;
        if (curr_atom_type != ATOM_TYPE::NO_ATOM) {
          // is this atom type diffundable by the simulator config?
          if (static_cast<int>(curr_atom_type) >= m_diffusion_mode) {
            n_at_diff++;
            m_atoms_diffused[n_at_diff].type = curr_atom_type;
            m_atoms_diffused[n_at_diff].idx = i;
            m_atoms_diffused[n_at_diff].idy = j;

            // calculate number of neighbors (nearest + next nearest)
            unsigned neighbors_total = 0;
            for (int i_neigh = -1; i_neigh <= 1; i_neigh++) {
              for (int j_neigh = -1; j_neigh <= 1; j_neigh++) {
                if ((std::abs(i_neigh) + std::abs(j_neigh)) > 0) {
                  std::size_t i2 =
                      static_cast<std::size_t>(static_cast<int>(i) + i_neigh +
                                               static_cast<int>(m_grid_x)) %
                      static_cast<std::size_t>(m_grid_x);
                  std::size_t j2 =
                      static_cast<std::size_t>(static_cast<int>(j) + j_neigh);
                  if (m_grid[i2][j2].type != ATOM_TYPE::NO_ATOM)
                    neighbors_total++;  // neighbor found
                }
              }
            }
            m_atoms_diffused[n_at_diff].neighbors = neighbors_total;
            m_atoms_diffused[n_at_diff].neighbors_copy = neighbors_total;

            double ep = get_elastic_energy(i, j);
            m_atoms_diffused[n_at_diff].el_energy = ep;

            double r0 = 12. * m_D /
                        (static_cast<double>(m_diffusion_range) + 1.0) /
                        (2 * static_cast<double>(m_diffusion_range) + 2.);

            // change effective energy based on no of neighbors
            double delta_w = ep;
            if (neighbors_total == 3) {
              delta_w = ep * 1.5;
            } else if (neighbors_total == 4) {
              delta_w = ep * 2.0;
            } else if (neighbors_total >= 5) {
              delta_w = ep * 3.5;
            }

            // diffusion rate of current atom
            double ri =
                r0 *
                std::exp((-m_bond_energy * neighbors_total + delta_w + m_E) /
                         m_kbt);

            m_atoms_diffused[n_at_diff].r_diff = ri;
            m_atoms_diffused[n_at_diff].r_diff_copy = ri;

            m_atoms_diffused[n_at_diff].delta_w = delta_w;
          }
          break;
        }
      }
    }

    // increment time
    double sum_ri = 0;
    for (std::size_t i = 0; i <= n_at_diff; i++) {
      sum_ri += m_atoms_diffused[i].r_diff;
    }
    sim_time_passed += 1. / sum_ri;

    // calculate probability array for diffusions
    for (std::size_t i = 1; i <= n_at_diff; i++) {
      m_atoms_diffused[i].r_diff += m_atoms_diffused[i - 1].r_diff;
    }

    // choose a random process
    double u1 = rngDriver.gen_uniform() * sum_ri;
    std::size_t chosen_action = n_at_diff;
    for (std::size_t i = 0; i <= n_at_diff; i++) {
      if (u1 <= m_atoms_diffused[i].r_diff) {
        chosen_action = i;
        break;
      }
    }

    // chosen process - deposition
    if (chosen_action == 0) {
      added_atoms_count++;

      // get random x coordinate
      std::size_t i_added = rngDriver.gen_discrete_1_grid_x() - 1;

      // prevent atoms from overflowing upwards from grid
      if (m_grid[i_added][m_grid_y - 1].type != ATOM_TYPE::NO_ATOM) {
        std::cerr << "Too many atoms in system, omitting point in i = "
                  << i_added << '\n';
      }

      // get y coordinate of added atom by scanning for nearest atom in chosen x
      // coordinate, top to down, add adatom one y upward than found
      std::size_t j_added = 0;
      std::size_t j{m_grid_y - 2};
      while (j != 0) {
        if (m_grid[i_added][j].type != ATOM_TYPE::NO_ATOM) {
          j_added = j + 1;
          m_grid[i_added][j_added].type = ATOM_TYPE::ADATOM;
          m_grid[i_added][j_added].u = 0.;
          m_grid[i_added][j_added].v = 0.;
          break;
        }
        j--;
      }

      conduct_local_relaxation(i_added, j_added);
    }

    // chosen process - diffusion
    else {
      std::size_t i_old = m_atoms_diffused[chosen_action].idx;
      std::size_t j_old = m_atoms_diffused[chosen_action].idy;

      // calculate energy with atom in, in range
      // [i_old +/- local_en_range, j_old +/- local_en_range]
      int local_en_range = 25;
      double en_old = 0.;
      for (int i = -local_en_range; i <= local_en_range; i++) {
        for (int j = -local_en_range; j <= local_en_range; j++) {
          std::size_t ii = static_cast<std::size_t>(
              (static_cast<int>(i_old) + i + static_cast<int>(m_grid_x)) %
              static_cast<int>(m_grid_x));
          std::size_t jj =
              static_cast<std::size_t>(static_cast<int>(j_old) + j);
          // if we didn't leave grid bounds in y direction and atom exist on
          // position [ii][jj]
          if (jj > 1 && jj < (m_grid_y - 1) &&
              m_grid[ii][jj].type != ATOM_TYPE::NO_ATOM) {
            en_old += get_elastic_energy(ii, jj);
          }
        }
      }

      // calculate energy w/out atom in range
      // [i_old +/- local_en_range, j_old +/- local_en_range],
      // use copy of grid for that
      m_grid_copy = m_grid;
      m_grid_copy[i_old][j_old].type = ATOM_TYPE::NO_ATOM;

      // recalculate displacements with no atom
      conduct_local_relaxation(i_old, j_old, true);

      // calculate energy with no atom, from copy of grid
      double en_new = 0.;
      for (int i = -local_en_range; i <= local_en_range; i++) {
        for (int j = -local_en_range; j <= local_en_range; j++) {
          int ii = static_cast<int>(
              (i_old + static_cast<std::size_t>(i) + m_grid_x) % m_grid_x);
          int jj = static_cast<int>(j_old + static_cast<std::size_t>(j));
          if (jj > 1 && static_cast<std::size_t>(jj) < (m_grid_y - 1) &&
              m_grid_copy[static_cast<std::size_t>(ii)]
                         [static_cast<std::size_t>(jj)]
                             .type != ATOM_TYPE::NO_ATOM) {
            en_new += get_elastic_energy(static_cast<std::size_t>(ii),
                                         static_cast<std::size_t>(jj), true);
          }
        }
      }

      // sprawdzamy czy prawdopodobienstwo r_new<r_old=r_approx - jesli
      // tak to atom dyfunduje

      // calculate probability of diffusion
      unsigned int neigh_d = m_atoms_diffused[chosen_action].neighbors_copy;
      double r0 = 12. * m_D / (static_cast<double>(m_diffusion_range) + 1.0) /
                  (2 * static_cast<double>(m_diffusion_range) + 2.);
      double ri_new =
          r0 *
          std::exp((-m_bond_energy * neigh_d + (en_old - en_new) / 2 + m_E) /
                   m_kbt);  // exact atom diffusion rate
      double ri_old = m_atoms_diffused[chosen_action]
                          .r_diff_copy;  // estimated diffusion probability

      proposed_diffusions++;

      // calculate if diffusion is accepted based on rates calculated above
      double random_uniform{rngDriver.gen_uniform()};
      double probability_border{ri_new / ri_old};
      if (random_uniform < (probability_border)) {
        accepted_diffusions++;

        // generate random shift in x direction
        int ii_shift = rngDriver.gen_discrete_plus_minus_diffusion_range();

        // new position, mind the sign in x dimension
        std::size_t i_new =
            (static_cast<std::size_t>(static_cast<int>(i_old) + ii_shift) +
             m_grid_x) %
            (m_grid_x);
        std::size_t j_new;

        // m_grid = m_grid_copy;

        // atom type shall not change
        enum ATOM_TYPE atom_type_tmp = m_grid[i_old][j_old].type;

        // delete atom in old position
        m_grid[i_old][j_old].type = ATOM_TYPE::NO_ATOM;

        // get new y coordinate of diffunded atom by scanning for nearest atom
        // in chosen x coordinate, top to down, add atom one j up
        for (std::size_t j = m_grid_y - 1; j >= 1; j--) {
          if (m_grid[i_new][j].type != ATOM_TYPE::NO_ATOM) {
            if (j < (m_grid_y - 1)) {
              j_new = j + 1;  // move atom to new (free) position
              m_grid[i_new][j_new].type = atom_type_tmp;
              break;
            } else {
              std::cerr << "No space to move atom - no space in crystal "
                           "structure\n";
              break;
            }
          }
        }

        // remember to recalculate displacements after moving atom
        conduct_local_relaxation(i_new, j_new);
      }
    }  // chosen_action: diffusion
  }  // main loop
}

double Simulator::get_elastic_energy(const std::size_t x_pos,
                                     const std::size_t y_pos,
                                     bool performOnCopy) {
  Grid& chosen_grid = performOnCopy ? m_grid_copy : m_grid;
  if (chosen_grid[x_pos][y_pos].type == ATOM_TYPE::NO_ATOM) {
    return 0.;
  }

  // looks like we only care about 3x3 neighborhood
  static Matrix3x3I atom_mask;
  // helper matrices, define as in Smereka 1+1
  static Matrix3x3D d1;  // 0 if si-si, a_g-a_s otherwise
  static Matrix3x3D d2;  // 0 if si-si, a_g-a_l otherwise
  static Matrix3x3D u;   // local neighborhoods' U_x
  static Matrix3x3D v;   // local neighborhoods' U_y

  // filling local helper matrices
  // iterate over neighborhood
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      // retrieve offsets for crystal access

      // we iterate over (i-1, i, i+1), with wrapping
      std::size_t i3 = (x_pos + m_grid_x + i - 1) % m_grid_x;
      std::size_t j3 = (y_pos + m_grid_y + j - 1) % m_grid_y;
      if (chosen_grid[i3][j3].type != ATOM_TYPE::NO_ATOM)
        atom_mask[i][j] = 1;
      else
        atom_mask[i][j] = 0;

      u[i][j] = chosen_grid[i3][j3].u;  // copy U and V from crystal
      v[i][j] = chosen_grid[i3][j3].v;

      int bond_type = static_cast<int>(chosen_grid[x_pos][y_pos].type) *
                      static_cast<int>(chosen_grid[i3][j3].type);
      if (bond_type == 1) {
        // substrate-substrate
        d1[i][j] = 0;
        d2[i][j] = 0;
      } else if (bond_type == 2 || bond_type == 4) {
        // adatom-adatom (bond_type=4) or
        // adatom-substrate (bond_type=2)
        d1[i][j] = m_adatom_lattice_constant - m_substrate_lattice_constant;
        d2[i][j] = m_adatom_lattice_constant - m_vertical_lat_spacing;
      }
    }
  }

  std::size_t ii = 1;  // atom centralny
  std::size_t jj = 1;

  // from smereka 1+1
  double wxx =
      m_spring_const_neighbors / 2. * atom_mask[ii][jj] *
          atom_mask[ii + 1][jj] *
          std::pow(u[ii + 1][jj] - u[ii][jj] - d1[ii + 1][jj], 2) +
      m_spring_const_neighbors / 2. * atom_mask[ii][jj] *
          atom_mask[ii - 1][jj] *
          std::pow(u[ii - 1][jj] - u[ii][jj] + d1[ii - 1][jj], 2) +
      m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
          atom_mask[ii + 1][jj + 1] *
          std::pow(u[ii + 1][jj + 1] - u[ii][jj] - d1[ii + 1][jj + 1], 2) +
      m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
          atom_mask[ii - 1][jj - 1] *
          std::pow(u[ii - 1][jj - 1] - u[ii][jj] + d1[ii - 1][jj - 1], 2) +
      m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
          atom_mask[ii + 1][jj - 1] *
          std::pow(u[ii + 1][jj - 1] - u[ii][jj] - d1[ii + 1][jj - 1], 2) +
      m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
          atom_mask[ii - 1][jj + 1] *
          std::pow(u[ii - 1][jj + 1] - u[ii][jj] + d1[ii - 1][jj + 1], 2);

  double wyy =
      m_spring_const_neighbors / 2. * atom_mask[ii][jj] *
          atom_mask[ii][jj + 1] *
          std::pow(v[ii][jj + 1] - v[ii][jj] - d2[ii][jj + 1], 2) +
      m_spring_const_neighbors / 2. * atom_mask[ii][jj] *
          atom_mask[ii][jj - 1] *
          std::pow(v[ii][jj - 1] - v[ii][jj] + d2[ii][jj - 1], 2) +
      m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
          atom_mask[ii + 1][jj + 1] *
          std::pow(v[ii + 1][jj + 1] - v[ii][jj] - d2[ii + 1][jj + 1], 2) +
      m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
          atom_mask[ii - 1][jj - 1] *
          std::pow(v[ii - 1][jj - 1] - v[ii][jj] + d2[ii - 1][jj - 1], 2) +
      m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
          atom_mask[ii + 1][jj - 1] *
          std::pow(v[ii + 1][jj - 1] - v[ii][jj] + d2[ii + 1][jj - 1], 2) +
      m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
          atom_mask[ii - 1][jj + 1] *
          std::pow(v[ii - 1][jj + 1] - v[ii][jj] - d2[ii - 1][jj + 1], 2);

  // god bless Schwart's theorem of symmetry of second derivatives
  double wxy = m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
                   atom_mask[ii - 1][jj - 1] *
                   (u[ii - 1][jj - 1] - u[ii][jj] + d1[ii - 1][jj - 1]) *
                   (v[ii - 1][jj - 1] - v[ii][jj] + d2[ii - 1][jj - 1]) +
               m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
                   atom_mask[ii + 1][jj + 1] *
                   (u[ii + 1][jj + 1] - u[ii][jj] - d1[ii + 1][jj + 1]) *
                   (v[ii + 1][jj + 1] - v[ii][jj] - d2[ii + 1][jj + 1]) -
               m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
                   atom_mask[ii + 1][jj - 1] *
                   (u[ii + 1][jj - 1] - u[ii][jj] - d1[ii + 1][jj - 1]) *
                   (v[ii + 1][jj - 1] - v[ii][jj] + d2[ii + 1][jj - 1]) -
               m_spring_const_next_neighbors / 4. * atom_mask[ii][jj] *
                   atom_mask[ii - 1][jj + 1] *
                   (u[ii - 1][jj + 1] - u[ii][jj] + d1[ii - 1][jj + 1]) *
                   (v[ii - 1][jj + 1] - v[ii][jj] - d2[ii - 1][jj + 1]);

  return wxx + wyy + 2 * wxy;
}

void Simulator::conduct_local_relaxation(const std::size_t x_pos,
                                         const std::size_t y_pos,
                                         bool performOnCopy) {
  /*************************************************************************************
   * lokalna zmiana polozen atomow w otoczeniu irange_min - minimalizacja
   * naprezen nie ma sensu zwiekszac rozmiaru otoczenia i liczenia bledu bmax
   * wielokrotnie bo wyznaczenie Acsr trwa zawsze 2.5 ms
   *************************************************************************************/
  const std::size_t imin =
      (x_pos + m_grid_x - m_local_relaxation_range_min) % m_grid_x;
  const std::size_t jmin = std::max(y_pos - m_local_relaxation_range_min, 1ul);
  const std::size_t jmax =
      std::min(y_pos + m_local_relaxation_range_min, m_grid_y - 2ul);
  // roznica: imin->imax
  const std::size_t i_nodes = 2 * m_local_relaxation_range_min;
  Grid& chosen_grid = performOnCopy ? m_grid_copy : m_grid;

  double biggest_abs_bi =
      mathDriver->compute_displacements(chosen_grid, imin, i_nodes, jmin, jmax);

  // wartosc bledu lokalnego
  double mu = (m_adatom_lattice_constant - m_substrate_lattice_constant) /
              m_substrate_lattice_constant;
  // is tolerance so high that we need to do global relaxation?
  double tol_local = biggest_abs_bi / mu / m_substrate_lattice_constant /
                     m_spring_const_neighbors;
  if (tol_local > m_local_relaxation_tolerance) {
    conduct_global_relaxation(performOnCopy);
  }
}

void Simulator::conduct_global_relaxation(bool performOnCopy) {
  const std::size_t imin = 0;
  const std::size_t jmin = 1;
  const std::size_t jmax = m_grid_y - 2;
  const std::size_t imax = m_grid_x - 1;
  Grid& chosen_grid = performOnCopy ? m_grid_copy : m_grid;
  mathDriver->compute_displacements(chosen_grid, imin, imax, jmin, jmax);
}

double Simulator::calculate_dEduv() {
  double gradient_norm = 0;

  for (std::size_t i = 0; i < m_grid_x; i++) {
    for (size_t j = 0; j < m_grid_y; j++) {
      if (m_grid[i][j].type != ATOM_TYPE::NO_ATOM) {
        double u_old = m_grid[i][j].u;
        double v_old = m_grid[i][j].v;
        double delta = 0.01;  // krok do liczenia pochodnych

        m_grid[i][j].u = u_old + delta;
        double epx = get_elastic_energy(i, j);

        m_grid[i][j].u = u_old - delta;
        double emx = get_elastic_energy(i, j);

        m_grid[i][j].u = u_old;

        m_grid[i][j].v = v_old + delta;
        double epy = get_elastic_energy(i, j);

        m_grid[i][j].v = v_old - delta;
        double emy = get_elastic_energy(i, j);

        m_grid[i][j].v = v_old;

        m_grid[i][j].grad_x = (epx - emx) / 2 / delta;  // pochodna w x
        m_grid[i][j].grad_y = (epy - emy) / 2 / delta;  // pochodna w y

        gradient_norm +=
            std::pow(m_grid[i][j].grad_x, 2) + std::pow(m_grid[i][j].grad_y, 2);
      } else {
        m_grid[i][j].grad_x = 0.;
        m_grid[i][j].grad_y = 0.;
      }
    }
  }

  gradient_norm = std::sqrt(gradient_norm);

  return gradient_norm;
}
