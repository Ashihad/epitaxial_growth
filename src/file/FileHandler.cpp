#include "FileHandler.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

FileHandler::FileHandler()
    : grid_filename{"grid.dat"},
      energy_filename{"energies.dat"},
      tmp_filename{"tmp.dat"},
      grid_fd{},
      energy_fd{},
      tmp_fd{} {}

FileHandler::~FileHandler() {}

// save grid above deposited atoms
void FileHandler::save_grid(const Grid& grid) {
  grid_fd.open(grid_filename, std::ios::out | std::ios::trunc);
  if (grid_fd.fail()) {
    std::cerr << "Unable to open grid file " << grid_filename << " for write"
              << std::endl;
    throw std::runtime_error("Cannot open grid file " + grid_filename);
  }
  grid_fd << std::scientific << std::setprecision(5);

  // header
  // clang-format off
  grid_fd << "#      i"
          << std::setw(8) << "j"
          << std::setw(8) << "atom"
          << std::setw(15) << "elastic_e" << '\n';
  // clang-format on
  std::size_t x_pos = 0;
  for (const auto& col : grid) {
    std::size_t y_pos = 0;
    for (const auto& atom : col) {
      if (atom.type != ATOM_TYPE::NO_ATOM) {
        // clang-format off
        grid_fd << std::setw(8) << x_pos
                << std::setw(8) << y_pos
                << std::setw(8) << static_cast<int>(atom.type)
                << std::setw(15) << atom.el_energy << '\n';
        // clang-format on
      }
      y_pos++;
    }
    x_pos++;
  }
  grid_fd.close();
}

// save elastic energy
void FileHandler::save_elastic_energy(const Grid& grid) {
  energy_fd.open(energy_filename, std::ios::out | std::ios::trunc);
  if (energy_fd.fail()) {
    std::cerr << "Unable to open energy file " << energy_filename
              << " for write" << std::endl;
    throw std::runtime_error("Cannot open energy file " + energy_filename);
  }
  energy_fd << std::scientific << std::setprecision(5);

  // header
  // clang-format off
  energy_fd << "#      i"
            << std::setw(8) << "j"
            << std::setw(8) << "atom"
            << std::setw(15) << "elastic_e"
            << std::setw(15) << "grad_loc_2" << '\n';
  // clang-format on

  // data
  std::size_t x_pos = 0;
  for (const auto& col : grid) {
    std::size_t y_pos = 0;
    for (const auto& atom : col) {
      const double grad_loc_2 =
          std::pow(atom.grad_x, 2) + std::pow(atom.grad_y, 2);
      // clang-format off
      energy_fd << std::setw(8) << x_pos
                << std::setw(8) << y_pos
                << std::setw(8) << static_cast<int>(atom.type)
                << std::setw(15) << atom.el_energy
                << std::setw(15) << grad_loc_2 << '\n';
      // clang-format on
      y_pos++;
    }
    energy_fd << '\n';
    x_pos++;
  }
  energy_fd.close();
}

// dump all data in grid
void FileHandler::save_tmp(const Grid& grid) {
  tmp_fd.open(tmp_filename, std::ios::out | std::ios::trunc);
  if (tmp_fd.fail()) {
    std::cerr << "Unable to open tmp file " << tmp_filename << " for write"
              << std::endl;
    throw std::runtime_error("Cannot open tmp file " + tmp_filename);
  }
  tmp_fd << std::scientific << std::setprecision(5);
  for (const auto& col : grid) {
    for (const auto& atom : col) {
      // clang-format off
      tmp_fd << std::setw(15)
             << static_cast<int>(atom.type) << "  "
             << atom.u << " "
             << atom.v << " "
             << atom.boundary1 << " "
             << atom.boundary2 << " "
             << atom.param_5 << " "
             << atom.param_6 << " "
             << atom.grad_x << " "
             << atom.grad_y << " "
             << atom.el_energy << " "
             << "\n";
      // clang-format on
    }
  }
  tmp_fd.close();
}
