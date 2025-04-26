#include "FileHandler.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

FileHandler::FileHandler()
    : grid_filename{"grid.dat"},
      energy_filename{"energies.dat"},
      tmp_filename{"tmp.dat"},
      grid_fd{grid_filename},
      energy_fd{energy_filename},
      tmp_fd{tmp_filename} {
  if (grid_fd.fail()) {
    std::cerr << "Unable to open grid file " << grid_filename << " for write"
              << std::endl;
    throw std::runtime_error("Cannot open grid file " + grid_filename);
  }
  if (energy_fd.fail()) {
    std::cerr << "Unable to open energy file " << energy_filename
              << " for write" << std::endl;
    throw std::runtime_error("Cannot open energy file " + energy_filename);
  }
  if (tmp_fd.fail()) {
    std::cerr << "Unable to open tmp file " << tmp_filename << " for write"
              << std::endl;
    throw std::runtime_error("Cannot open tmp file " + tmp_filename);
  }
  grid_fd << std::scientific << std::setprecision(5);
  energy_fd << std::scientific << std::setprecision(5);
  tmp_fd << std::scientific << std::setprecision(5);
}

FileHandler::~FileHandler() {
  // close file descriptors
  grid_fd.close();
  energy_fd.close();
  tmp_fd.close();
}

// save grid above deposited atoms
void FileHandler::save_grid(const Grid& grid) {
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
  std::flush(grid_fd);
}

// save elastic energy???
void FileHandler::save_elastic_energy(const Grid& grid) {
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
  std::flush(energy_fd);
}

void FileHandler::save_tmp(const Grid& grid) {
  // TODO: Numbers Jason, what do they mean?
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
  std::flush(tmp_fd);
}
