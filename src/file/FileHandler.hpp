#pragma once

#include "AtomContainers.hpp"

#include <fstream>
#include <string>

class FileHandler {
 public:
  FileHandler();
  ~FileHandler();
  FileHandler(const FileHandler&) = delete;
  FileHandler(FileHandler&&) = delete;
  FileHandler& operator=(const FileHandler&) = delete;
  FileHandler&& operator=(FileHandler&&) = delete;

  void save_grid(const Grid&);
  void save_elastic_energy(const Grid&);
  void save_tmp(const Grid&);

  // filenames
  const std::string grid_filename;
  const std::string energy_filename;
  const std::string tmp_filename;
  // file descriptors
  std::ofstream grid_fd;
  std::ofstream energy_fd;
  std::ofstream tmp_fd;
};
