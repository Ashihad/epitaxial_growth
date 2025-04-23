#pragma once

#include <vector>

#include "Atom.hpp"
#include "DiffusedAtom.hpp"

// main data structure
using Grid = std::vector<std::vector<Atom>>;
using DiffStructure = std::vector<DiffusedAtom>;
