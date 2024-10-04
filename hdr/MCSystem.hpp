#include "Atom.hpp"

#include <array>

class MCSystem2D {
  public:
    MCSystem2D();
  private:
    std::array<std::array<Atom, 5>, 5> space;
};