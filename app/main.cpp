#include "CGDriver.hpp"
#include "ConfigCG.hpp"
#include "ConfigPhysics.hpp"
#include "ConfigSimulation.hpp"
#include "Simulator.hpp"

int main() {
  Simulator sim{ConfigPhysics(), ConfigSimulation(), ConfigCG()};
  sim.run_loop();
}