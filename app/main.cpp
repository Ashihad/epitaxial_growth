#include "CGDriver.hpp"
#include "ConfigMathDriver.hpp"
#include "ConfigPhysics.hpp"
#include "ConfigSimulation.hpp"
#include "Simulator.hpp"

int main() {
  Simulator sim{ConfigPhysics(), ConfigSimulation(), ConfigMathDriver()};
  sim.run_loop();
}