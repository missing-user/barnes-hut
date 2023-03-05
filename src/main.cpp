#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Distributions.h"
#include "Order.h"
#include "Particle.h"
#include "Simulation.h"

int main(int argc, char *argv[]) {
  bool output_csv = false;
  bool brute_force = false;
  bool no_reorder = true;
  int num_particles = 1000;
  myfloat theta = 1.5;
  double duration = 10;
  myfloat timestep = .1;


  set_seed(42);
  std::vector<Particle> particles =
      make_universe(Distribution::UNIVERSE4, num_particles);

  if (!no_reorder && !brute_force)
    computeAndOrder(particles);


  simulate(particles, duration, timestep, brute_force, theta);

  return 0;
}