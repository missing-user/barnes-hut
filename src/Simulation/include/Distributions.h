#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

#include <vector>

#include "Simulation.h"

enum class Distribution {
  UNIVERSE1,
  UNIVERSE2,
  UNIVERSE3,
  UNIVERSE4,
  COLLISION,
  BIGBANG,
  STABLE_ORBIT
};

void set_seed(unsigned int seed);

std::vector<Particle> &set_mass(std::vector<Particle> &particles, myfloat m);

/*******************************************************************/

std::vector<Particle> make_universe(Distribution dist, int num_particles);

#endif
