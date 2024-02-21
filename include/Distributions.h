#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

#include "Particle.h"
#include <vector>

enum class Distribution {
  UNIVERSE1,
  UNIVERSE2,
  UNIVERSE3,
  UNIVERSE4,
  COLLISION,
  BIGBANG,
  STABLE_ORBIT,
  SPHERE,
  CRYSTALLINE,
  PLUMMER,
  DEBUG_CUBE
};

void set_seed(unsigned int seed);

std::vector<Particle> &set_mass(std::vector<Particle> &particles, myfloat m);

/*******************************************************************/

std::vector<Particle> make_universe(Distribution dist, size_t num_particles);

#endif
