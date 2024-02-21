#ifndef BARNESHUT
#define BARNESHUT

#include "Cuboid.h"
#include "Particle.h"
#include "Forces.h"
void stepSimulation(Particles& particles, myfloat dt, myfloat theta2);

struct debug_information
{
  int depth;
  int max_particles_in_leaf;
  std::vector<DrawableCuboid> debug_boxes;
};
debug_information bh_superstep_debug(myvec3 position, Particles& particles, size_t count, myfloat theta2);

#endif