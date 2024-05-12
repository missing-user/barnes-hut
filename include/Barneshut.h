#ifndef BARNESHUT
#define BARNESHUT

#include "Cuboid.h"
#include "Particle.h"
#include "Forces.h"

std::vector<DrawableCuboid> losing_my_mind(Particles& particles, myfloat evalx, myfloat evaly, myfloat evalz, myfloat theta2,
  int* max_particles_in_leaf, int* depth);
void stepSimulation(Particles& particles, myfloat dt, myfloat theta2);

#endif