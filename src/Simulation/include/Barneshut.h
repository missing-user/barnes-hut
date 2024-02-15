#ifndef BARNESHUT
#define BARNESHUT

#include "Cuboid.h"
#include "Particle.h"
#include <atomic>
std::vector<DrawableCuboid> stepSimulation(Particles& particles, myfloat dt, myfloat theta);
extern std::atomic<int> force_count;

#endif