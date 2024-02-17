#ifndef BARNESHUT
#define BARNESHUT

#include "Cuboid.h"
#include "Particle.h"
#include "Forces.h"
std::vector<DrawableCuboid> stepSimulation(Particles& particles, myfloat dt, myfloat theta2);
#endif