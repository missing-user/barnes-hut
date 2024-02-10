#ifndef ORDER
#define ORDER

#include "Particle.h"
#include "Cuboid.h"

void computeAndOrder(std::vector<Particle> &particles);
void computeAndOrder(Particles &particles, Cuboid bb);
#endif
