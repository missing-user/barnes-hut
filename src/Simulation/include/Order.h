#ifndef ORDER
#define ORDER

#include "Particle.h"
#include "Cuboid.h"

std::vector<uint_fast64_t> computeMortonCodes(const Particles &particles, Cuboid bb);
void computeAndOrder(std::vector<Particle> &particles);
void reorderByCodes(Particles &particles,std::vector<uint_fast64_t>& mortonCodes);
void computeAndOrder(Particles &particles, Cuboid bb);
#endif
