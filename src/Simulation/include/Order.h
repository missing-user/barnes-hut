#ifndef ORDER
#define ORDER

#include "Particle.h"
#include "Cuboid.h"

std::vector<uint_fast64_t> computeMortonCodes(const Particles &particles, const Cuboid &bb);
void reorderByCodes(Particles &particles, const std::vector<uint_fast64_t>& mortonCodes);
void computeAndOrder(Particles &particles, const Cuboid &bb);
#endif
