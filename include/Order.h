#ifndef ORDER
#define ORDER

#include "Particle.h"
#include "Cuboid.h"

std::vector<uint_fast64_t> computeMortonCodes(const Particles &particles, const Cuboid &bb);



void indirect_sort(Particles &particles, const std::vector<int>& indices);
void indirect_sort(std::vector<uint_fast64_t>& mortonCodes, const std::vector<int>& indices);
std::vector<int> sort_indices(const std::vector<uint_fast64_t>& mortonCodes);

void reorderByCodes(Particles &particles, const std::vector<uint_fast64_t>& mortonCodes);
void computeAndOrder(Particles &particles, const Cuboid &bb);
#endif
