#ifndef ORDER
#define ORDER

#include "Particle.h"
#include "Tree.h"

void reorder(std::vector<Particle> &data,
             std::vector<std::size_t> const &order);

void computeAndOrder(std::vector<Particle> &particles);
#endif
