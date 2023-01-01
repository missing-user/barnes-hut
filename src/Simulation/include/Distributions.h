#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

#include <vector>

#include "Simulation.h"

void set_seed(unsigned int seed);

std::vector<Particle> &set_mass(std::vector<Particle> &particles, myfloat m);

/*******************************************************************/
std::vector<Particle> universe2();

std::vector<Particle> universe3();

std::vector<Particle> universe1();

std::vector<Particle> universe4();

std::vector<Particle> bigbang(int);

std::vector<Particle> stable_orbit();

#endif
