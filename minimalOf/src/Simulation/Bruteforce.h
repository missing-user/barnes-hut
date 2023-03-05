#ifndef BRUTEFORCE
#define BRUTEFORCE

#include "Forces.h"
#include "Particle.h"
#include <vector>

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt);
#endif