#ifndef SIMULATION
#define SIMULATION

#include <iostream>
#include <vector>

#include "Bruteforce.h"
#include "QuadTree.h"

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt, double theta);
std::string makeCsvHeader(size_t numberOfParticles);
void simulate(std::vector<Particle> &particles, double duration, myfloat dt,
              std::ostream *outputwriter = nullptr, bool brute_force = true,
              myfloat theta = 0);

#endif