#ifndef SIMULATION
#define SIMULATION

#include <iostream>
#include <vector>

#include "Bruteforce.h"
#include "Tree.h"

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt, double theta);
void simulate(std::vector<Particle> &particles, double duration, myfloat dt,
              bool outputwriter = false, bool brute_force = true,
              myfloat theta = 0);
#endif