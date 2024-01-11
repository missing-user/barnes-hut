#ifndef SIMULATION
#define SIMULATION

#include <iostream>
#include <vector>
#include <functional>

#include "Bruteforce.h"
#include "Tree.h"

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt, double theta);
void simulate(std::vector<Particle> &particles, double duration, myfloat dt,
              bool brute_force = true, myfloat theta = 0, 
              std::function<void(const std::vector<Particle>&, size_t)> writeCallback = nullptr);
#endif