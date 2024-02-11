#ifndef SIMULATION
#define SIMULATION
#include <vector>
#include <functional>
#include "Bruteforce.h"
#include "Barneshut.h"

void simulate(std::vector<Particle> &particles, double duration, myfloat dt,
              bool brute_force = true, myfloat theta = 0, 
              std::function<void(const Particles&, size_t)> writeCallback = nullptr);
#endif