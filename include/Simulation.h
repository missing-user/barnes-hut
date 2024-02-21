#ifndef SIMULATION
#define SIMULATION
#include <functional>
#include "Particle.h"
#include "Bruteforce.h"
#include "Barneshut.h"

void simulate(Particles &particles, myfloat duration, myfloat dt,
              bool brute_force = true, myfloat theta = 0, 
              std::function<void(const Particles&, size_t)> writeCallback = nullptr);
#endif