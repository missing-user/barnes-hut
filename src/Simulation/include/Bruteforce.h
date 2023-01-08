#ifndef BRUTEFORCE
#define BRUTEFORCE

#include "Forces.h"
#include "Particle.h"
#include <vector>
myvec3 bruteForceAcc(const Particle &p1,
                     const std::vector<Particle> &particles);

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt);
#endif