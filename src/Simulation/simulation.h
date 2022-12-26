#ifndef SIMULATION
#define SIMULATION

#include <vector>
#include <string>
#include <iostream>

#include "../QuadTree/QuadTree.h"

#include <glm/glm.hpp>    // pvec3
#include <glm/gtx/io.hpp> // Allows us to easily std::cout << pvec3;

myvec3 getTotalAcceleration(myvec3 position, const std::vector<Particle> &particles);
std::vector<Particle> stepSimulation(const std::vector<Particle> &particles, myfloat dt);
std::vector<Particle> stepSimulation(const std::vector<Particle> &particles, myfloat dt, double theta);
std::string makeCsvHeader(size_t numberOfParticles);
void simulate(std::vector<Particle> &particles, double duration, myfloat dt, std::ostream *outputwriter = nullptr, bool brute_force = true, myfloat theta = 0);

#endif