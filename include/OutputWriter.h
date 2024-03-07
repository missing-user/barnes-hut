#include "Particle.h"

void writeToCsvFile(const Particles &particles,
                    const std::size_t timestep);

void writeToBinaryFile(const Particles &particles,
                    const std::size_t timestep);