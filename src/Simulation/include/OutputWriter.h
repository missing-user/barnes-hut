#include "Particle.h"

void writeToCsvFile(const std::vector<Particle> &particles,
                    const std::size_t timestep);
void writeToBinaryFile(const std::vector<Particle> &particles,
                    const std::size_t timestep);