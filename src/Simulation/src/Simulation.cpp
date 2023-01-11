#include "Simulation.h"
#include <chrono>
#include <string>

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt, double theta) {
  // barnes hut optimized step
  // Buffer for the new state vector of particles
  std::vector<Particle> particles_next{particles};
  std::vector<myvec3> accelerations{particles.size()};

  const Tree mytree(particles);

  // We have constructed the tree, use it to efficiently compute the
  // accelerations. Far away particles get grouped and their contribution is
  // approximated using their center of mass (Barnes Hut algorithm)

  // Use velocity verlet algorithm to update the particle positions and
  // velocities https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
#pragma omp parallel
  {
#pragma omp for
    for (size_t i = 0; i < particles.size(); i++) {
      // First update the positions
      const auto &p1 = particles[i];
      auto &p2 = particles_next[i];

      const auto acc = mytree.computeAcc(p1, theta);
      p2.p = p1.p + p1.v * dt + acc * dt * dt / 2.;
      accelerations[i] = acc;
    }

    // New Octree using the updated particle positions
    // Since the following loop will not change the positions of particles_next,
    // we can safely reuse the vector for iteration and acceleration calculation
    // TODO: Reuse this octree for the next timestep
    const Tree mytree2(particles_next);

#pragma omp for
    for (size_t i = 0; i < particles.size(); i++) {
      // Then update the velocities using v(t+1) = dt*(a(t) + a(t+dt))/2
      auto &p2 = particles_next[i];

      const auto acc2 = mytree2.computeAcc(p2, theta);
      p2.v += (acc2 + accelerations[i]) * dt / 2.;
    }
  }

  return particles_next;
}

std::string makeCsvHeader(size_t numberOfParticles) {
  std::string outputBuffer;

  // Create column names px_0,py_0,pz_0,px_1,py_1,pz_1,...
  // That contain the respective xyz positions of particle 0,1,2,...
  for (size_t i = 0; i < numberOfParticles; i++) {
    std::string istr = std::to_string(i);
    outputBuffer += "px_" + istr + ",py_" + istr + ",pz_" + istr + ",";
  }
  return outputBuffer + "time\n";
}

void simulate(std::vector<Particle> &particles, double duration, myfloat dt,
              std::ostream *outputwriter, bool brute_force, myfloat theta) {
  // The pointer to the outputwriter is optional and will receive the positions
  // of all particles at each timestep if passes
  if (outputwriter != nullptr)
    *outputwriter << makeCsvHeader(particles.size());

  std::cout << "Starting " << duration << "s simulation with " << duration / dt
            << " steps at dt =" << dt << "\n";

  // The timestep must start at 1 or we will simulate one timestep more than
  // necessary
  for (size_t timestep = 1; timestep <= duration / dt; timestep++) {
    if (brute_force)
      particles = stepSimulation(particles, dt);
    else
      particles = stepSimulation(particles, dt, theta);

    if (outputwriter != nullptr)
      *outputwriter << particles << timestep * dt << "\n";
  }

  myfloat residualTimestep = duration - dt * static_cast<size_t>(duration / dt);
  if (residualTimestep != 0.0) {
    std::cout << "completing the residual timestep " << residualTimestep
              << "\n";

    if (brute_force)
      particles = stepSimulation(particles, residualTimestep);
    else
      particles = stepSimulation(particles, residualTimestep, theta);

    if (outputwriter != nullptr)
      *outputwriter << particles << duration << "\n";
  }
}