#include "Simulation.h"
#include "Progress.h"
#include <chrono>
#include <string>

std::vector<Particle> stepSimulation(const std::vector<Particle>& particles,
                                     myfloat dt, double theta2) {
  // barnes hut optimized step
  // Buffer for the new state vector of particles
  std::vector<Particle> particles_next{particles};
  std::vector<myvec3> accelerations{particles.size()};
  const Tree mytree(particles);

  // We have constructed the tree, use it to efficiently compute the
  // accelerations. Far away particles get grouped and their contribution is
  // approximated using their center of mass (Barnes Hut algorithm)

  // Use velocity verlet algorithm to update the particle positions and
  // velocities
  // https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < particles.size(); i++) {
    // First update the positions
    const auto &p1 = particles[i];
    auto &p2 = particles_next[i];

    const auto acc = mytree.computeAcc(p1, theta2);
    p2.p = p1.p + p1.v * dt + acc * dt * dt / 2.;
    accelerations[i] = acc;
  }

  // New Octree using the updated particle positions
  // Since the following loop will not change the positions of particles_next,
  // we can safely reuse the vector for iteration and acceleration calculation
  // TODO: Reuse this octree for the next timestep

  const Tree mytree2(particles_next);
 
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < particles.size(); i++) {
    // Then update the velocities using v(t+1) = dt*(a(t) + a(t+dt))/2
    auto &p2 = particles_next[i];

    const auto acc2 = mytree2.computeAcc(p2, theta2);
    p2.v += (acc2 + accelerations[i]) * dt / 2.;
  }

  return particles_next;
}

void simulate(std::vector<Particle> &particles, double duration, myfloat dt, bool brute_force, 
              myfloat theta, std::function<void(const std::vector<Particle>&, size_t)> writeCallback) {
  // The pointer to the outputwriter is optional and will receive the
  // positions of all particles at each timestep if passes
  auto theta2 = theta*theta; // theta is squared in the algorithm

  std::cout << "Starting " << duration << "s simulation with " << duration / dt
            << " steps at dt=" << dt << "\nwith the "<< (brute_force ? "brute force" : "Barnes Hut") 
            << " solver and n="<< particles.size() << " particles\n";

  boost::timer::progress_display show_progress(duration / dt);

  // The timestep must start at 1 or we will simulate one timestep more than
  // necessary
  for (size_t timestep = 1; timestep <= duration / dt; timestep++) {

    if (brute_force) // Brute force step
      particles = stepSimulation(particles, dt);
    else // Barnes Hut step
      particles = stepSimulation(particles, dt, theta2);

    if (writeCallback){
      writeCallback(particles, timestep);
    }

    ++show_progress;
  }

  myfloat residualTimestep = duration - dt * static_cast<size_t>(duration / dt);
  if (residualTimestep != 0.0) {
    if (brute_force)
      particles = stepSimulation(particles, residualTimestep);
    else
      particles = stepSimulation(particles, residualTimestep, theta2);
  }
}