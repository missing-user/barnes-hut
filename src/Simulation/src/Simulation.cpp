#include "Simulation.h"
#include "Progress.h"
#include <iostream>
#include <chrono>
#include <string>
#include "Order.h"

void simulate(std::vector<Particle> &particles, double duration, myfloat dt, bool brute_force, 
              myfloat theta, std::function<void(const Particles&, size_t)> writeCallback) {
  // The pointer to the outputwriter is optional and will receive the
  // positions of all particles at each timestep if passes
  auto theta2 = theta * theta;

  std::cout << "Starting " << duration << "s simulation with " << duration / dt
            << " steps at dt=" << dt << "\nwith the "<< (brute_force ? "brute force" : "Barnes Hut") 
            << " solver and n="<< particles.size() << " particles\n";

  boost::timer::progress_display show_progress(duration / dt);

  // The timestep must start at 1 or we will simulate one timestep more than
  // necessary

  Particles bodies{particles.size()};
#pragma omp parallel for // First touch initialization
  for (size_t i = 0; i < bodies.size(); i++) {
    bodies.p[i] = particles[i].p;
    bodies.v[i] = particles[i].v;
    bodies.m[i] = particles[i].m;
  }


  if (writeCallback){
    writeCallback(bodies, 0); // Write initial state
  }
  for (size_t timestep = 1; timestep <= duration / dt; timestep++) {

    if (brute_force) // Brute force step
      stepSimulation(bodies, dt);
    else // Barnes Hut step
      stepSimulation(bodies, dt, theta2);

    if (writeCallback){
      writeCallback(bodies, timestep);
    }

    ++show_progress;
  }

  myfloat residualTimestep = duration - dt * static_cast<size_t>(duration / dt);
  if (residualTimestep != 0.0) {
    if (brute_force)
      stepSimulation(bodies, residualTimestep);
    else
      stepSimulation(bodies, residualTimestep, theta2);
  }
}