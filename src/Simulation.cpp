#include "Simulation.h"
#include "Progress.h"
#include <iostream>
#include <chrono>
#include <string>
#include "Order.h"

void simulate(Particles &particles, myfloat duration, myfloat dt, bool brute_force, 
              myfloat theta, std::function<void(const Particles&, size_t)> writeCallback) {
  // The pointer to the outputwriter is optional and will receive the
  // positions of all particles at each timestep if passes
  auto theta2 = theta * theta;

  std::cout << "Starting " << duration << "s simulation with " << duration / dt
            << " steps at dt=" << dt << "\nwith the "<< (brute_force ? "brute force" : "Barnes Hut") 
            << " solver and n="<< particles.size() << " particles\n";

  boost::timer::progress_display show_progress(duration / dt);

  if (writeCallback){
    writeCallback(particles, 0); // Write initial state
  }
  // The timestep must start at 1 or we will simulate one timestep more than
  // necessary
  for (size_t timestep = 1; timestep <= duration / dt; timestep++) {

    if (brute_force) // Brute force step
      stepSimulation(particles, dt);
    else // Barnes Hut step
      stepSimulation(particles, dt, theta2);

    if (writeCallback){
      writeCallback(particles, timestep);
    }

    ++show_progress;
  }

  myfloat residualTimestep = duration - dt * static_cast<size_t>(duration / dt);
  if (residualTimestep != 0.0) {
    if (brute_force)
      stepSimulation(particles, residualTimestep);
    else
      stepSimulation(particles, residualTimestep, theta2);
  }
}