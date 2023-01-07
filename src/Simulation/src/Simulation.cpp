#include "Simulation.h"
#include <chrono>
#include <string>

inline myvec3 accelFunc(const myvec3 &diff, myfloat mass) {
  const myfloat softening_param = 0.025;
  return glm::normalize(diff) * mass / (glm::length2(diff) + softening_param);
}

myvec3 getTotalAcceleration(myvec3 position,
                            const std::vector<Particle> &particles) {
  // Compute the total acceleration acting on a particle at position p.
  // the particle itself is excluded from the computation by its particle id

  myvec3 acc{0, 0, 0};

  for (const auto &p2 : particles) {
    if (p2.p != position) {
      const auto diff = p2.p - position;
      acc += accelFunc(diff, p2.m);
    }
  }

  return acc;
}

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt) {
  // brute force step
  // Computes the state vector of all particles for the next timestep,
  // integrating position and velocity
  std::vector<Particle> particles_next{particles};

  for (size_t i = 0; i < particles.size(); i++) {
    const auto &p = particles[i];

    auto acc = getTotalAcceleration(p.p, particles);

    particles_next[i].v = p.v + acc * dt;
    particles_next[i].p = p.p + p.v * dt;
  }

  return particles_next;
}

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt, double theta) {
  // barnes hut optimized step
  // Computes the state vector of all particles for the next timestep,
  // integrating position and velocity
  std::vector<Particle> particles_next{particles};

  // Construct the tree by first computing the bouning box of all particles and
  // then inserting them one by one

  auto begin = std::chrono::steady_clock::now();

  Tree mytree(particles);

  auto end = std::chrono::steady_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
  std::cout << "construct tree: " << elapsed.count() << " ms\n";

  begin = std::chrono::steady_clock::now();
  // Now that we have constructed the tree, use it to efficiently compute the
  // accelerations
  for (size_t i = 0; i < particles.size(); i++) {
    // TODO: iterate in a different order to avoid cache misses
    const auto &p = particles[i];

    myvec3 acc = mytree.computeAcc(p.p, theta);

    particles_next[i].v = p.v + acc * dt;
    particles_next[i].p = p.p + p.v * dt;
  }
  end = std::chrono::steady_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
  std::cout << "calc forces: " << elapsed.count() << " ms\n";

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