#include "Bruteforce.h"

myvec3 bruteForceAcc(const Particle &p1,
                     const std::vector<Particle> &particles) {
  myvec3 acc{0, 0, 0};
  const auto position = p1.p;

  /* Tried manual vectorization using simd accumulate here, no measurable
   * difference at -Ofast -march=native. GLM already uses vectorized simd
   * instructions for the accumulation.
   */
  for (const auto &p2 : particles) {
    if (position == p2.p)
      continue;
    acc += accelFunc(p2.p - position, p2.m);
  }

  /* subtract self contribution at the end. This is in some cases more
   * performant than checking if the interacting particle is itself at every
   * iteration, seems to be hardware dependent
   */
  // acc -= accelFunc({0, 0, 0}, p1.m);

  return acc;
}

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt) {

  /* Computes the state vector of all particles for the next timestep,
   * integrating position and velocity by summing up all the individual particle
   * contributions (n^2 brute force algorithm)
   */
  std::vector<Particle> particles_next{particles};
  std::vector<myvec3> accelerations{particles.size()};

/*Multithreading the outer loop instead of the inner loop is more than 10x
 * faster. This is because the inner loop is very short and the overhead of
 * creating threads is too high.
 */
#pragma omp parallel
  {
#pragma omp for
    for (size_t i = 0; i < particles.size(); i++) {
      // First update the positions
      const auto &p1 = particles[i];
      auto &p2 = particles_next[i];

      const auto acc = bruteForceAcc(p1, particles);
      p2.p = p1.p + p1.v * dt + acc * dt * dt / 2.;
      accelerations[i] = acc;
    }

#pragma omp for
    for (size_t i = 0; i < particles.size(); i++) {
      // Then update the velocities using v(t+1) = dt*(a(t) + a(t+dt))/2
      auto &p2 = particles_next[i];

      const auto acc2 = bruteForceAcc(p2, particles_next);
      p2.v += (acc2 + accelerations[i]) * dt / 2.;
    }
  }

  return particles_next;
}
