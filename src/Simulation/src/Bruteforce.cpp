#include "Bruteforce.h"

myvec3 bruteForceAcc(const Particle &p1,
                     const std::vector<Particle> &particles)
{
  // Compute the total acceleration acting on particle p1.
  // the particle itself must be excluded from the computation
  myvec3 acc{0, 0, 0};
  const auto position = p1.p;

  // Tried using simd accumulate here, no measurable difference
  // at -Ofast -march=native. GLM already uses simd instructions
  for (const auto &p2 : particles)
  {
    acc += accelFunc(p2.p - position, p2.m);
  }
  acc -= accelFunc({0, 0, 0}, p1.m); // subtract self

  return acc;
}

std::vector<Particle> stepSimulation(const std::vector<Particle> &particles,
                                     myfloat dt)
{
  // brute force step
  // Computes the state vector of all particles for the next timestep,
  // integrating position and velocity
  std::vector<Particle> particles_next{particles};
#pragma omp parallel for
  for (size_t i = 0; i < particles.size(); i++)
  {
    const auto &p = particles[i];

    auto acc = bruteForceAcc(p, particles);

    particles_next[i].v = p.v + acc * dt;
    particles_next[i].p = p.p + particles_next[i].v * dt;
  }

  return particles_next;
}
