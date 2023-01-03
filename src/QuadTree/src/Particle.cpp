
#include "Particle.h"
// class to hold the position and velocity of a particle

// A function to calculate the bounding box of a group of particles
std::pair<myvec3, myvec3> bounding_box(const std::vector<Particle> &particles)
{
  myvec3 bmin = particles.at(0).p;
  myvec3 bmax = particles.at(0).p;

  for (const auto &p : particles)
  {
    bmin = glm::min(bmin, p.p);
    bmax = glm::max(bmax, p.p);
  }
  // Bounding box containment is defined as a half open interval: [min, max)
  // To actually make the particles be contained within bmax, we need to add a
  // small offset
  bmax += myvec3{1e-6, 1e-6, 1e-6};

  return {bmin, bmax};
}

Particle operator+(const Particle &P1, const Particle &P2)
{ // Overload the + operator to add two particles position together
  // This is used to calculate the center of mass of a group of particles
  // The velocity of the center of mass is not calculated
  Particle p;
  p.p = P1.p + P2.p;
  p.m = P1.m + P2.m;
  return p;
}