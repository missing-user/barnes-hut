
#include "Particle.h"

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
  // To actually make the particles contained within bmax, we need to add a small offset
  bmax += myvec3{1e-6, 1e-6, 1e-6};

  return {bmin, bmax};
}

Particle operator+(const Particle &P1, const Particle &P2)
{
  Particle p;
  p.p = P1.p + P2.p;
  // p.v = P1.v + P2.v; // Only used for center of mass calc, no v
  p.m = P1.m + P2.m;
  return p;
}