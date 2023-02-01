#ifndef CUBOID
#define CUBOID
// This class is used to store the coords and boundaries of a cuboid in 3D space
// For this class to work, the cuboid's edges must be parallel to the x, y, and
// z
#include "Particle.h"
#include <utility>

class Cuboid
{

public:
  const myvec3 center;     // center point
  const myvec3 dimension;  // length of each dimension
  const myfloat diagonal2; // length squared of the diagonal of the cuboid

  Cuboid(const myvec3& center, const myvec3& dimension);
  Cuboid minMaxCuboid(const myvec3& min, const myvec3& max) const;
  std::array<Cuboid, 8> subdivide() const;            // returns an array of 8 cuboids, splitting the parent
                                                      // cuboid in half along each dimension (i.e. Octant)
  std::array<Cuboid, 8> subdivideAtP(const myvec3& P) const; // returns an array of 8 cuboids, splitting the parent
                                                      // cuboid at the given vector P

  std::string print() const;
};

// A function to calculate the bounding box of a group of particles
template <typename T>
Cuboid bounding_box(const std::vector<T> &particles)
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

  return Cuboid((bmin + bmax) / 2.0, bmax - bmin);
}


class DrawableCuboid{
public:
  myvec3 center;
  myvec3 dimension;
  int level;
  DrawableCuboid(const Cuboid &cuboid, int level);
};

#endif
