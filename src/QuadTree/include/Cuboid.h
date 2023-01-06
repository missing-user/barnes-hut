#ifndef CUBOID
#define CUBOID
// This class is used to store the coords and boundaries of a cuboid in 3D space
// For this class to work, the cuboid's edges must be parallel to the x, y, and
// z
#include "Particle.h"
#include <utility>

class Cuboid {

public:
  const myvec3 center;     // center point
  const myvec3 dimension;  // length of each dimension
  const myfloat diagonal2; // length squared of the diagonal of the cuboid

  Cuboid(myvec3 center, myvec3 dimension);

  std::array<Cuboid, 8>
  subdivide() const; // returns an array of 8 cuboids, splitting the parent
                     // cuboid in half along each dimension (i.e. Octant)
  std::string print() const;
};

Cuboid bounding_box(const std::vector<Particle> &);

#endif
