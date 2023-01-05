#ifndef CUBOID
#define CUBOID
// This class is used to store the coords and boundaries of a cuboid in 3D space
// For this class to work, the cuboid's edges must be parallel to the x, y, and
// z
#include "Particle.h"
#include <utility>

class Cuboid {
  const myvec3 min_extent; // minimum corner point
  const myvec3 max_extent; // maximum corner point
public:
  Cuboid(myvec3 min, myvec3 max);
  Cuboid(std::pair<myvec3, myvec3> min_max);
  myvec3 dimension() const;
  inline myvec3 center() const {
    return (min_extent + max_extent) / static_cast<myfloat>(2);
  }
  bool contains(const Particle &P)
      const; // returns whether a particle is insiside cuboid boundary
  std::array<Cuboid, 8>
  subdivide() const; // returns an array of 8 cuboids, splitting the parent
                     // cuboid in half along each dimension (i.e. Octant)
  std::string print() const;
};

#endif
