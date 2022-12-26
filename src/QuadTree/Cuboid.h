#ifndef CUBOID
#define CUBOID

#include "Particle.h"

class Cuboid {
public:
  const myvec3 min_extent; // minimum corner point
  const myvec3 max_extent; // maximum corner point

  Cuboid(myvec3 min, myvec3 max);
  bool contains(const Particle &P)
      const; // returns whether a particle is insiside cuboid boundary
  std::array<Cuboid, 8>
  subdivide() const; // returns an array of 8 cuboids, splitting the parent
                     // cuboid in half along each dimension (i.e. Octant)
  std::string print() const;
};

#endif
