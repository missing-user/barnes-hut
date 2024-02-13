#ifndef CUBOID
#define CUBOID
// This class is used to store the coords and boundaries of a cuboid in 3D space
// For this class to work, the cuboid's edges must be parallel to the x, y, and
// z
#include "Particle.h"
#include <utility>
#include <array>
#include <algorithm>


class Cuboid
{ 
public:
  myvec3 center;     // center point
  myvec3 dimension;  // length of each dimension
  myfloat diagonal2; // length squared of the diagonal of the cuboid

  Cuboid(const myvec3& center, const myvec3& dimension);
  std::array<Cuboid, 8> subdivideAtP(const myvec3& P) const; // returns an array of 8 cuboids, splitting the parent
                                                      // cuboid at the given vector P
  std::string print() const;
  myvec3 min() const { return center - dimension*static_cast<myfloat>(0.5); }
  myvec3 max() const { return center + dimension*static_cast<myfloat>(0.5); }
};

// A function to calculate the bounding box of a group of particles
Cuboid bounding_box(const Vectors &positions, const size_t count);
Cuboid minMaxCuboid(const myvec3& min, const myvec3& max); // returns a cuboid with the given min and max coords

class DrawableCuboid{
  // Like a cuboid, but with a depth level field, so that we can draw it colored by depth
public:
  myvec3 center;
  myvec3 dimension;
  int level; // depth level of the cuboid
  DrawableCuboid(const Cuboid &cuboid, int level);
};

#endif
