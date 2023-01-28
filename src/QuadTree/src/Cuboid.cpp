#include "Cuboid.h"

// This class is used to store the coords and boundaries of a cuboid in 3D space
// For this class to work, the cuboids edges must be parallel to the x, y, and z
Cuboid::Cuboid(myvec3 center, myvec3 dimension)
    : center(center), dimension(dimension), diagonal2(glm::length2(dimension)) {
}

std::array<Cuboid, 8>
Cuboid::subdivide() const // returns an array of 8 cuboids, splitting the parent
                          // cuboid in half along each dimension (i.e. Octant)
{
  const myvec3 newDimension = dimension / 2.0;
  // Calculating the coordinates of the new cuboid divisions
  std::array<Cuboid, 8> subcuboids{
      Cuboid(center + myvec3(-.5, -.5, -.5) * newDimension, newDimension),
      Cuboid(center + myvec3(.5, -.5, -.5) * newDimension, newDimension),
      Cuboid(center + myvec3(-.5, .5, -.5) * newDimension, newDimension),
      Cuboid(center + myvec3(.5, .5, -.5) * newDimension, newDimension),
      Cuboid(center + myvec3(-.5, -.5, .5) * newDimension, newDimension),
      Cuboid(center + myvec3(.5, -.5, .5) * newDimension, newDimension),
      Cuboid(center + myvec3(-.5, .5, .5) * newDimension, newDimension),
      Cuboid(center + myvec3(.5, .5, .5) * newDimension, newDimension),
  };

  return subcuboids;
}

std::string Cuboid::print() const {
  std::ostringstream str;
  str << center << " diagonal: " << diagonal2;
  std::string s = str.str();
  return s;
}