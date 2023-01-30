#include "Cuboid.h"

// This class is used to store the coords and boundaries of a cuboid in 3D space
// For this class to work, the cuboids edges must be parallel to the x, y, and z
Cuboid::Cuboid(const myvec3& center, const myvec3& dimension)
    : center(center), dimension(dimension), diagonal2(glm::length2(dimension))
{
}

Cuboid Cuboid::minMaxCuboid(myvec3 min, myvec3 max) const
{
  return Cuboid(0.5 * (min + max), max - min);
}

std::array<Cuboid, 8> Cuboid::subdivide() const // returns an array of 8 cuboids, splitting the parent
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

std::array<Cuboid, 8>
Cuboid::subdivideAtP(myvec3 P) const // returns an array of 8 cuboids, splitting the parent
                                     // cuboid in half along each dimension (i.e. Octant)
{
  myvec3 minP = center - dimension/2.0;
  myvec3 maxP = center + dimension/2.0;
  // Calculating the coordinates of the new cuboid divisions
  std::array<Cuboid, 8> subcuboids{
      minMaxCuboid(myvec3(minP.x, minP.y, minP.z), myvec3(P.x, P.y, P.z)),
      minMaxCuboid(myvec3(P.x, minP.y, minP.z), myvec3(maxP.x, P.y, P.z)),
      minMaxCuboid(myvec3(minP.x, P.y, minP.z), myvec3(P.x, maxP.y, P.z)),
      minMaxCuboid(myvec3(P.x, P.y, minP.z), myvec3(maxP.x, maxP.y, P.z)),
      minMaxCuboid(myvec3(minP.x, minP.y, P.z), myvec3(P.x, P.y, maxP.z)),
      minMaxCuboid(myvec3(P.x, minP.y, P.z), myvec3(maxP.x, P.y, maxP.z)),
      minMaxCuboid(myvec3(minP.x, P.y, P.z), myvec3(P.x, maxP.y, maxP.z)),
      minMaxCuboid(myvec3(P.x, P.y, P.z), myvec3(maxP.x, maxP.y, maxP.z)),
  };

  return subcuboids;
}

std::string Cuboid::print() const
{
  std::ostringstream str;
  str << center << " diagonal: " << diagonal2;
  std::string s = str.str();
  return s;
}


DrawableCuboid::DrawableCuboid(const Cuboid &cuboid, int level) : 
  center(cuboid.center), dimension(cuboid.dimension), level(level) {}
