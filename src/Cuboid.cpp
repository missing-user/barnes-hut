#include "Cuboid.h"

// This class is used to store the coords and boundaries of a cuboid in 3D space
// For this class to work, the cuboids edges must be parallel to the x, y, and z
Cuboid::Cuboid(const myvec3& center, const myvec3& dimension)
    : center(center), dimension(dimension), diagonal2(glm::length2(dimension))
{}

Cuboid minMaxCuboid(const myvec3& min, const myvec3& max) // returns a cuboid with the given min and max coords
{
  return Cuboid(static_cast<myfloat>(0.5) * (min + max), max - min);
}

std::array<Cuboid, 8>
Cuboid::subdivideAtP(const myvec3& P) const // returns an array of 8 cuboids, splitting the parent
                                     // cuboid in half along each dimension (i.e. Octant)
{
  myvec3 minP = center - dimension/static_cast<myfloat>(2.0);
  myvec3 maxP = center + dimension/static_cast<myfloat>(2.0);
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

Cuboid bounding_box(const Vectors &positions)
{
  std::pair<myfloat, myfloat> xx{std::numeric_limits<myfloat>::max(), std::numeric_limits<myfloat>::min()};
  std::pair<myfloat, myfloat> yy{std::numeric_limits<myfloat>::max(), std::numeric_limits<myfloat>::min()};
  std::pair<myfloat, myfloat> zz{std::numeric_limits<myfloat>::max(), std::numeric_limits<myfloat>::min()};
  // This already operates in the memory bound regime, no gain from parallelization
  #pragma omp simd
  for (size_t i = 0; i < positions.size(); i++)
  {
    xx.first =  std::min(xx.first,  positions.x[i]);
    xx.second = std::max(xx.second, positions.x[i]);
    yy.first =  std::min(yy.first,  positions.y[i]);
    yy.second = std::max(yy.second, positions.y[i]);
    zz.first =  std::min(zz.first,  positions.z[i]);
    zz.second = std::max(zz.second, positions.z[i]);
  }
  
  return minMaxCuboid(myvec3(xx.first, yy.first, zz.first), myvec3(xx.second, yy.second, zz.second));
}

std::string Cuboid::print() const
{
  std::ostringstream str;
  str << center << " diagonal: " << diagonal2;
  std::string s = str.str();
  return s;
}


DrawableCuboid::DrawableCuboid(const Cuboid &cuboid, int level) : 
  center(cuboid.center), dimension(cuboid.dimension), level(level), isLeaf(false) {}

DrawableCuboid::DrawableCuboid(const myvec3 &pos, int level):center(pos), level(level), isLeaf(true){}
