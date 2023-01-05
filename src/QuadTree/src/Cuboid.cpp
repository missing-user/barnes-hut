#include "Cuboid.h"

// This class is used to store the coords and boundaries of a cuboid in 3D space
// For this class to work, the cuboids edges must be parallel to the x, y, and z
Cuboid::Cuboid(myvec3 min, myvec3 max)
    : min_extent(min), max_extent(max)
// center((min_extent + max_extent) / 2.0)
{}
Cuboid::Cuboid(std::pair<myvec3, myvec3> min_max)
    : min_extent(min_max.first), max_extent(min_max.second)
// center((min_extent + max_extent) / 2.0)
{}

bool Cuboid::contains(const Particle &P)
    const // returns whether a particle is insiside cuboid boundary
{
  return glm::all(glm::greaterThanEqual(P.p, min_extent)) &&
         glm::all(glm::lessThan(P.p, max_extent));
}

myvec3 Cuboid::dimension() const { return max_extent - min_extent; }

std::array<Cuboid, 8>
Cuboid::subdivide() const // returns an array of 8 cuboids, splitting the parent
                          // cuboid in half along each dimension (i.e. Octant)
{

  const myfloat x1 = min_extent.x, x2 = max_extent.x, y1 = min_extent.y,
                y2 = max_extent.y, z1 = min_extent.z,
                z2 = max_extent.z; // renaming to facilitate manipulation

  // Calculating the coordinates of the new cuboid divisions
  std::array<Cuboid, 8> subcuboids{
      Cuboid(min_extent, center()),
      Cuboid(myvec3(center().x, y1, z1), myvec3(x2, center().y, center().z)),
      Cuboid(myvec3(x1, center().y, z1), myvec3(center().x, y2, center().z)),
      Cuboid(myvec3(center().x, center().y, z1), myvec3(x2, y2, center().z)),
      Cuboid(myvec3(x1, y1, center().z), myvec3(center().x, center().y, z2)),
      Cuboid(myvec3(center().x, y1, center().z), myvec3(x2, center().y, z2)),
      Cuboid(myvec3(x1, center().y, center().z), myvec3(center().x, y2, z2)),
      Cuboid(center(), max_extent),
  };

  return subcuboids;
}

std::string Cuboid::print() const {
  std::ostringstream str;
  str << min_extent << " " << max_extent;
  std::string s = str.str();
  return s;
}
