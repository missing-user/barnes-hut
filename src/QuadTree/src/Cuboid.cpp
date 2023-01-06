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

// A function to calculate the bounding box of a group of particles
Cuboid bounding_box(const std::vector<Particle> &particles) {
  myvec3 bmin = particles.at(0).p;
  myvec3 bmax = particles.at(0).p;

  for (const auto &p : particles) {
    bmin = glm::min(bmin, p.p);
    bmax = glm::max(bmax, p.p);
  }
  // Bounding box containment is defined as a half open interval: [min, max)
  // To actually make the particles be contained within bmax, we need to add a
  // small offset
  bmax += myvec3{1e-6, 1e-6, 1e-6};

  return Cuboid((bmin + bmax) / 2.0, bmax - bmin);
}