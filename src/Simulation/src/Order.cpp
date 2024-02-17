#include "Order.h"
#include <libmorton/morton.h>
#include <chrono>
#include <algorithm>
#include <numeric>

void computeAndOrder(std::vector<Particle> &particles)
{ 
  Particles bodies{particles.size()};

  for (size_t i = 0; i < bodies.size(); i++) {
    bodies.p[i] = particles[i].p;
    bodies.v[i] = particles[i].v;
    bodies.m[i] = particles[i].m;
  }

  computeAndOrder(bodies, bounding_box(bodies.p, bodies.size()));

  for (size_t i = 0; i < bodies.size(); i++) {
    particles[i].p.x = bodies.p.x[i];
    particles[i].p.y = bodies.p.y[i];
    particles[i].p.z = bodies.p.z[i];
    particles[i].v.x = bodies.v.x[i];
    particles[i].v.y = bodies.v.y[i];
    particles[i].v.z = bodies.v.z[i];
    particles[i].m = bodies.m[i];
  }
}

uint_fast64_t positionToCode(const myfloat& x, const myfloat& y, const myfloat& z, 
                            const myvec3 &min, const myvec3 &invdimension)
{
  uint_fast32_t xi = (x - min.x) * invdimension.x;
  uint_fast32_t yi = (y - min.y) * invdimension.y;
  uint_fast32_t zi = (z - min.z) * invdimension.z;
  return libmorton::morton3D_64_encode(xi,yi,zi);
}

std::vector<uint_fast64_t> computeMortonCodes(const Particles &particles, const Cuboid &bb)
{
  std::vector<uint_fast64_t> mortonCodes(particles.size());
  const myvec3 invRange = static_cast<myfloat>(std::pow(2, 21)-1) / bb.dimension;
  #pragma omp parallel for
  for (size_t i = 0; i < particles.size(); i++) {
    mortonCodes[i] = positionToCode(particles.p.x[i],particles.p.y[i],particles.p.z[i], bb.min(), invRange);
  }
  return mortonCodes;
}

void reorderByCodes(Particles &particles, const std::vector<uint_fast64_t>& mortonCodes){
  // morton order allows for 21 bits per dimension = 63 bits, scale all entries to this size
  // Although libmorton uses unsigned integers, it seemingly expects a range of [-2^20,2^20] for each dimension
  // Sort using morton order. This is parallelized if _GLIBCXX_PARALLEL is defined
  std::vector<int> indices(particles.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&mortonCodes](const int a, const int b) {
    return mortonCodes[a] < mortonCodes[b];
  });

  // Reorder the particles
  Particles reorderedParticles{particles.size()};
  #pragma omp parallel for
  for (size_t i = 0; i < particles.size(); i++) {
    reorderedParticles.p.x[i] = particles.p.x[indices[i]];
    reorderedParticles.p.y[i] = particles.p.y[indices[i]];
    reorderedParticles.p.z[i] = particles.p.z[indices[i]];
    reorderedParticles.v.x[i] = particles.v.x[indices[i]];
    reorderedParticles.v.y[i] = particles.v.y[indices[i]];
    reorderedParticles.v.z[i] = particles.v.z[indices[i]];
    reorderedParticles.m[i] = particles.m[indices[i]];
  }
  // Swap pointers of the reordered particles with the original particles
  std::swap(particles.p.x, reorderedParticles.p.x);
  std::swap(particles.p.y, reorderedParticles.p.y);
  std::swap(particles.p.z, reorderedParticles.p.z);
  std::swap(particles.v.x, reorderedParticles.v.x);
  std::swap(particles.v.y, reorderedParticles.v.y);
  std::swap(particles.v.z, reorderedParticles.v.z);
  std::swap(particles.m, reorderedParticles.m);
}

void computeAndOrder(Particles &particles, const Cuboid &bb)
{
  std::vector<uint_fast64_t> mortonCodes = computeMortonCodes(particles, bb);
  reorderByCodes(particles, mortonCodes);
}
