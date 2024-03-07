#ifndef ORDER
#define ORDER

#include "Particle.h"
#include "Cuboid.h"
#include <algorithm>
#include <libmorton/morton.h>

template <typename T>
T positionToCode(const myfloat& x, const myfloat& y, const myfloat& z, 
                            const myvec3 &min, const myvec3 &invdimension)
{
  auto xi = (x - min.x) * invdimension.x;
  auto yi = (y - min.y) * invdimension.y;
  auto zi = (z - min.z) * invdimension.z;
  if constexpr(std::is_same<T, uint_fast64_t>::value)
    return libmorton::morton3D_64_encode(xi,yi,zi);
  else
    return libmorton::morton3D_32_encode(xi,yi,zi);
}

template <typename T>
std::vector<T> computeMortonCodes(const Particles &particles, const Cuboid &bb)
{
  std::vector<T> mortonCodes(particles.size());
  // morton order allows for 21 bits per dimension = 63 bits, scale all entries to this size
  // Although libmorton uses unsigned integers, it seemingly expects a range of [-2^20,2^20] for each dimension
  const myvec3 invRange = static_cast<myfloat>(std::pow(2, sizeof(T)*8/3)-1) / bb.dimension;
  assert( (sizeof(T)*8/3) == 21);
  #pragma omp parallel for
  for (int i = 0; i < particles.size(); i++) {
    mortonCodes[i] = positionToCode<T>(particles.p.x[i],particles.p.y[i],particles.p.z[i], bb.min(), invRange);
  }
  return mortonCodes;
}

template <typename T>
std::vector<int> sort_indices(const T& mortonCodes)
{
  std::vector<int> indices(mortonCodes.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&mortonCodes](const int a, const int b) {
    return mortonCodes[a] < mortonCodes[b];
  });
  return indices;
}

void indirect_sort(Particles &particles, const std::vector<int>& indices);
void indirect_sort(std::vector<uint_fast64_t>& mortonCodes, const std::vector<int>& indices);

void reorderByCodes(Particles &particles, const std::vector<uint_fast64_t>& mortonCodes);
void computeAndOrder(Particles &particles, const Cuboid &bb);
#endif
