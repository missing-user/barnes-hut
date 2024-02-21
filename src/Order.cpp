#include "Order.h"
#include <libmorton/morton.h>
#include <chrono>
#include <algorithm>
#include <numeric>

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
  // morton order allows for 21 bits per dimension = 63 bits, scale all entries to this size
  // Although libmorton uses unsigned integers, it seemingly expects a range of [-2^20,2^20] for each dimension
  const myvec3 invRange = static_cast<myfloat>(std::pow(2, 21)-1) / bb.dimension;
  #pragma omp parallel for
  for (int i = 0; i < particles.size(); i++) {
    mortonCodes[i] = positionToCode(particles.p.x[i],particles.p.y[i],particles.p.z[i], bb.min(), invRange);
  }
  return mortonCodes;
}

void reorderByCodes(Particles &particles, const std::vector<uint_fast64_t>& mortonCodes){
  // Sort using morton order. This is parallelized if _GLIBCXX_PARALLEL is defined
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<int> indices(particles.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&mortonCodes](const int a, const int b) {
    return mortonCodes[a] < mortonCodes[b];
  });

  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time to sort: " << elapsed.count() << "ms\n";
  start = std::chrono::high_resolution_clock::now();

  /* Reordering the particles component by component is faster than having a merged loop
  *  in which we copy all 7 components of the particles struct at once. Especially for 
  *  large particle counts, only allocating the temporary vector for a single component,
  *  and copying the reordered entries, is significantly more efficient than allocating 
  *  a temporary Particles struct (7 vectors), copying all components, then swapping all pointers.
  */
  vector_type reorderedComponent(particles.size());
  // Iterate over all member variables of the Particles struct (x, y, z, vx, vy, vz, m)
  for (int member_var = 0; member_var < 7; member_var++)
  {
    #pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
      reorderedComponent[i] = particles.get(member_var)[indices[i]];
    }
    // Swap pointers of the reordered particles with the original particles for the next iteration
    std::swap(particles.get(member_var), reorderedComponent);
  }
  
  end = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time to reorder: " << elapsed.count() << " ms\n";
}

void computeAndOrder(Particles &particles, const Cuboid &bb)
{
  std::vector<uint_fast64_t> mortonCodes = computeMortonCodes(particles, bb);
  reorderByCodes(particles, mortonCodes);
}
