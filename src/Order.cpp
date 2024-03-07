#include "Order.h"
#include <algorithm>
#include <numeric>

#ifdef MEASURE_TIME
#include <chrono>
#endif

void indirect_sort(Particles &particles, const std::vector<int>& indices){
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
}

void indirect_sort(std::vector<uint_fast64_t>& mortonCodes, const std::vector<int>& indices){
  std::vector<uint_fast64_t> reorderedComponent(mortonCodes.size());
  #pragma omp parallel for
  for (int i = 0; i < mortonCodes.size(); i++) {
    reorderedComponent[i] = mortonCodes[indices[i]];
  }
  std::swap(mortonCodes, reorderedComponent);
}

void reorderByCodes(Particles &particles, const std::vector<uint_fast64_t>& mortonCodes){
  // Sort using morton order. This is parallelized if _GLIBCXX_PARALLEL is defined
  
  #ifdef MEASURE_TIME
  auto start = std::chrono::high_resolution_clock::now();
  #endif
  std::vector<int> indices = sort_indices(mortonCodes);

  #ifdef MEASURE_TIME
  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time to sort: " << elapsed.count() << "ms\n";
  start = std::chrono::high_resolution_clock::now();
  #endif
  indirect_sort(particles, indices);
  
  #ifdef MEASURE_TIME
  end = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time to reorder: " << elapsed.count() << " ms\n";
  #endif
}

void computeAndOrder(Particles &particles, const Cuboid &bb)
{
  std::vector<uint_fast64_t> mortonCodes = computeMortonCodes<uint_fast64_t>(particles, bb);
  reorderByCodes(particles, mortonCodes);
}
