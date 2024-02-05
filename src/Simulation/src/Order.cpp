#include "Order.h"
#include <libmorton/morton.h>
#include <chrono>
#include <parallel/algorithm>

void computeAndOrder(std::vector<Particle> &particles)
{
  // auto time1 = std::chrono::high_resolution_clock::now();

  // Bounds of the universe
  auto xx = std::minmax_element(particles.begin(), particles.end(), [](const Particle &a, const Particle &b) { return a.p.x < b.p.x; });
  auto yy = std::minmax_element(particles.begin(), particles.end(), [](const Particle &a, const Particle &b) { return a.p.y < b.p.y; });
  auto zz = std::minmax_element(particles.begin(), particles.end(), [](const Particle &a, const Particle &b) { return a.p.z < b.p.z; });

  // morton order allows for 21 bits per dimension = 63 bits, scale all entries to this size
  myvec3 invRange = static_cast<myfloat>(1<<21) / myvec3(xx.second - xx.first, yy.second - yy.first, zz.second - zz.first);

  // Sort using morton order. This is parallelized if _GLIBCXX_PARALLEL is defined
  std::sort(particles.begin(), particles.end(), [invRange](const Particle &a, const Particle &b) {
    uint_fast64_t mortonA = libmorton::morton3D_64_encode(a.p.x * invRange.x, a.p.y * invRange.y, a.p.z * invRange.z);
    uint_fast64_t mortonB = libmorton::morton3D_64_encode(b.p.x * invRange.x, b.p.y * invRange.y, b.p.z * invRange.z);
    return mortonA < mortonB;
  });

  // std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - time1;
  // std::cout << "Reordering took " << elapsed.count()<<"s"<< "\n";
}
