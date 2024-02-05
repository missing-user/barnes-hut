#include "Order.h"
#include <libmorton/morton.h>
#include <chrono>
#include <parallel/algorithm>

void computeAndOrder(std::vector<Particle> &particles)
{
  // auto time1 = std::chrono::high_resolution_clock::now();

  // Bounds of the universe
  myfloat minx, miny, minz;
  myfloat maxx, maxy, maxz;
  minx = miny = minz = std::numeric_limits<myfloat>::max();
  maxx = maxy = maxz = std::numeric_limits<myfloat>::min();

  // Observed neither speedup nor slowdown with the parallelization
  #pragma omp parallel for reduction(min:minx, miny, minz) reduction(max:maxx, maxy, maxz)
  for (const auto &p : particles)
  {
    minx = std::min(minx, p.p.x);
    miny = std::min(miny, p.p.y);
    minz = std::min(minz, p.p.z);
    maxx = std::max(maxx, p.p.x);
    maxy = std::max(maxy, p.p.y);
    maxz = std::max(maxz, p.p.z);
  }

  // morton order allows for 21 bits per dimension = 63 bits, scale all entries to this size
  myvec3 invRange = static_cast<myfloat>(1<<21) / (myvec3(maxx, maxy, maxz) - myvec3(minx, miny, minz));

  // Morton order
  std::__parallel::sort(particles.begin(), particles.end(), [invRange](const Particle &a, const Particle &b) {
    uint_fast64_t mortonA = libmorton::morton3D_64_encode(a.p.x * invRange.x, a.p.y * invRange.y, a.p.z * invRange.z);
    uint_fast64_t mortonB = libmorton::morton3D_64_encode(b.p.x * invRange.x, b.p.y * invRange.y, b.p.z * invRange.z);
    return mortonA < mortonB;
  });

  // std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - time1;
  // std::cout << "Reordering took " << elapsed.count()<<"s"<<(1<<21) << "\n";
}
