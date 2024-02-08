#include "Order.h"
#include <libmorton/morton.h>
#include <chrono>
#include <algorithm>
#include <numeric>

void computeAndOrder(std::vector<Particle> &particles)
{
  // auto time1 = std::chrono::high_resolution_clock::now();

  // // Bounds of the universe
  // auto xx = std::minmax_element(particles.begin(), particles.end(), [](const Particle &a, const Particle &b) { return a.p.x < b.p.x; });
  // auto yy = std::minmax_element(particles.begin(), particles.end(), [](const Particle &a, const Particle &b) { return a.p.y < b.p.y; });
  // auto zz = std::minmax_element(particles.begin(), particles.end(), [](const Particle &a, const Particle &b) { return a.p.z < b.p.z; });

  // // morton order allows for 21 bits per dimension = 63 bits, scale all entries to this size
  // myvec3 dimension = myvec3(xx.second->p.x - xx.first->p.x, yy.second->p.y - yy.first->p.y, zz.second->p.z - zz.first->p.z);
  // myvec3 invRange = static_cast<myfloat>(1<<21) / dimension;
  // // Sort using morton order. This is parallelized if _GLIBCXX_PARALLEL is defined
  // std::sort(particles.begin(), particles.end(), [invRange, xx,yy,zz](const Particle &a, const Particle &b) {
  //   uint_fast64_t mortonA = libmorton::morton3D_64_encode(a.p.x * invRange.x, a.p.y * invRange.y, a.p.z * invRange.z);
  //   uint_fast64_t mortonB = libmorton::morton3D_64_encode(b.p.x * invRange.x, b.p.y * invRange.y, b.p.z * invRange.z);
  //   // uint_fast64_t mortonA = libmorton::morton3D_64_encode((a.p.x - xx.first->p.x) * invRange.x, (a.p.y - yy.first->p.y) * invRange.y, (a.p.z - zz.first->p.z) * invRange.z);
  //   // uint_fast64_t mortonB = libmorton::morton3D_64_encode((b.p.x - xx.first->p.x) * invRange.x, (b.p.y - yy.first->p.y) * invRange.y, (b.p.z - zz.first->p.z) * invRange.z);
  //   return mortonA < mortonB;
  // });

  // std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - time1;
  // std::cout << "Reordering particles took " << elapsed.count()<<"s "<<invRange<<" size "<< dimension<< "\n";

  
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


void computeAndOrder(Particles &particles, Cuboid bb)
{
  auto time1 = std::chrono::high_resolution_clock::now();

  std::vector<uint_fast64_t> mortonCodes(particles.size());
  myvec3 invRange = (std::pow(2, 21)-1) / bb.dimension;
  #pragma omp parallel for
  for (size_t i = 0; i < particles.size(); i++) {
    uint_fast32_t x = (particles.p.x[i] - bb.min().x) * invRange.x;
    uint_fast32_t y = (particles.p.y[i] - bb.min().y) * invRange.y;
    uint_fast32_t z = (particles.p.z[i] - bb.min().z) * invRange.z;
    mortonCodes[i] = libmorton::morton3D_64_encode(x,y,z);
  }

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
  std::swap(particles.p.x, reorderedParticles.p.x);
  std::swap(particles.p.y, reorderedParticles.p.y);
  std::swap(particles.p.z, reorderedParticles.p.z);
  std::swap(particles.v.x, reorderedParticles.v.x);
  std::swap(particles.v.y, reorderedParticles.v.y);
  std::swap(particles.v.z, reorderedParticles.v.z);
  std::swap(particles.m, reorderedParticles.m);

  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Index reordering took " << elapsed.count()<<"s"<< "\n";
}
