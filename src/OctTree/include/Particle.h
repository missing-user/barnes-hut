#ifndef PARTICLE
#define PARTICLE

#include <sstream>
#include <vector>
#include <iostream>
#define GLM_FORCE_INTRINSICS 
//#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES // Slowdowns on my machine
//#define GLM_FORCE_INLINE // Didn't make a measurable difference
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp> // Allows us to easily std::cout << pvec3;
#include "xsimd/xsimd.hpp"
namespace xs = xsimd;

typedef glm::dvec3 myvec3; // We can easily switch the entire implementation to
                           // float precision by adjusting these two variables
typedef double myfloat;
using vector_type = std::vector<myfloat, xsimd::aligned_allocator<myfloat>>;

struct Particle { // A particle with position, velocity and unique id
  myvec3 p;
  myvec3 v;
  myfloat m;
};

template <typename T>
inline T length2(T x, T y, T z) {
  return x * x + y * y + z * z;
}

//template <typename T>
struct Vectors{
  vector_type x, y, z;
  size_t size() const noexcept {
    // assert(x.size() == y.size());
    // assert(x.size() == z.size());
    return x.size();
  }
  Vectors(size_t size) : x(vector_type(size)), y(vector_type(size)), z(vector_type(size)) {}
};

class Particles{
public:
  Vectors p, v;
  vector_type m;
public:
  size_t size() const noexcept {
    // assert(p.size() == v.size());
    // assert(p.size() == m.size());
    return p.size();
  }

// Allow for iteration over all member variables
  vector_type &get(int i) noexcept {
    vector_type* components[] = {&p.x, &p.y, &p.z, &v.x, &v.y, &v.z, &m};
    return *components[i];
  }

  Particles(size_t size) : p(size), v(size), m(vector_type(size)) {}
  Particles() = default;
};
#endif