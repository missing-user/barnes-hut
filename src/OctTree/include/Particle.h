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

typedef glm::dvec3 myvec3; // We can easily switch the entire implementation to
                           // float precision by adjusting these two variables
typedef double myfloat;

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
  myfloat *x, *y, *z;

  Vectors(size_t size) {
    x = (myfloat*)aligned_alloc(64, sizeof(myfloat)*size);
    y = (myfloat*)aligned_alloc(64, sizeof(myfloat)*size);
    z = (myfloat*)aligned_alloc(64, sizeof(myfloat)*size);
  }

  ~Vectors() {
    free(x);
    free(y);
    free(z);
  }
};

class Particles{
public:
  Vectors p, v;
  myfloat *m;
  Vectors p2; // For the next timestep, should this really be here?
private:
  size_t count;
public:
  size_t size() const {
    return count;
  }

  Particles(size_t size) : p(size), v(size), p2(size), count(size) {
    m = (myfloat*)aligned_alloc(64, sizeof(myfloat)*count);
  }

  ~Particles() {
    free(m);
  }
};
#endif