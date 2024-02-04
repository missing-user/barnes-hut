#ifndef PARTICLE
#define PARTICLE

#include <sstream>
#include <vector>

#define GLM_FORCE_INTRINSICS 
//#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES // Slowdowns on my machine
//#define GLM_FORCE_INLINE // Didn't make a measurable difference
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp> // Allows us to easily std::cout << pvec3;

typedef glm::dvec3 myvec3; // We can easily switch the entire implementation to
                           // float precision by adjusting these two variables
typedef double myfloat;

/*
template <typename T>
concept HasPosition = requires(T t) {
  { t.p } -> std::convertible_to<myvec3>;
};

template <typename T>
concept HasMassAndPosition = requires(T t) {
  { t.m } -> std::convertible_to<myfloat>;
}
&&HasPosition<T>;
*/

struct CenterOfMass {
  myvec3 p;
  myfloat m;

  template <typename T>
  // Adding two centers of mass creates the center of mass between the two
  CenterOfMass &operator+=(const T &rhs) {
    const auto total_m = this->m + rhs.m;
    if (total_m == 0)
      return *this; // If the mass is zero, skip division, return original

    this->p = (this->p * this->m + rhs.p * rhs.m) / total_m;
    this->m = total_m;
    return *this;
  }
};

struct Particle { // A particle with position, velocity and unique id
  myvec3 p;
  myvec3 v;
  myfloat m;
};

struct Particles{
  myfloat *x, *y, *z;
  myfloat *vx, *vy, *vz;
  myfloat *x2, *y2, *z2; // Temp storage for the second step of the integration
  myfloat *m;
  size_t count;

  Particle operator[](size_t i) const {
    return {{x[i], y[i], z[i]}, {vx[i], vy[i], vz[i]}, m[i]};
  }
};


// Helper function for the particle class, so we can print and debug it easily
inline std::ostream &operator<<(std::ostream &out, const CenterOfMass &p) {
  return out << p.p << "\tm=" << p.m;
}
inline std::ostream &operator<<(std::ostream &out, const Particle &p) {
  return out << p.p << "\tv=" << p.v << "\tm=" << p.m;
}

// Print the resulting values of all particles in a .csv format
inline std::ostream &operator<<(std::ostream &out,
                                const std::vector<Particle> &particles) {
  for (auto &p : particles) {
    out << p.p[0] << "," << p.p[1] << "," << p.p[2] << ","<< p.m<< "\n";
  }

  return out;
}

#endif