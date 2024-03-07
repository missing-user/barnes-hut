#ifndef FORCES
#define FORCES

#include "Particle.h"
#include "xsimd/xsimd.hpp"

namespace xs = xsimd;
myfloat potentialFunc(const myvec3 &diff, myfloat mass);

// #define LENNARD_JONES

constexpr myfloat epsilon = 1e2; // depth
constexpr myfloat delta = 0.5;   // optimal distance from 0
const myfloat A = 4 * epsilon * std::pow(delta, 12);
const myfloat B = 4 * epsilon * std::pow(delta, 6);

const myfloat softening_param = 0.025;

template <typename T>
void accelFunc(T *__restrict  accx, T *__restrict  accy, T *__restrict  accz, 
              T dx, T dy, T dz, T mass) {
  /*
  * This function is used to compute the acceleration of a particle
  * given the difference vector between the particle attracting particles mass.
  * It must be able to evaluate at a myfloat of zero for the brute force
  * method to work, since we are not checking for divide by zero (performance
  * reasons).
  */

  #ifndef LENNARD_JONES
  auto r2 = length2(dx, dy, dz) + softening_param;
  T r;
  if constexpr(xs::is_batch<T>::value) {
    r = xs::sqrt(r2);
  }else{
    r = std::sqrt(r2);
  }
  auto mOverDist3 = mass / (r2 * r);
  
  #else
  mOverDist3 += mass*(-12 * A / (r2*r2*r2*r) + 6 * B / (r2*r2))
  #endif

  *accx += dx * mOverDist3;
  *accy += dy * mOverDist3;
  *accz += dz * mOverDist3;
}


template <typename T>
void bruteForceAcc(T *__restrict accx, T  *__restrict accy, T   *__restrict accz,
                  const myfloat  *__restrict xin, const myfloat  *__restrict yin, const myfloat  *__restrict zin,
                  const T  x, const T  y, const T  z,
                  const myfloat  *__restrict m, const size_t n){
  T dvx = 0, dvy = 0, dvz = 0;
  #pragma omp simd
  for (size_t i = 0; i < n; i++)
  {
    accelFunc(&dvx, &dvy, &dvz, xin[i] - x, yin[i] - y, zin[i] - z, static_cast<T>(m[i]));
  }
  *accx += dvx;
  *accy += dvy;
  *accz += dvz;
}

#endif