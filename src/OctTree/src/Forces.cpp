#include "Forces.h"

/* Force functions
******************************************************************/
// Parameters for the lennard jones potential
constexpr myfloat epsilon = 1e2; // depth
constexpr myfloat delta = 0.5;   // optimal distance from 0
const myfloat A = 4 * epsilon * std::pow(delta, 12);
const myfloat B = 4 * epsilon * std::pow(delta, 6);

const myfloat lj_softening_param = 0.05;

// #define LENNARD_JONES

myvec3 lennardJonesForce(const myvec3 &diff, myfloat mass) {
  const myfloat r = glm::length2(diff) + lj_softening_param;

  /* !!! diff is not normalized (|diff| is r) !!!
   *
   * The force is the derivative of the potential energy:
   * F = -dU/dr * rhat = -dU/dr * diff/r
   * F = rhat * (12*A*r^(-13) - 6*B*r^(-7))     (rhat is the unit vector)
   * F = diff * (12*A*r^(-14) - 6*B*r^(-8))     (diff is the vector p1-p2)
   * And since we are using r2=r^2 to save a sqrt() operation:
   * F = diff * (12*A*r2^(-7) - 6*B*r2^(-4))
   */

  // For some reason the compiler does not optimize pow(r,7) and pow(r,4) in this case
  // So writing the expression out explicitly is about 70% faster
  return diff * (-12 * A / (r*r*r*r*r*r*r) + 6 * B / (r*r*r*r));
}

#pragma omp declare simd
void gravityForce(myfloat* accx, myfloat* accy, myfloat* accz, 
  myfloat diffx, myfloat diffy, myfloat diffz, myfloat mass) {
  constexpr myfloat softening_param = 0.025;
  auto r2 = length2(diffx, diffy, diffz) + softening_param;
  auto mOverDist3 = mass / (r2 * std::sqrt(r2));
  
  *accx += diffx * mOverDist3;
  *accy += diffy * mOverDist3;
  *accz += diffz * mOverDist3;
  // about 5% slower than the above. same using glm::fastNormalize
  // This function would also fail at r=0, due to glm::normalize()
  // return glm::normalize(diff) * mass / (glm::length2(diff) +
  // softening_param);
}

#pragma omp declare simd
void accelFunc(myfloat* accx, myfloat* accy, myfloat* accz, 
  myfloat diffx, myfloat diffy, myfloat diffz, myfloat mass) {
  /*
   * This function is used to compute the acceleration of a particle
   * given the difference vector between the particle attracting particles mass.
   * It must be able to evaluate at a myfloat of zero for the brute force
   * method to work, since we are not checking for divide by zero (performance
   * reasons).
   */

  #ifndef LENNARD_JONES
    gravityForce(accx, accy, accz, diffx, diffy, diffz, mass);
  #else
    return gravityForce(diff, mass) + lennardJonesForce(diff, mass);
  #endif

}

void bruteForceAcc(double  *accx, double  *accy, double   *accz,
                  const double  *xin, const double  *yin, const double  *zin,
                  const double  x, const double  y, const double  z,
                  const double  *m, const size_t n) {
  /* Tried manual vectorization using simd accumulate here, no measurable
   * difference at -Ofast -march=native. GLM already uses vectorized simd
   * instructions for the accumulation.
   */
  double dvx = 0, dvy = 0, dvz = 0;
  #pragma omp simd
  for (size_t i = 0; i < n; i++)
  {
    accelFunc(&dvx, &dvy, &dvz, xin[i] - x, yin[i] - y, zin[i] - z, m[i]);
  }
  *accx += dvx;
  *accy += dvy;
  *accz += dvz;
}

/* Potential energy
******************************************************************/

myfloat gravityPotential(const myvec3 &diff, myfloat mass) {
  return mass / glm::length(diff);
}

myfloat lennardJonesPotential(const myvec3 &diff, myfloat mass) {
  const myfloat r = glm::length2(diff) + lj_softening_param;

  return A / pow(r, 6) - B / pow(r, 3);
}

myfloat potentialFunc(const myvec3 &diff, myfloat mass) {

  #ifndef LENNARD_JONES
    return gravityPotential(diff, mass);
  #else
    return gravityPotential(diff, mass) + lennardJonesPotential(diff, mass);
  #endif  
}