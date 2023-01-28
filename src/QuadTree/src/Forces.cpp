#include "Forces.h"

/* Force functions
******************************************************************/
// Parameters for the lennard jones potential
constexpr myfloat epsilon = 1e2; // depth
constexpr myfloat delta = 0.5;   // optimal distance from 0
constexpr myfloat A = 4 * epsilon * pow(delta, 12);
constexpr myfloat B = 4 * epsilon * pow(delta, 6);

const myfloat lj_softening_param = 0.05;

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

  return diff * (-12 * A / pow(r, 7) + 6 * B / pow(r, 4));
}

myvec3 gravityForce(const myvec3 &diff, myfloat mass) {
  const myfloat softening_param = 0.025;
  return diff * mass /
         (glm::length2(diff) * glm::length(diff) + softening_param);

  // about 5% slower than the above. same using glm::fastNormalize
  // This function would also fail at r=0, due to glm::normalize()
  // return glm::normalize(diff) * mass / (glm::length2(diff) +
  // softening_param);
}

myvec3 accelFunc(const myvec3 &diff, myfloat mass) {
  /*
   * This function is used to compute the acceleration of a particle
   * given the difference vector between the particle attracting particles mass.
   * It must be able to evaluate at a distance of zero for the brute force
   * method to work, since we are not checking for divide by zero (performance
   * reasons).
   */

  return gravityForce(diff, mass); // + lennardJonesForce(diff, mass);
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
  return gravityPotential(diff, mass) + lennardJonesPotential(diff, mass);
}