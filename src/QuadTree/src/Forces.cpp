#include "Forces.h"

myvec3 lennardJones(const myvec3 &diff, myfloat mass)
{
  const myfloat softening_param = 0.2;

  constexpr myfloat epsilon = 1e3; // depth
  constexpr myfloat delta = 0.6;   // distance from 0
  constexpr myfloat A = 4 * epsilon * pow(delta, 11);
  constexpr myfloat B = 4 * epsilon * pow(delta, 5);
  const myfloat r = glm::length2(diff) + softening_param;

  // only ^6 and ^3 terms, length is already squared
  return diff * (-A / pow(r, 6) + B / pow(r, 3));
}

myvec3 gravity(const myvec3 &diff, myfloat mass)
{
  const myfloat softening_param = 0.025;
  return diff * mass /
         (glm::length2(diff) * glm::length(diff) + softening_param);

  // about 5% slower than the above. same using glm::fastNormalize
  // return glm::normalize(diff) * mass / (glm::length2(diff) +
  // softening_param);
}

// This function must be able to evaluate at a distance of zero
myvec3 accelFunc(const myvec3 &diff, myfloat mass)
{
  return gravity(diff, mass) + lennardJones(diff, mass);
}

myfloat energyAtDistance(myfloat distance, myfloat mass)
{
  constexpr myfloat G = 1;
  constexpr myfloat M = 1;
  auto gravityContrib = -G * M / distance;

  return gravityContrib;
}