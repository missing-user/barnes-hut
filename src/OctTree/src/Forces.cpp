#include "Forces.h"

/* Potential energy
******************************************************************/

myfloat gravityPotential(const myvec3 &diff, myfloat mass) {
  return mass / glm::length(diff);
}

myfloat lennardJonesPotential(const myvec3 &diff, myfloat mass) {
  const myfloat r = glm::length2(diff) + softening_param;

  return A / pow(r, 6) - B / pow(r, 3);
}

myfloat potentialFunc(const myvec3 &diff, myfloat mass) {

  #ifndef LENNARD_JONES
    return gravityPotential(diff, mass);
  #else
    return gravityPotential(diff, mass) + lennardJonesPotential(diff, mass);
  #endif  
}