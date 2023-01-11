#ifndef FORCES
#define FORCES

#include "Particle.h"
myfloat potentialFunc(const myvec3 &diff, myfloat mass);
myvec3 accelFunc(const myvec3 &diff, myfloat mass);

#endif