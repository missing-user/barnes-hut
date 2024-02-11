#ifndef FORCES
#define FORCES

#include "Particle.h"
myfloat potentialFunc(const myvec3 &diff, myfloat mass);
void accelFunc(myfloat* accx, myfloat* accy, myfloat* accz, 
              myfloat, myfloat, myfloat, 
              myfloat mass);

#endif