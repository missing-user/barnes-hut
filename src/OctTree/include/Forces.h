#ifndef FORCES
#define FORCES

#include "Particle.h"
myfloat potentialFunc(const myvec3 &diff, myfloat mass);
void accelFunc(myfloat *__restrict  accx, myfloat *__restrict  accy, myfloat *__restrict  accz, 
              myfloat, myfloat, myfloat, 
              myfloat mass);


void bruteForceAcc(double *__restrict accx, double  *__restrict accy, double   *__restrict accz,
                  const double  *__restrict xin, const double  *__restrict yin, const double  *__restrict zin,
                  const double  x, const double  y, const double  z,
                  const double  *__restrict m, const size_t n);
#endif