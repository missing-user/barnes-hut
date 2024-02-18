#include "Bruteforce.h"
#include "Forces.h"

void stepSimulation(Particles &particles, myfloat dt) {
/*Multithreading the outer loop instead of the inner loop is more than 10x
 * faster. This is because the inner loop is very short and the overhead of
 * creating threads is too high.
 */
#pragma omp parallel if (particles.size() > 1000)
  {
    #pragma omp for
    for (size_t i = 0; i < particles.size(); i++) {
      myfloat dvx = 0, dvy = 0, dvz = 0;
      bruteForceAcc<myfloat>(&dvx, &dvy, &dvz, 
                    particles.p.x.data(), particles.p.y.data(), particles.p.z.data(), 
                    particles.p.x[i], particles.p.y[i], particles.p.z[i], 
                    particles.m.data(), particles.size());
      particles.v.x[i] += dvx*dt;
      particles.v.y[i] += dvy*dt; 
      particles.v.z[i] += dvz*dt;
    }

    for (size_t i = 0; i < particles.size(); i++)
    {
      particles.p.x[i] += particles.v.x[i] * dt;
      particles.p.y[i] += particles.v.y[i] * dt;
      particles.p.z[i] += particles.v.z[i] * dt;
    }
  }
}
