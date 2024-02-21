#include "Bruteforce.h"
#include "Forces.h"

void stepSimulation(Particles &particles, myfloat dt) {
/*Multithreading the outer loop instead of the inner loop is more than 10x
 * faster. This is because the inner loop is very short and the overhead of
 * creating threads is too high.
 */
    myfloat half_dt = 0.5 * dt;
    // v_i+1/2 Velocity half-step
#pragma omp parallel for if (particles.size() > 1000)
    for (int i = 0; i < particles.size(); i++) {
      myfloat dvx = 0, dvy = 0, dvz = 0;
      bruteForceAcc<myfloat>(&dvx, &dvy, &dvz, 
                    particles.p.x.data(), particles.p.y.data(), particles.p.z.data(), 
                    particles.p.x[i], particles.p.y[i], particles.p.z[i], 
                    particles.m.data(), particles.size());
      particles.v.x[i] += dvx*half_dt;
      particles.v.y[i] += dvy*half_dt; 
      particles.v.z[i] += dvz*half_dt;
    }

    // x_i+1 Update positions
#pragma omp simd
    for (size_t i = 0; i < particles.size(); i++)
    {
      particles.p.x[i] += particles.v.x[i] * dt;
      particles.p.y[i] += particles.v.y[i] * dt;
      particles.p.z[i] += particles.v.z[i] * dt;
    }

    //  v_i+1 Velocity full-step
#pragma omp parallel for simd if (particles.size() > 1000)
    for (size_t i = 0; i < particles.v.x.size(); ++i) {
      myfloat dvx = 0, dvy = 0, dvz = 0;
      bruteForceAcc<myfloat>(&dvx, &dvy, &dvz, 
                    particles.p.x.data(), particles.p.y.data(), particles.p.z.data(), 
                    particles.p.x[i], particles.p.y[i], particles.p.z[i], 
                    particles.m.data(), particles.size());
        particles.v.x[i] += dvx*half_dt;
        particles.v.y[i] += dvy*half_dt;
        particles.v.z[i] += dvz*half_dt;
    }
}
