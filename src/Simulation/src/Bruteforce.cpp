#include "Bruteforce.h"
#include "Forces.h"

void bruteForceAcc(myfloat  *accx, myfloat  *accy, myfloat   *accz,
                  const myfloat  *xin, const myfloat  *yin, const myfloat  *zin,
                  const myfloat  x, const myfloat  y, const myfloat  z,
                  const myfloat  *m, const size_t n) {
  myfloat dvx = 0, dvy = 0, dvz = 0;
  #pragma omp simd
  for (size_t i = 0; i < n; i++)
  {
    auto diffx = xin[i] - x;
    auto diffy = yin[i] - y;
    auto diffz = zin[i] - z;

    constexpr myfloat softening_param = 0.025;
    auto r2 = length2(diffx, diffy, diffz)+softening_param;
    myfloat mOverDist3 = m[i] / (r2 * std::sqrt(r2));
    
    dvx += diffx * mOverDist3;
    dvy += diffy * mOverDist3;
    dvz += diffz * mOverDist3;
  }
  *accx = dvx;
  *accy = dvy;
  *accz = dvz;
}

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
      bruteForceAcc(&dvx, &dvy, &dvz, 
                    particles.p.x, particles.p.y, particles.p.z, 
                    particles.p.x[i], particles.p.y[i], particles.p.z[i], 
                    particles.m, particles.size());
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
