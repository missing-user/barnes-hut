#include "Bruteforce.h"

//#pragma omp declare simd
void bruteForceAcc(myfloat __restrict__ *accx, myfloat __restrict__ *accy, myfloat  __restrict__ *accz,
                    const myfloat __restrict__ *xin, const myfloat  __restrict__ *yin, const myfloat __restrict__ *zin,
                    const myfloat __restrict__ *m, const size_t count, const size_t me) {
  /* Tried manual vectorization using simd accumulate here, no measurable
   * difference at -Ofast -march=native. GLM already uses vectorized simd
   * instructions for the accumulation.
   */

  const myfloat softening_param = 0.025;
  // #pragma omp simd reduction(+:accx[me], accy[me], accz[me])
  // simd reduction made it slower
  for (size_t i = 0; i < count; i++)
  {
    myfloat diffx = xin[i] - xin[me];
    myfloat diffy = yin[i] - yin[me];
    myfloat diffz = zin[i] - zin[me];
    myfloat r2 = diffx*diffx + diffy*diffy + diffz*diffz;
    myfloat r = std::sqrt(r2);
    accx[me] += diffx * m[i] / (r2*r + softening_param);
    accy[me] += diffy * m[i] / (r2*r + softening_param);
    accz[me] += diffz * m[i] / (r2*r + softening_param);
  }
  
  accx[me] -= 0* m[me] / (softening_param);
  accy[me] -= 0* m[me] / (softening_param);
  accz[me] -= 0* m[me] / (softening_param);
}

void stepSimulation(Particles &particles,
                                     myfloat dt) {

  /* Computes the state vector of all particles for the next timestep,
   * integrating position and velocity by summing up all the individual particle
   * contributions (n^2 brute force algorithm)
   */

/*Multithreading the outer loop instead of the inner loop is more than 10x
 * faster. This is because the inner loop is very short and the overhead of
 * creating threads is too high.
 */
  myfloat *accx, *accy, *accz;
  accx = (myfloat*)aligned_alloc(64, sizeof(myfloat)*particles.count);
  accy = (myfloat*)aligned_alloc(64, sizeof(myfloat)*particles.count);
  accz = (myfloat*)aligned_alloc(64, sizeof(myfloat)*particles.count);
  memset(accx, 0, particles.count * sizeof(myfloat));
  memset(accy, 0, particles.count * sizeof(myfloat));
  memset(accz, 0, particles.count * sizeof(myfloat));
#pragma omp parallel if (particles.count > 1000)
  {
#pragma omp for
    for (size_t i = 0; i < particles.count; i++) {

      bruteForceAcc(accx, accy, accz,
                    particles.x, particles.y, particles.z, 
                    particles.m, particles.count, i);
      particles.x2[i] = particles.x[i] + particles.vx[i] * dt + accx[i] * dt * dt / 2.;
      particles.y2[i] = particles.y[i] + particles.vy[i] * dt + accy[i] * dt * dt / 2.;
      particles.z2[i] = particles.z[i] + particles.vz[i] * dt + accz[i] * dt * dt / 2.;
    }

#pragma omp for
    for (size_t i = 0; i < particles.count; i++) {
      // Then update the velocities using v(t+1) = dt*(a(t) + a(t+dt))/2
      bruteForceAcc(accx, accy, accz,
                    particles.x2, particles.y2, particles.z2, 
                    particles.m, particles.count, i);
      particles.vx[i] += accx[i] * dt / 2.;
      particles.vy[i] += accy[i] * dt / 2.;
      particles.vz[i] += accz[i] * dt / 2.;
    }
  }

  delete[] accx;
  delete[] accy;
  delete[] accz;
}
