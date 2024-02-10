#include "Bruteforce.h"
#include "Forces.h"
//#pragma omp declare simd

void bruteForceAcc(myfloat  *accx, myfloat  *accy, myfloat   *accz,
                    const myfloat  *xin, const myfloat   *yin, const myfloat  *zin,
                    const myfloat  *m, const size_t n, const size_t me) {
  /* Tried manual vectorization using simd accumulate here, no measurable
   * difference at -Ofast -march=native. GLM already uses vectorized simd
   * instructions for the accumulation.
   */
  //#pragma omp simd reduction(+:accx[me], accy[me], accz[me])
  for (size_t i = 0; i < n; i++)
  {
    myfloat diffx = xin[i] - xin[me];
    myfloat diffy = yin[i] - yin[me];
    myfloat diffz = zin[i] - zin[me];
    accelFunc(&accx[me], &accy[me], &accz[me], diffx, diffy, diffz, m[i]);
  }
}

void stepSimulation(Particles &particles, myfloat dt) {

  /* Computes the state vector of all particles for the next timestep,
   * integrating position and velocity by summing up all the individual particle
   * contributions (n^2 brute force algorithm)
   */

/*Multithreading the outer loop instead of the inner loop is more than 10x
 * faster. This is because the inner loop is very short and the overhead of
 * creating threads is too high.
 */
  Vectors acc{particles.size(), 0};
  Vectors p2{particles.size()};
#pragma omp parallel if (particles.size() > 1000)
  {
#pragma omp for 
    for (size_t i = 0; i < particles.size(); i++) {

      bruteForceAcc(acc.x, acc.y, acc.z,
                    particles.p.x, particles.p.y, particles.p.z, 
                    particles.m, particles.size(), i);
      p2.x[i] = particles.p.x[i] + particles.v.x[i] * dt + acc.x[i] * dt * dt / 2.;
      p2.y[i] = particles.p.y[i] + particles.v.y[i] * dt + acc.y[i] * dt * dt / 2.;
      p2.z[i] = particles.p.z[i] + particles.v.z[i] * dt + acc.z[i] * dt * dt / 2.;
    }

#pragma omp for
    for (size_t i = 0; i < particles.size(); i++) {
      // Then update the velocities using v(t+1) = dt*(a(t) + a(t+dt))/2
      bruteForceAcc(acc.x, acc.y, acc.z,
                    p2.x,p2.y,p2.z,
                    particles.m, particles.size(), i);
      particles.v.x[i] += acc.x[i] * dt / 2.;
      particles.v.y[i] += acc.y[i] * dt / 2.;
      particles.v.z[i] += acc.z[i] * dt / 2.;
    }

#pragma omp for
    for (size_t i = 0; i < particles.size(); i++) {
      particles.p.x[i] = p2.x[i];
      particles.p.y[i] = p2.y[i];
      particles.p.z[i] = p2.z[i];
    }
  }
}
