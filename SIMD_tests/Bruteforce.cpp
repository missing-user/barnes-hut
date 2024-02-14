#include <math.h>
#include <iostream>
#include <omp.h>

#include "Forces.h"

void stepSimulation(double *x, double *y, double *z, 
                    double *vx, double *vy, double *vz, 
                    double *m, size_t n, double dt)
{
/*Multithreading the outer loop instead of the inner loop is more than 10x
 * faster. This is because the inner loop is very short and the overhead of
 * creating threads is too high.
 */
  
#pragma omp parallel for if (n > 1000)
  for (size_t i = 0; i < n; i++) {
    double dvx = 0, dvy = 0, dvz = 0;
    bruteForceAcc(&dvx, &dvy, &dvz, x, y, z, x[i], y[i], z[i], m, n);
    vx[i] += dvx*dt;
    vy[i] += dvy*dt; 
    vz[i] += dvz*dt;
  }

  for (size_t i = 0; i < n; i++)
  {
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
    z[i] += vz[i] * dt;
  }
}


int main() {
  std::cout << "Running on "<<omp_get_max_threads()<<" threads\n";
  size_t n = 10000;
  double dt = 0.1;
  
  double *x = new double[n];
  double *y = new double[n];
  double *z = new double[n];
  double *vx = new double[n];
  double *vy = new double[n];
  double *vz = new double[n];
  double *m = new double[n];

#pragma omp parallel for // First touch initialization
  for (size_t i = 0; i < n; i++)
  {
    x[i] = i*10;
    y[i] = (i%3)*10;
    z[i] = 500-i*10;
    vx[i] = -i;
    vy[i] = 0.2;
    vz[i] = 0.4;
    m[i] = 1;
  }

  for (size_t i = 0; i < 100; i++)
  {
    stepSimulation(x, y, z, vx, vy, vz, m, n, dt);
    std::cout<<"."<<std::flush;
  }
  std::cout << "x: "<<x[0]<<", y: "<<y[0]<<", z: "<<z[0]<<", vx: "<<vx[0]<<", vy: "<<vy[0]<<", vz: "<<vz[0]<<", m: "<<m[0]<<"\n";

}