#include <math.h>
#include <iostream>
#include <omp.h>

const size_t block = 100ul;

inline double length2(double x, double y, double z) {
  return x * x + y * y + z * z;
}

void bruteForceAcc(double  *accx, double  *accy, double   *accz,
                  const double  *xin, const double  *yin, const double  *zin,
                  const double  x, const double  y, const double  z,
                  const double  *m, const size_t n) {
  /* Tried manual vectorization using simd accumulate here, no measurable
   * difference at -Ofast -march=native. GLM already uses vectorized simd
   * instructions for the accumulation.
   */
  double dvx = 0, dvy = 0, dvz = 0;
  #pragma omp simd
  for (size_t i = 0; i < n; i++)
  {
    double diffx = xin[i] - x;
    double diffy = yin[i] - y;
    double diffz = zin[i] - z;

    constexpr double softening_param = 0.025;
    auto r2 = length2(diffx, diffy, diffz)+softening_param;
    double mOverDist3 = m[i] / (r2 * std::sqrt(r2));
    
    dvx += diffx * mOverDist3;
    dvy += diffy * mOverDist3;
    dvz += diffz * mOverDist3;
  }
  *accx = dvx;
  *accy = dvy;
  *accz = dvz;
}

void stepSimulation(double *x, double *y, double *z, 
                    double *vx, double *vy, double *vz, 
                    double *m, size_t n, double dt)
{
/*Multithreading the outer loop instead of the inner loop is more than 10x
 * faster. This is because the inner loop is very short and the overhead of
 * creating threads is too high.
 */
  
#pragma omp parallel for if (n > 1000)
  for (size_t ii = 0; ii < n-block; ii+=block) {
    // Blocking to improve caching
    for (size_t i = ii; i < ii + block; i++)
    {
      for (size_t jj = 0; jj < n-block; jj+=block)
      {
        double dvx = 0, dvy = 0, dvz = 0;
        bruteForceAcc(&dvx, &dvy, &dvz, x+jj, y+jj, z+jj, x[i], y[i], z[i], m+jj, block);
        vx[i] += dvx*dt;
        vy[i] += dvy*dt; 
        vz[i] += dvz*dt;
      }
    }
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
  size_t n = 200000;
  double dt = 0.1;
  
  double *x = new double[n];
  double *y = new double[n];
  double *z = new double[n];
  double *vx = new double[n];
  double *vy = new double[n];
  double *vz = new double[n];
  double *m = new double[n];

#pragma omp parallel for // First touch initialization
  for (size_t ii = 0; ii < n-block; ii+=block) {
    // Blocking to improve caching
    for (size_t i = ii; i < ii + block; i++)
    {
      x[i] = i*10;
      y[i] = (i%3)*10;
      z[i] = 500-i*10;
      vx[i] = -i;
      vy[i] = 0.2;
      vz[i] = 0.4;
      m[i] = 1;
    }
  }

  for (size_t i = 0; i < 1; i++)
  {
    stepSimulation(x, y, z, vx, vy, vz, m, n, dt);
    std::cout<<"."<<std::flush;
  }
  std::cout << "x: "<<x[0]<<", y: "<<y[0]<<", z: "<<z[0]<<", vx: "<<vx[0]<<", vy: "<<vy[0]<<", vz: "<<vz[0]<<", m: "<<m[0]<<"\n";

}