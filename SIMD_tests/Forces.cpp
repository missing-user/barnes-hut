#include "Forces.h"
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