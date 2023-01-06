
#include "Particle.h"
// class to hold the position and velocity of a particle

Particle operator+(const Particle &P1,
                   const Particle &P2) { // Overload the + operator to add two
                                         // particles position together
  // This is used to calculate the center of mass of a group of particles
  // The velocity of the center of mass is not calculated
  Particle p;
  p.m = P1.m + P2.m;
  if (p.m == 0)
    return p; // If the mass is zero, skip division, return default particle

  p.p = (P1.p * P1.m + P2.p * P2.m) / (p.m);
  return p;
}