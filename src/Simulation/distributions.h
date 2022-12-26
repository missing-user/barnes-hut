#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

#include <numbers>
#include <random>
#include <vector>


#include "Simulation.h"

void set_seed(unsigned int seed);

/*******************************************************************/
std::vector<Particle> normal_distribution(int);

std::vector<Particle> sphere_distribution(int);

std::vector<Particle> box_distribution(int);

std::vector<Particle> exponential_disk_distribution(int);

/*******************************************************************/

std::vector<Particle> &scale(std::vector<Particle> &particles, myfloat x,
                             myfloat y, myfloat z);

std::vector<Particle> &translate(std::vector<Particle> &particles, myfloat x,
                                 myfloat y, myfloat z);

std::vector<Particle> &add_velocity(std::vector<Particle> &particles, myfloat x,
                                    myfloat y, myfloat z);

std::vector<Particle> &add_angular_momentum(std::vector<Particle> &particles,
                                            myvec3);

std::vector<Particle> &add_random_velocity(std::vector<Particle> &, myfloat,
                                           myfloat, myfloat);

std::vector<Particle> &set_mass(std::vector<Particle> &, myfloat m);

std::vector<Particle> &set_maxwell_v_dist(std::vector<Particle> &,
                                          myfloat temperature = 0.1);

std::vector<Particle> &add_radial_velocity(std::vector<Particle> &, myfloat);

/*******************************************************************/

std::vector<Particle> universe2();

std::vector<Particle> universe3();

std::vector<Particle> universe1();

std::vector<Particle> universe4();

std::vector<Particle> bigbang();

std::vector<Particle> stable_orbit();

#endif
