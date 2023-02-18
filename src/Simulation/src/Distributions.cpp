#include "Distributions.h"
#include <glm/gtc/random.hpp>
#include <numbers>
#include <random>

thread_local std::mt19937 mt{std::random_device{}()};
std::uniform_real_distribution uniform_dist{-1.0, 1.0};
std::normal_distribution normal_dist{0.0, 1.0};

void set_seed(unsigned int seed) { mt.seed(seed); }

std::vector<Particle> normal_distribution(int num_particles)
{
  std::vector<Particle> particles(num_particles);
  #pragma omp parallel for
  for (size_t i = 0; i < num_particles; i++)
  {
    particles[i].p.x = normal_dist(mt);
    particles[i].p.y = normal_dist(mt);
    particles[i].p.z = normal_dist(mt);
  }
  return particles;
}

std::vector<Particle> ball_dist(int num_particles)
{
  std::vector<Particle> particles(num_particles);
  #pragma omp parallel for
  for (size_t i = 0; i < num_particles; i++)
  {
    particles[i].p = glm::ballRand(1.0);
  }
  return particles;
}

std::vector<Particle> sphere_dist(int num_particles)
{
  std::vector<Particle> particles(num_particles);
  #pragma omp parallel for
  for (size_t i = 0; i < num_particles; i++)
  {
    particles[i].p = glm::sphericalRand(1.0);
  }
  return particles;
}

std::vector<Particle> box_distribution(int num_particles)
{
  std::vector<Particle> particles(num_particles);
  #pragma omp parallel for
  for (size_t i = 0; i < num_particles; i++)
  {
    particles[i].p.x = uniform_dist(mt);
    particles[i].p.y = uniform_dist(mt);
    particles[i].p.z = uniform_dist(mt);
  }
  return particles;
}

std::vector<Particle> exponential_disk_distribution(int num_particles)
{
  std::vector<Particle> particles(num_particles);

  std::exponential_distribution radial_dist{1.0};
  std::normal_distribution vertical_dist{0.0, 1.0};
  std::uniform_real_distribution angle_dist{0.0, 2 * 3.14159265358};

  // Since we use this distribution for testing and benchmarks, it has to be deterministic. 
  // thread_local random numbers DO NOT GUARANTEE DETERMINISM ACROSS DIFFERENT RUNS! 
  //#pragma omp parallel for
  for (size_t i = 0; i < num_particles; i++)
  {
    auto r = radial_dist(mt);
    auto phi = angle_dist(mt);
    auto z = vertical_dist(mt);

    particles[i].p.x = r * cos(phi);
    particles[i].p.y = r * sin(phi);
    particles[i].p.z = z;
  }
  return particles;
}

/*******************************************************************/

std::vector<Particle> &scale(std::vector<Particle> &particles, myfloat x,
                             myfloat y, myfloat z)
{
  myvec3 scale{x, y, z};
  for (auto &p : particles)
  {
    p.p *= scale;
  }
  return particles;
}

std::vector<Particle> &translate(std::vector<Particle> &particles, myfloat x,
                                 myfloat y, myfloat z)
{
  myvec3 scale{x, y, z};
  for (auto &p : particles)
  {
    p.p += scale;
  }
  return particles;
}

std::vector<Particle> &add_velocity(std::vector<Particle> &particles, myfloat x,
                                    myfloat y, myfloat z)
{
  myvec3 scale{x, y, z};
  for (auto &p : particles)
  {
    p.v += scale;
  }
  return particles;
}

std::vector<Particle> &add_angular_momentum(std::vector<Particle> &particles,
                                            myvec3 axis)
{
  const auto mag = glm::length(axis);
  for (auto &p : particles)
  {
    p.v += glm::normalize(glm::cross(p.p, axis)) * mag;
  }
  return particles;
}

std::vector<Particle> &set_mass(std::vector<Particle> &particles, myfloat m)
{
  for (auto &p : particles)
  {
    p.m = m;
  }
  return particles;
}

std::vector<Particle> &set_maxwell_v_dist(std::vector<Particle> &particles,
                                          myfloat temperature)
{
  // https://milianw.de/code-snippets/maxwell-distribution-in-c11.html
  // Boltzmann factor times temperature
  const myfloat k_T = 0.1;
  // setup the Maxwell distribution, i.e. gamma distribution with alpha = 3/2
  std::gamma_distribution<myfloat> maxwell(3. / 2., k_T);

  #pragma omp parallel for
  for (auto &p : particles)
  {
    myfloat x, y, z;
    x = uniform_dist(mt);
    y = uniform_dist(mt);
    z = uniform_dist(mt);
    myvec3 velo{x, y, z};

    velo = glm::normalize(velo) * maxwell(mt);

    p.v = velo;
  }
  return particles;
}

std::vector<Particle> &add_radial_velocity(std::vector<Particle> &particles,
                                           myfloat scale)
{
  for (auto &p : particles)
  {
    p.v = p.p * scale;
  }
  return particles;
}

/*******************************************************************/

std::vector<Particle> universe2()
{
  auto initial_dist = ball_dist(100);
  scale(initial_dist, 100, 100, 100);
  set_mass(initial_dist, 10);

  add_angular_momentum(initial_dist, myvec3(0, .1, 0));
  initial_dist.push_back(Particle{{0, 0, 0}, {0, 0, 0}, 1e4});
  return initial_dist;
}

std::vector<Particle> universe1()
{
  // A disk shaped universe with just 32 particles, useful for simple testing
  auto initial_dist = exponential_disk_distribution(32);
  const myfloat diameter = 100;

  scale(initial_dist, diameter, diameter, 10); // flat disk
  set_mass(initial_dist, 200);

  return initial_dist;
}

std::vector<Particle> universe4(int n)
{
  // A disk shaped universe with just 5k particles, will be used for performance
  // eval eventually
  auto initial_dist = exponential_disk_distribution(n);
  const myfloat diameter = 100;

  scale(initial_dist, diameter, diameter, diameter / 10); // flat disk
  set_mass(initial_dist, 50);

  add_angular_momentum(initial_dist, myvec3(0, .0, 50));
  return initial_dist;
}

std::vector<Particle> bigbang(int n)
{
  auto initial_dist = ball_dist(n);

  const myfloat diameter = 10;
  scale(initial_dist, diameter, diameter, diameter);
  add_radial_velocity(initial_dist, 100);

  set_mass(initial_dist, 30);

  return initial_dist;
}

std::vector<Particle> stable_orbit()
{
  Particle p1 = {{0, 0, 100}, {0, 10, 0}, 0};
  Particle p2 = {{0, 0, 0}, {0, 0, 0}, 10000};
  Particle p3 = {{0, 0, -100}, {0, -10, 0}, 0};

  std::vector<Particle> particles{p1, p2, p3};

  return particles;
}

std::vector<Particle> collision(int n)
{
  auto first_half = universe4(n / 2);
  auto second_half = universe4(n / 2);

  translate(second_half, 0, 1000, 200);
  add_velocity(second_half, 0, -10, -2);
  add_velocity(first_half, 0, 10, 2);
  first_half.insert(first_half.end(), second_half.begin(), second_half.end());
  return first_half;
}

std::vector<Particle> plummer(int n){
  // https://en.wikipedia.org/wiki/Plummer_model

  std::vector<Particle> particles{n};

  const myfloat M = 1;
  const myfloat r_s = 1;
  const myfloat m = M / static_cast<myfloat>(n);
  
  std::uniform_real_distribution<myfloat> uniform(0, 1);

  for(int i = 0; i < n; i++){
    const myfloat r = r_s / std::sqrt(std::pow(uniform(mt), -2.0/3.0) - 1);
    const myfloat theta = std::acos(uniform_dist(mt));
    const myfloat phi = 2 * glm::pi<myfloat>() * uniform(mt);

    const myfloat px = r * std::sin(theta) * std::cos(phi);
    const myfloat py = r * std::sin(theta) * std::sin(phi);
    const myfloat pz = r * std::cos(theta);

    myfloat x = 0.0;
    myfloat y = 0.1;
    while (y > x*x*std::pow(1.0-x*x,3.5)){
      x = uniform(mt);
      y = uniform(mt)*0.1;
    }
    const myfloat velocity = x * std::sqrt(2.0) * std::pow(1.0 + r*r, -0.25);

    const myfloat vx = velocity /r * px;
    const myfloat vy = velocity /r * py;
    const myfloat vz = velocity /r * pz;

    particles[i] = Particle({px,py,pz}, {vx,vy,vz}, m, i);
  }

  return particles;
}

std::vector<Particle> make_universe(Distribution dist, int num_particles)
{
  switch (dist)
  {
  case Distribution::UNIVERSE1:
    return universe1();
  case Distribution::UNIVERSE2:
    return universe2();
  case Distribution::COLLISION:
    return collision(num_particles);
  case Distribution::UNIVERSE4:
    return universe4(num_particles);
  case Distribution::PLUMMER:
    return plummer(num_particles);
  case Distribution::BIGBANG:
    return bigbang(num_particles);
  case Distribution::STABLE_ORBIT:
    return stable_orbit();
  case Distribution::SPHERE:
  {
    auto particles = sphere_dist(num_particles);
    scale(particles, 100, 100, 100);
    set_mass(particles, 10);
    return particles;
  }
  case Distribution::CRYSTALLINE:
  {
    auto particles = ball_dist(num_particles);
    scale(particles, 100, 100, 100);
    set_mass(particles, 10);
    return particles;
  }
  default:
    return {};
  }
}
