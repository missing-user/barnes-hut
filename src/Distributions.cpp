#include "Distributions.h"
#include <glm/gtc/random.hpp>
#include <numbers>
#include <random>

thread_local std::mt19937 mt{std::random_device{}()};
std::uniform_real_distribution uniform_dist{-1.0, 1.0};
std::normal_distribution normal_dist{0.0, 1.0};

void set_seed(unsigned int seed) { mt.seed(seed); }

std::vector<Particle> normal_distribution(size_t num_particles)
{
  std::vector<Particle> particles(num_particles);
#pragma omp parallel for
  for (int i = 0; i < num_particles; i++)
  {
    particles[i].p.x = normal_dist(mt);
    particles[i].p.y = normal_dist(mt);
    particles[i].p.z = normal_dist(mt);
  }
  return particles;
}

std::vector<Particle> ball_dist(size_t num_particles)
{
  std::vector<Particle> particles(num_particles);
  #pragma omp parallel for
  for (int i = 0; i < num_particles; i++)
  {
    particles[i].p = glm::ballRand(1.0);
  }
  return particles;
}

std::vector<Particle> sphere_dist(size_t num_particles)
{
  std::vector<Particle> particles(num_particles);
  #pragma omp parallel for
  for (int i = 0; i < num_particles; i++)
  {
    particles[i].p = glm::sphericalRand(1.0);
  }
  return particles;
}

std::vector<Particle> box_distribution(size_t num_particles)
{
  std::vector<Particle> particles(num_particles);
  #pragma omp parallel for
  for (int i = 0; i < num_particles; i++)
  {
    particles[i].p.x = uniform_dist(mt);
    particles[i].p.y = uniform_dist(mt);
    particles[i].p.z = uniform_dist(mt);
  }
  return particles;
}

std::vector<Particle> exponential_disk_distribution(size_t num_particles)
{
  std::vector<Particle> particles(num_particles);

  std::exponential_distribution radial_dist{1.0};
  std::normal_distribution vertical_dist{0.0, 1.0};
  std::uniform_real_distribution angle_dist{0.0, 2 * 3.14159265358};

  // Since we use this distribution for testing and benchmarks, it has to be deterministic. 
  // thread_local random numbers DO NOT GUARANTEE DETERMINISM ACROSS DIFFERENT RUNS! 
  //#pragma omp parallel for
  for (int i = 0; i < num_particles; i++)
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

Particles &set_mass(Particles &particles, myfloat m)
{
  for (auto &mass : particles.m)
  {
    mass = m;
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

  return initial_dist;
}

std::vector<Particle> universe4(int n)
{
  // A disk shaped universe with just 5k particles, will be used for performance
  // eval eventually
  auto initial_dist = exponential_disk_distribution(n);
  const myfloat diameter = 100;

  scale(initial_dist, diameter, diameter, diameter / 10); // flat disk

  add_angular_momentum(initial_dist, myvec3(0, .0, 50));
  return initial_dist;
}

std::vector<Particle> bigbang(int n)
{
  auto initial_dist = ball_dist(n);

  const myfloat diameter = 10;
  scale(initial_dist, diameter, diameter, diameter);
  add_radial_velocity(initial_dist, 100);

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
  // https://en.wikipedia.org/wiki/Plummer_sphere

  std::vector<Particle> particles;
  particles.resize(n);

  const myfloat M = 1;
  const myfloat r_s = 1;

  std::uniform_real_distribution<myfloat> uniform(0, 1);

  #pragma omp parallel for
  for(int i = 0; i < n; i++){
    myfloat r = r_s / std::pow(uniform(mt), 2.0/3.0);
    myfloat theta = 2 * glm::pi<myfloat>() * uniform(mt);
    myfloat phi = glm::pi<myfloat>() * uniform(mt);

    myfloat x = r * std::sin(phi) * std::cos(theta);
    myfloat y = r * std::sin(phi) * std::sin(theta);
    myfloat z = r * std::cos(phi);

    myfloat vx = std::sqrt(M / r) * std::sin(phi) * std::cos(theta);
    myfloat vy = std::sqrt(M / r) * std::sin(phi) * std::sin(theta);
    myfloat vz = std::sqrt(M / r) * std::cos(phi);

    myfloat m = M / static_cast<myfloat>(n);
    
    particles[i] = {{x,y,z}, {vx,vy,vz}, m};
  }

  return particles;
}

Particles make_universe(Distribution dist, size_t num_particles)
{
  std::vector<Particle> particles;
  switch (dist)
  {
  case Distribution::UNIVERSE1:
    particles = universe1();
  break;
  case Distribution::UNIVERSE2:
    particles = universe2();
  break;
  case Distribution::COLLISION:
    particles = collision(num_particles);
  break;
  case Distribution::UNIVERSE4:
    particles = universe4(num_particles);
  break;
  case Distribution::PLUMMER:
    particles = plummer(num_particles);
  break;
  case Distribution::BIGBANG:
    particles = bigbang(num_particles);
  break;
  case Distribution::STABLE_ORBIT:
    particles = stable_orbit();
  break;
  case Distribution::SPHERE:
  {
    particles = sphere_dist(num_particles);
    scale(particles, 100, 100, 100);
  }
  break;
  case Distribution::CRYSTALLINE:
  {
    particles = ball_dist(num_particles);
    scale(particles, 100, 100, 100);
  }
  break;
  case Distribution::DEBUG_CUBE:
  {
    int numx = std::pow(static_cast<double>(num_particles), 0.34);
    int numy = std::pow(static_cast<double>(num_particles), 0.34);
    int numz = std::pow(static_cast<double>(num_particles), 0.34);
    particles.resize(numx*numy*numz);
    for (int x = 0; x < numx; x++)
    {
      for (int y = 0; y < numy; y++)
      {
        for (int z = 0; z < numz; z++)
        {
          particles[x*numz*numy+ y*numz + z] = Particle{{x, y, z}, {0, 0, 0}, 1};
        }
      }
    }
  }
  break;
  default:
    throw std::runtime_error("Unknown distribution");
    return {num_particles};
  }

  Particles bodies{particles.size()};
#pragma omp parallel for // First touch initialization
  for (int i = 0; i < bodies.size(); i++) {
    bodies.p.x[i] = particles[i].p.x;
    bodies.p.y[i] = particles[i].p.y;
    bodies.p.z[i] = particles[i].p.z;
    bodies.v.x[i] = particles[i].v.x;
    bodies.v.y[i] = particles[i].v.y;
    bodies.v.z[i] = particles[i].v.z;
    bodies.m[i] = particles[i].m;
  }
  set_mass(bodies, 50e4/num_particles);

  return bodies;
}
