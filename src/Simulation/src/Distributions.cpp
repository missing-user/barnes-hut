#include "Distributions.h"
#include <numbers>
#include <random>

std::mt19937 mt{std::random_device{}()};
std::uniform_real_distribution uniform_dist{-1.0, 1.0};
std::normal_distribution normal_dist{0.0, 1.0};

void set_seed(unsigned int seed) { mt.seed(seed); }

std::vector<Particle> normal_distribution(int num_particles) {
  std::vector<Particle> particles(num_particles);
  for (size_t i = 0; i < num_particles; i++) {
    particles[i].p.x = normal_dist(mt);
    particles[i].p.y = normal_dist(mt);
    particles[i].p.z = normal_dist(mt);
  }
  return particles;
}

std::vector<Particle> sphere_distribution(int num_particles) {
  // https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability/87238#87238
  std::vector<Particle> particles(num_particles);
  for (size_t i = 0; i < num_particles; i++) {
    myfloat w;
    w = cbrt(uniform_dist(mt));
    myvec3 p{normal_dist(mt), normal_dist(mt), normal_dist(mt)};
    p = glm::normalize(p) * w;

    particles[i].p = p;
  }
  return particles;
}

std::vector<Particle> box_distribution(int num_particles) {
  std::vector<Particle> particles(num_particles);
  for (size_t i = 0; i < num_particles; i++) {
    particles[i].p.x = uniform_dist(mt);
    particles[i].p.y = uniform_dist(mt);
    particles[i].p.z = uniform_dist(mt);
  }
  return particles;
}

std::vector<Particle> exponential_disk_distribution(int num_particles) {
  std::vector<Particle> particles(num_particles);

  std::exponential_distribution radial_dist{1.0};
  std::normal_distribution vertical_dist{0.0, 1.0};
  std::uniform_real_distribution angle_dist{0.0, 2 * std::numbers::pi};

  for (size_t i = 0; i < num_particles; i++) {
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
                             myfloat y, myfloat z) {
  myvec3 scale{x, y, z};
  for (auto &p : particles) {
    p.p *= scale;
  }
  return particles;
}

std::vector<Particle> &translate(std::vector<Particle> &particles, myfloat x,
                                 myfloat y, myfloat z) {
  myvec3 scale{x, y, z};
  for (auto &p : particles) {
    p.p += scale;
  }
  return particles;
}

std::vector<Particle> &add_velocity(std::vector<Particle> &particles, myfloat x,
                                    myfloat y, myfloat z) {
  myvec3 scale{x, y, z};
  for (auto &p : particles) {
    p.v += scale;
  }
  return particles;
}

std::vector<Particle> &add_angular_momentum(std::vector<Particle> &particles,
                                            myvec3 axis) {
  for (auto &p : particles) {
    p.v += glm::cross(p.p, axis);
  }
  return particles;
}

std::vector<Particle> &add_random_velocity(std::vector<Particle> &particles,
                                           myvec3 scale) {

  for (auto &p : particles) {
    myvec3 velo{normal_dist(mt), normal_dist(mt), normal_dist(mt)};
    p.v += scale * velo;
  }
  return particles;
}

std::vector<Particle> &set_mass(std::vector<Particle> &particles, myfloat m) {
  for (auto &p : particles) {
    p.m = m;
  }
  return particles;
}

std::vector<Particle> &set_maxwell_v_dist(std::vector<Particle> &particles,
                                          myfloat temperature) {
  // https://milianw.de/code-snippets/maxwell-distribution-in-c11.html
  // Boltzmann factor times temperature
  const myfloat k_T = 0.1;
  // setup the Maxwell distribution, i.e. gamma distribution with alpha = 3/2
  std::gamma_distribution<myfloat> maxwell(3. / 2., k_T);

  for (auto &p : particles) {
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
                                           myfloat scale) {
  for (auto &p : particles) {
    p.v = p.p * scale;
  }
  return particles;
}

/*******************************************************************/

std::vector<Particle> universe2() {
  auto initial_dist = sphere_distribution(100);
  scale(initial_dist, 100, 100, 100);
  set_mass(initial_dist, 10);

  add_angular_momentum(initial_dist, myvec3(0, .1, 0));
  initial_dist.push_back(Particle{{0, 0, 0}, {0, 0, 0}, 1e4});
  return initial_dist;
}

std::vector<Particle> universe1() {
  // A disk shaped universe with just 32 particles, useful for simple testing
  auto initial_dist = exponential_disk_distribution(32);
  const myfloat diameter = 100;

  scale(initial_dist, diameter, diameter, 10); // flat disk
  set_mass(initial_dist, 200);

  // add_angular_momentum(initial_dist, myvec3(0, .0, .1) / diameter);
  return initial_dist;
}

std::vector<Particle> universe4(int n) {
  // A disk shaped universe with just 5k particles, will be used for performance
  // eval eventually
  auto initial_dist = exponential_disk_distribution(n);
  const myfloat diameter = 100;

  scale(initial_dist, diameter, diameter, diameter / 10); // flat disk
  set_mass(initial_dist, 50);

  add_angular_momentum(initial_dist, myvec3(0, .0, 10) / diameter);
  return initial_dist;
}

std::vector<Particle> bigbang(int n) {
  auto initial_dist = sphere_distribution(n);

  const myfloat diameter = 10;
  scale(initial_dist, diameter, diameter, diameter);
  add_radial_velocity(initial_dist, 100);
  // add_random_velocity(initial_dist, {0.1, 0.1, 0.1});

  set_mass(initial_dist, 30);

  return initial_dist;
}

std::vector<Particle> stable_orbit() {
  Particle p1 = {{0, 0, 100}, {0, 10, 0}, 0};
  Particle p2 = {{0, 0, 0}, {0, 0, 0}, 10000};
  Particle p3 = {{0, 0, -100}, {0, -10, 0}, 0};

  std::vector<Particle> particles{p1, p2, p3};

  return particles;
}

std::vector<Particle> collision(int n) {
  auto first_half = universe4(n / 2);
  auto second_half = universe4(n / 2);

  translate(second_half, 0, 0, 200);
  first_half.insert(first_half.end(), second_half.begin(), second_half.end());
  return first_half;
}

std::vector<Particle> make_universe(Distribution dist, int num_particles) {
  switch (dist) {
  case Distribution::UNIVERSE1:
    return universe1();
  case Distribution::UNIVERSE2:
    return universe2();
  case Distribution::COLLISION:
    return collision(num_particles);
  case Distribution::UNIVERSE4:
    return universe4(num_particles);
  case Distribution::BIGBANG:
    return bigbang(num_particles);
  case Distribution::STABLE_ORBIT:
    return stable_orbit();
  default:
    return {};
  }
}