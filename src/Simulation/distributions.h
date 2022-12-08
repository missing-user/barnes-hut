#pragma once

#include <random>
#include <vector>

#include "simulation.h"

const double PI = 3.141592653589793238462643383279502884197169;

std::mt19937 mt{std::random_device{}()};
std::uniform_real_distribution uniform_dist{-1.0, 1.0};
std::normal_distribution normal_dist{0.0, 1.0};

static std::vector<Particle> normal_distribution(int num_particles)
{
	std::vector<Particle> particles(num_particles);
	for (size_t i = 0; i < num_particles; i++)
	{
		particles[i].p.x = normal_dist(mt);
		particles[i].p.y = normal_dist(mt);
		particles[i].p.z = normal_dist(mt);
	}
	return particles;
}

static std::vector<Particle> sphere_distribution(int num_particles)
{
	// https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability/87238#87238
	std::vector<Particle> particles(num_particles);
	for (size_t i = 0; i < num_particles; i++)
	{
		myfloat w, x, y, z;
		w = cbrt(uniform_dist(mt));
		x = normal_dist(mt);
		y = normal_dist(mt);
		z = normal_dist(mt);
		auto mag = sqrt(x * x + y * y + z * z);

		x = w * x / mag;
		y = w * y / mag;
		z = w * z / mag;

		particles[i].p.x = x;
		particles[i].p.y = y;
		particles[i].p.z = z;
	}
	return particles;
}

static std::vector<Particle> box_distribution(int num_particles)
{
	std::vector<Particle> particles(num_particles);
	for (size_t i = 0; i < num_particles; i++)
	{
		particles[i].p.x = uniform_dist(mt);
		particles[i].p.y = uniform_dist(mt);
		particles[i].p.z = uniform_dist(mt);
	}
	return particles;
}

static std::vector<Particle> exponential_disk_distribution(int num_particles)
{
	std::vector<Particle> particles(num_particles);

	std::exponential_distribution radial_dist{1.0};
	std::normal_distribution vertical_dist{0.0, 1.0};
	std::uniform_real_distribution angle_dist{0.0, 2 * PI};

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

static std::vector<Particle> &scale(std::vector<Particle> &particles, myfloat x, myfloat y, myfloat z)
{
	myvec3 scale{x, y, z};
	for (auto &p : particles)
	{
		p.p *= scale;
	}
	return particles;
}

static std::vector<Particle> &translate(std::vector<Particle> &particles, myfloat x, myfloat y, myfloat z)
{
	myvec3 scale{x, y, z};
	for (auto &p : particles)
	{
		p.p += scale;
	}
	return particles;
}

static std::vector<Particle> &add_velocity(std::vector<Particle> &particles, myfloat x, myfloat y, myfloat z)
{
	myvec3 scale{x, y, z};
	for (auto &p : particles)
	{
		p.v += scale;
	}
	return particles;
}

static std::vector<Particle> &add_angular_momentum(std::vector<Particle> &particles, myvec3 axis)
{
	for (auto &p : particles)
	{
		p.v += glm::cross(p.p, axis);
	}
	return particles;
}

static std::vector<Particle> &add_random_velocity(std::vector<Particle> &particles, myfloat x, myfloat y, myfloat z)
{
	myvec3 scale{x, y, z};
	for (auto &p : particles)
	{
		myfloat x, y, z;
		x = normal_dist(mt);
		y = normal_dist(mt);
		z = normal_dist(mt);
		myvec3 velo{x, y, z};

		p.v += scale * velo;
	}
	return particles;
}

static std::vector<Particle> &set_mass(std::vector<Particle> &particles, myfloat m)
{
	for (auto &p : particles)
	{
		p.m = m;
	}
	return particles;
}

static std::vector<Particle> &set_maxwell_v_dist(std::vector<Particle> &particles, myfloat temperature = 0.1)
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

static std::vector<Particle> &add_radial_velocity(std::vector<Particle> &particles, myfloat scale)
{
	for (auto &p : particles)
	{
		p.v = p.p * scale;
	}
	return particles;
}

/*******************************************************************/

static std::vector<Particle> universe2()
{
	auto initial_dist = sphere_distribution(50);
	scale(initial_dist, 100, 100, 100);
	set_mass(initial_dist, 100);

	add_angular_momentum(initial_dist, myvec3(0, .25, 0));
	initial_dist.push_back(Particle{{0, 0, 0}, {0, 0, 0}, 1e4});
	// add_random_velocity(initial_dist, 10, 10, 10);
	return initial_dist;
}

static std::vector<Particle> universe3()
{
	auto initial_dist = normal_distribution(100);
	scale(initial_dist, 80, 80, 80);
	set_mass(initial_dist, 100);

	// add_angular_momentum(initial_dist, pvec3(0, .25, 0));
	// initial_dist.push_back(particle{ 1e4, {0,0,0}, {0,0,0} });
	// add_random_velocity(initial_dist, 10, 10, 10);
	return initial_dist;
}

static std::vector<Particle> universe1()
{
	auto initial_dist = exponential_disk_distribution(6);
	const myfloat diameter = 100;

	scale(initial_dist, diameter, diameter, 10); // flat disk
	set_mass(initial_dist, 200);

	//add_angular_momentum(initial_dist, myvec3(0, .0, .1) / diameter);
	return initial_dist;
}


static std::vector<Particle> universe4()
{
	auto initial_dist = exponential_disk_distribution(6000);
	const myfloat diameter = 100;

	scale(initial_dist, diameter, diameter, diameter/10); // flat disk
	set_mass(initial_dist, 10);

	add_angular_momentum(initial_dist, myvec3(0, .0, .1) / diameter);
	return initial_dist;
}


static std::vector<Particle> bigbang()
{
	auto initial_dist = sphere_distribution(1000);

	const myfloat diameter = 0.1;
	scale(initial_dist, diameter, diameter, diameter);
	add_radial_velocity(initial_dist, 1e3);

	set_mass(initial_dist, 50);

	return initial_dist;
}

static std::vector<Particle> stable_orbit()
{
	Particle p1 = {{0, 0, 100}, {0, 10, 0}, 0};
	Particle p2 = {{0, 0, 0}, {0, 0, 0}, 10000};
	Particle p3 = {{0, 0, -100}, {0, -10, 0}, 0};

	std::vector<Particle> particles{p1, p2, p3};

	return particles;
}