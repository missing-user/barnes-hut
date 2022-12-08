#include <gtest/gtest.h>
#include "QuadTree/QuadTree.h"
#include "Simulation/simulation.h"
#include "Simulation/distributions.h"

#include <fstream>
#include <fstream>

TEST(Simulation, Analytical100Rotations) {
	// Stable orbit solution
	// v = sqrt(G*m/r)
	Particle p1 = {{0,0,100}, {0,10,0},0 };
	Particle p2 = {{0,0,0}, {0,0,0}, 10000};
	Particle p3 = {{0,0,-100}, {0,-10,0},0 };

	std::vector<Particle> particles;
	particles.push_back(p1);
	particles.push_back(p2);
	particles.push_back(p3);

	//Simulate 10 rotations around a massive object + one quarter rotation
	const double simDuration = 10 * p1.p.z * (2 * PI) / p1.v.y;
	const double THRESHOLD = 0.9;

	simulate(particles, simDuration, 0.05);

	// Particle 1 should return to its initial position 
	EXPECT_DOUBLE_EQ(particles[0].p.x, 0);
	EXPECT_NEAR(particles[0].p.y, p1.p.y, THRESHOLD);
	EXPECT_NEAR(particles[0].p.z, p1.p.z, THRESHOLD) << "diff: " << (particles[0].p - p1.p) << std::endl;

	// Particle 2 should not move at all
	EXPECT_DOUBLE_EQ(particles[1].p.y, 0);
	EXPECT_DOUBLE_EQ(particles[1].p.z, 0);
}

TEST(Simulation, AnalyticalQuarterRotation) {
	//Stable orbit solution
	// v = sqrt(G*m/r)
	Particle p1 = {{0,0,10}, {0,5,0}, 0};
	Particle p2 = {{0,0,0}, {0,0,0}, 250};

	std::vector<Particle> particles;
	particles.push_back(p1);
	particles.push_back(p2);

	// Simulate 100 rotations around a massive object + one quarter rotation
	// time = radius * pi/2 / (orbit velocity)
	const double simDuration = p1.p.z * (PI / 2) / p1.v.y;
	const double THRESHOLD = 0.8;

	simulate(particles, simDuration, 0.05);

	// Particle 1 should have completed a 90deg rotation 
	// (y and z coordinates have switched places)
	EXPECT_DOUBLE_EQ(particles[0].p.x, 0);
	EXPECT_NEAR(particles[0].p.z, p1.p.y, THRESHOLD);
	EXPECT_NEAR(particles[0].p.y, p1.p.z, THRESHOLD);

	// Particle 2 should not move at all
	EXPECT_DOUBLE_EQ(particles[1].p.y, 0);
	EXPECT_DOUBLE_EQ(particles[1].p.z, 0);
}


TEST(Simulation, DifferentTimesteps) {

	Particle p1 = {{0,0,0}, {1.99255722618,0.95797067952,0.4032403082}, 20};
	std::vector<Particle> particles(1);
	auto simDuration = 1;

	// Step through i orders of magnitude for the step size, starting at 1
	// the irrational base was chosen to test the residual timestepping capabilities
	for (double i = 1; i > -13; i--)
	{
		particles[0] = p1;
		auto timestep = std::exp(i);
		simulate(particles, simDuration, timestep);
		EXPECT_FLOAT_EQ(particles[0].p.x, p1.v.x * simDuration) << " at timestep size " << timestep;
		EXPECT_FLOAT_EQ(particles[0].p.y, p1.v.y * simDuration) << " at timestep size " << timestep;
		EXPECT_FLOAT_EQ(particles[0].p.z, p1.v.z * simDuration) << " at timestep size " << timestep;
	}
}

TEST(FileWriting, IsValidCsv) {
	Particle p1 = {{5,6,7},{1,2,3}, 20};
	std::vector<Particle> particles{ p1 };
	auto simDuration = 1;

	std::stringstream buffer;
	simulate(particles, simDuration, 0.5, &buffer);

	std::string expectedFileContent = "px_0,py_0,pz_0,time\n5.5,7,8.5,0.5\n6,8,10,1\n";
	EXPECT_EQ(buffer.str().compare(expectedFileContent), 0) << "The real file content was:\n" << buffer.str();
}

TEST(BarnesHut, CompareTheta0) {

	std::vector<Particle> particles = universe1();
    std::vector<Particle> particlesTree{particles};

	const auto simDuration = 10.0;
	const auto timestep = 0.1;


	simulate(particles, simDuration, timestep);

	simulate(particlesTree, simDuration, timestep, nullptr, false, 0);

	for (int i = 0; i<particles.size(); i++) {
		EXPECT_FLOAT_EQ(particles[i].p.x, particlesTree[i].p.x);
		EXPECT_FLOAT_EQ(particles[i].p.y, particlesTree[i].p.y);
		EXPECT_FLOAT_EQ(particles[i].p.z, particlesTree[i].p.z);	
	}
}

TEST(BarnesHut, CompareApproximation) {

	std::vector<Particle> particles = universe1(); 
	std::vector<Particle> particlesTree{ particles };

	const auto simDuration = 10.0;
	const auto timestep = 0.1;

	simulate(particles, simDuration, timestep);

	simulate(particlesTree, simDuration, timestep, nullptr, false, 1.5);

	for (int i = 0; i < particles.size(); i++) {
		// The approximate values should not be identical to the real ones
		EXPECT_NE(particles[i].p.x, particlesTree[i].p.x);
		EXPECT_NE(particles[i].p.y, particlesTree[i].p.y);
		EXPECT_NE(particles[i].p.z, particlesTree[i].p.z);

		// But they should be pretty close
		EXPECT_NEAR(particles[i].p.x, particlesTree[i].p.x, 0.1);
		EXPECT_NEAR(particles[i].p.y, particlesTree[i].p.y, 0.1);
		EXPECT_NEAR(particles[i].p.z, particlesTree[i].p.z, 0.1);
	}
}
