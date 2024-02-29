#include "Distributions.h"
#include "Simulation.h"
#include "Barneshut.h"
#include "OutputWriter.h"
#include <gtest/gtest.h>

#include <fstream>
#include <numbers>
#include <chrono>
#include <thread>
#include "Order.h"

TEST(Simulation, Analytical45Rotations) {
  // Simulate 4.5 rotations of two massless particles around a massive object.
  // After the simulation, both particles should have switched places.
  
  // Stable orbit solution
  // v = sqrt(G*m/r)
  Particles p0 = make_universe(Distribution::STABLE_ORBIT, 3);
  Particles particles = make_universe(Distribution::STABLE_ORBIT, 3);

  const double simDuration = 4.5 * p0.p.z[0] * (2 * std::numbers::pi) / p0.v.y[0];

  simulate(particles, simDuration, 0.1);

  for (int component = 0; component < 7; component++)
  {
    // Particles should have switched places
    EXPECT_NEAR(particles.get(component)[0], p0.get(component)[2], .12);
    EXPECT_NEAR(particles.get(component)[2], p0.get(component)[0], .12);
    // The heavy Particle 2 should not move at all
    EXPECT_DOUBLE_EQ(particles.get(component)[1], p0.get(component)[1]);
  }
}

TEST(Simulation, DifferentTimesteps) {
  // When simulating a certain time duration, the results should be identical,
  // regardless of the timestep size that was used.

  Particles particles{1};
  particles.m[0] = 20;
  particles.p.x[0] = 0;
  particles.p.y[0] = 0;
  particles.p.z[0] = 0;
  particles.v.x[0] = 1.99255722618;
  particles.v.y[0] = 0.95797067952;
  particles.v.z[0] = 0.4032403082;
  auto simDuration = 1.52;

  // Step through i orders of magnitude for the step size, starting at 1
  // the irrational base was chosen to test the residual timestepping
  // capabilities
  for (double i = 1; i > -12; i--) {
    particles.m[0] = 20;
    particles.p.x[0] = 0;
    particles.p.y[0] = 0;
    particles.p.z[0] = 0;
    particles.v.x[0] = 1.99255722618;
    particles.v.y[0] = 0.95797067952;
    particles.v.z[0] = 0.4032403082;
    auto timestep = std::exp(i);
    simulate(particles, simDuration, timestep);
    EXPECT_FLOAT_EQ(particles.p.x[0], particles.v.x[0] * simDuration)
        << " at timestep size " << timestep;
    EXPECT_FLOAT_EQ(particles.p.y[0], particles.v.y[0] * simDuration)
        << " at timestep size " << timestep;
    EXPECT_FLOAT_EQ(particles.p.z[0], particles.v.z[0] * simDuration)
        << " at timestep size " << timestep;
  }
}

TEST(QuadTree, DepthCalculation) {
  // Create a distribution of particles and check if the depth of the tree is as
  // expected
  set_seed(4756);
  auto particles = make_universe(Distribution::UNIVERSE4, 800);
  auto info = bh_superstep_debug({100,0,0}, particles, 1.5*1.5);
  EXPECT_EQ(info.depth, 7);
  EXPECT_EQ(info.max_particles_in_leaf, 4);
  auto force_calcs_15 = info.debug_boxes.size();
  
  info = bh_superstep_debug({100,0,0}, particles, 0.7*0.7);
  EXPECT_EQ(info.depth, 7);
  EXPECT_EQ(info.max_particles_in_leaf, 4);
  auto force_calcs_07 = info.debug_boxes.size();
  
  // Count force evaluations, increasing theta should reduce computational cost
  EXPECT_GT(force_calcs_15, 100);
  EXPECT_GT(force_calcs_07, force_calcs_15);
  EXPECT_LT(force_calcs_07, 800);
}

TEST(FileWriting, IsValidCsv) {
  // Validate our output file format

  Particles particles{1};
  particles.p.x[0] = 5;
  particles.p.y[0] = 6;
  particles.p.z[0] = 7;
  particles.v.x[0] = 1;
  particles.v.y[0] = 2;
  particles.v.z[0] = 3;
  particles.m[0] = 20;

  simulate(particles, 2, 0.5, true, 1.5, writeToCsvFile);

std::vector<std::string> files{"output0.csv", "output1.csv","output2.csv", "output3.csv"};
std::vector<std::string> expectedFileContents{
    "x,y,z,m\n5,6,7,20\n",
    "x,y,z,m\n5.5,7,8.5,20\n",
    "x,y,z,m\n6,8,10,20\n",
    "x,y,z,m\n6.5,9,11.5,20\n",
    "x,y,z,m\n7,10,13,20\n"};

  for(int i = 0; i<files.size(); i++){
    std::ifstream f(files[i]);
    // Wait for the file to be written, maximum 1 second
    for(int i = 0; i<10; i++){
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      if(f.good() && !f.is_open()){
        break;
      }
    }
    EXPECT_TRUE(f.good());

    // Read the file into a string
    std::string buffer;
    buffer.assign(std::istreambuf_iterator<char>(f),
                  std::istreambuf_iterator<char>());

    EXPECT_EQ(buffer.compare(expectedFileContents[i]), 0)
        << "The real file content was:\n"
        << buffer<<"\nbut we expected:\n"<<expectedFileContents[i]<<"\n";
    f.close();
  }
}

TEST(BarnesHut, CompareTheta0) {
  // With the multipole rejection parameter set to 0, the barnes hut algorithm
  // degenerates to a brute force solution. The results should be identical.

  // Fix the random seed, so test cases are reproducible
  set_seed(4756);
  auto particles = make_universe(Distribution::UNIVERSE1, 5);
  set_seed(4756);
  auto particlesTree = make_universe(Distribution::UNIVERSE1, 5);

  const auto simDuration = 10.0;
  const auto timestep = 0.1;

  for (int i = 0; i < particles.size(); i++) {
    EXPECT_FLOAT_EQ(particles.p.x[i], particlesTree.p.x[i]);
    EXPECT_FLOAT_EQ(particles.p.y[i], particlesTree.p.y[i]);
    EXPECT_FLOAT_EQ(particles.p.z[i], particlesTree.p.z[i]);
  }

  simulate(particles, simDuration, timestep, true, 0.0);
  simulate(particlesTree, simDuration, timestep, false, 0.0);

  // Since Barnes Hut reorders the particles, we also need to reorder the 
  // brute force results to compare them. 
  auto boundingbox = bounding_box(particles.p);
  computeAndOrder(particles, boundingbox);

  for (int i = 0; i < particles.size(); i++) {
    EXPECT_FLOAT_EQ(particles.p.x[i], particlesTree.p.x[i]);
    EXPECT_FLOAT_EQ(particles.p.y[i], particlesTree.p.y[i]);
    EXPECT_FLOAT_EQ(particles.p.z[i], particlesTree.p.z[i]);
  }
}

TEST(BarnesHut, CompareApproximation) {
  // Check if Barnes Hut approximation is working correctly.
  // The time evolution of the particles should be similiar, but not identical
  // to the brute force solution. Does not work at large timescales, because the
  // system is chaotic in nature.

  // Fix the random seed, so test cases are reproducible
  set_seed(4756);
  auto particles = make_universe(Distribution::UNIVERSE1, 64);
  set_seed(4756);
  auto particlesTree = make_universe(Distribution::UNIVERSE1, 64);
  
  const auto simDuration = 1.0;
  const auto timestep = 0.1;

  simulate(particles, simDuration, timestep);
  simulate(particlesTree, simDuration, timestep, false, 1.5, writeToCsvFile);

  auto boundingbox = bounding_box(particles.p);
  computeAndOrder(particles, boundingbox);

  for (int i = 0; i < particles.size(); i++) {
    // The approximate values should not be identical to the real ones
    EXPECT_NE(particles.p.x[i], particlesTree.p.x[i]);
    EXPECT_NE(particles.p.y[i], particlesTree.p.y[i]);
    EXPECT_NE(particles.p.z[i], particlesTree.p.z[i]);

    // But they should be pretty close
    EXPECT_NEAR(particles.p.x[i], particlesTree.p.x[i], .25);
    EXPECT_NEAR(particles.p.y[i], particlesTree.p.y[i], .2);
    EXPECT_NEAR(particles.p.z[i], particlesTree.p.z[i], .1);
  }
}
