#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Distributions.h"
#include "Particle.h"
#include "Simulation.h"


int main() {
  std::vector<Particle> particles =
      universe2(); // stable_orbit(); // universe1();

  // Create a CSV file for the particles and generate the header
  std::ofstream csvfile;
  csvfile.open("output.csv");

  // Run the simulation
  simulate(particles, 20, 0.1, &csvfile, false, 1.5);

  // Print the resulting values of all particles
  for (const auto &p : particles) {
    std::cout << p << "\n";
  }
  csvfile.close();
}