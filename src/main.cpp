#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Distributions.h"
#include "Particle.h"
#include "Simulation.h"

void reorder(std::vector<Particle> &data,
             std::vector<std::size_t> const &order) {
  // Reorder function from
  // https://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices
  std::vector<Particle> tmp; // create an empty vector
  tmp.reserve(data.size());  // ensure memory and avoid moves in the vector
  for (std::size_t i = 0; i < order.size(); ++i) {
    tmp.push_back(data[order[i]]);
  }
  data.swap(tmp); // swap vector contents
}

int main() {
  set_seed(42);
  std::vector<Particle> particles =
      make_universe(Distribution::UNIVERSE4, 30000);

  // Create a CSV file for the particles and generate the header
  std::ofstream csvfile;
  csvfile.open("output.csv");

  // Run the simulation
  // simulate(particles, .1, 0.1, &csvfile, false, 1.5);

  simulate(particles, 10, 0.1, nullptr, false, 1.5);

  // Print the resulting values of all particles
  /*for (const auto &p : particles) {
    std::cout << p << "\n";
  }*/
  csvfile.close();
}