#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Distributions.h"
#include "Order.h"
#include "Particle.h"
#include "Simulation.h"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  bool output_csv = false;
  bool brute_force = false;
  bool no_reorder = true;
  int num_particles = 1000;
  myfloat theta = 1.5;
  double duration = 10;
  myfloat timestep = .1;

  po::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "num_particles,n", po::value<int>(&num_particles),
      "number of particles in the simulation. Default is 1000")(
      "theta,T", po::value<myfloat>(&theta),
      "Multipole rejection threshold for barnes-hut. Default is 1.5")(
      "duration,d", po::value<myfloat>(&duration),
      "Simulation duration. Default is 10s")(
      "timestep,t", po::value<myfloat>(&timestep),
      "Simulation timestep. Default is 0.1s")(
      "csv", po::bool_switch(&output_csv),
      "output the results into a csv file")(
      "brute_force", po::bool_switch(&brute_force),
      "Enable this flag to use the brute force algorithm")(
      "noreorder", po::bool_switch(&no_reorder),
      "Reorders the particle array before starting the calculation to improve "
      "cache locality. This is only useful for the barnes-hut algorithm.");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  set_seed(42);
  std::vector<Particle> particles =
      make_universe(Distribution::UNIVERSE4, num_particles);

  if (!no_reorder && !brute_force)
    computeAndOrder(particles);

  if (output_csv) {
    // Create a CSV file for the particles and generate the header
    std::ofstream csvfile;
    csvfile.open("output.csv");
    simulate(particles, duration, timestep, &csvfile, brute_force, theta);
    csvfile.close();
  } else {
    simulate(particles, duration, timestep, nullptr, brute_force, theta);
  }

  // Check for nan and inf, throw if encountered
  for (const auto &p : particles) {
    if (std::isnan(p.p.x) || std::isnan(p.p.y) || std::isnan(p.p.z)) {
      throw std::runtime_error("Nan encountered");
    }
    if (std::isinf(p.p.x) || std::isinf(p.p.y) || std::isinf(p.p.z)) {
      throw std::runtime_error("Inf encountered");
    }
  }
}