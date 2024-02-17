#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Distributions.h"
#include "Order.h"
#include "Particle.h"
#include "Simulation.h"
#include "OutputWriter.h"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  bool output_csv = false;
  bool output_pcd = false;
  bool brute_force = false;
  int num_particles = 1000;
  double theta = 1.5;
  double duration = 10;
  double timestep = .1;

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce this help message")(
      "num_particles,n", po::value<int>(&num_particles),
      "number of particles in the simulation. (default 1000)")(
      "theta,T", po::value<double>(&theta),
      "Multipole rejection threshold for barnes-hut. (default 1.5)")(
      "duration,d", po::value<double>(&duration),
      "Simulation duration. (default 10s)")(
      "timestep,t", po::value<double>(&timestep),
      "Simulation timestep. (default 0.1s)")(
      "csv", po::bool_switch(&output_csv),
      "output the results into a .csv file")(
      "bin", po::bool_switch(&output_pcd),
      "output the results into a .particles file (binary format)")(
      "brute_force", po::bool_switch(&brute_force),
      "Enable this flag to use the brute force algorithm");
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
  
  if (output_csv)
  {
    simulate(particles, duration, timestep, brute_force, theta, writeToCsvFile);
  }else if(output_pcd){
    simulate(particles, duration, timestep, brute_force, theta, writeToBinaryFile);
  }else{
    simulate(particles, duration, timestep, brute_force, theta);
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