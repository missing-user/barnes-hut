#include <boost/program_options.hpp>
#include <iostream>

#include "Distributions.h"
#include "Simulation.h"
#include "OutputWriter.h"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  bool output_csv = false;
  bool output_pcd = false;
  bool brute_force = false;
  int num_particles = 1000;
  int universe = static_cast<int>(Distribution::UNIVERSE4);
  double theta = 1.5;
  double duration = 10;
  double timestep = .1;

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce this help message")(
      "num_particles,n", po::value<int>(&num_particles)->default_value(1000),
      "number of particles in the simulation.")(
      "theta,T", po::value<double>(&theta)->default_value(1.5),
      "Multipole rejection threshold for barnes-hut.")(
      "duration,d", po::value<double>(&duration)->default_value(10),
      "Simulation duration.")(
      "timestep,t", po::value<double>(&timestep)->default_value(0.1),
      "Simulation timestep.")(
      "csv", po::bool_switch(&output_csv)->default_value(false),
      "output the results into a .csv file")(
      "bin", po::bool_switch(&output_pcd)->default_value(false),
      "output the results into a .particles file (binary format)")(
      "brute_force", po::bool_switch(&brute_force)->default_value(false),
      "Enable this flag to use the brute force algorithm")(
      "universe,u", po::value<int>(&universe)->default_value(3),
      "Universe to simulate. UNIVERSE4=3, COLLISION=4, BIGBANG=5, CRYSTALLINE=8");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  set_seed(42);
  auto particles = make_universe(static_cast<Distribution>(universe), num_particles);

  if (output_csv)
  {
    simulate(particles, duration, timestep, brute_force, theta, writeToCsvFile);
  }else if(output_pcd){
    simulate(particles, duration, timestep, brute_force, theta, writeToBinaryFile);
  }else{
    simulate(particles, duration, timestep, brute_force, theta);
  }
}