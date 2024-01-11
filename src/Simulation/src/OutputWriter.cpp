#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include "OutputWriter.h"

void writeToCsvFile(const std::vector<Particle> &particles,
                    const std::size_t timestep) { 
  //std::thread([&filename, particles](){ 
  std::string filename = "output"+std::to_string(timestep)+".csv";
  // Open a csv file and write the positions of all particles
  std::ofstream csvfile;
  if (csvfile.is_open()) {
    std::cerr << "Could not open file " << filename << " for writing\n";
  }
  csvfile.open(filename);
  csvfile << "x,y,z,m\n";
  csvfile << particles;
  csvfile.close();
  //}).detach();
}

void writeToBinaryFile(const std::vector<Particle> &particles,
                    const std::size_t timestep) {
  std::string filename = "barnes_hut_"+std::to_string(timestep)+".particles";
  // Open a Binary file and write the positions and velocities of all particles
  std::ofstream vtkfile;
  if (vtkfile.is_open()) {
    std::cerr << "Could not open file " << filename << " for writing\n";
  }
  vtkfile.open(filename, std::ios::binary);
  std::vector<float> output;
  int i = 0;
  for(const auto& p : particles) {
    output.push_back(p.p.x);
    output.push_back(p.p.y);
    output.push_back(p.p.z);
    output.push_back(p.v.length());
  }
  vtkfile.write(reinterpret_cast<const char *>(output.data()), sizeof(float)*output.size());
  vtkfile.close();

}