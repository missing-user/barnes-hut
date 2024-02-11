#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include "OutputWriter.h"

void writeToCsvFile(const Particles &particles,
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
  for (size_t i = 0; i < particles.size(); i++) {
    csvfile << particles.p.x[i] << "," << particles.p.y[i] << "," << particles.p.z[i] << "," << particles.m[i] << "\n";
  }
  csvfile.close();
  //}).detach();
}

void writeToBinaryFile(const Particles &particles,
                    const std::size_t timestep) {
  std::string filename = "barnes_hut_"+std::to_string(timestep)+".particles";
  // Open a Binary file and write the positions and velocities of all particles
  std::ofstream vtkfile;
  if (vtkfile.is_open()) {
    std::cerr << "Could not open file " << filename << " for writing\n";
  }
  vtkfile.open(filename, std::ios::binary);
  std::vector<float> output;
  for (size_t i = 0; i < particles.size(); i++) {
    output.push_back(particles.p.x[i]);
    output.push_back(particles.p.y[i]);
    output.push_back(particles.p.z[i]);
    output.push_back(particles.m[i]);
  }

  vtkfile.write(reinterpret_cast<const char *>(output.data()), sizeof(float)*output.size());
  vtkfile.close();

}