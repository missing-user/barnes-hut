#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

#include "QuadTree/Particle.h"
#include "Simulation/simulation.h"
#include "Simulation/distributions.h"

int main()
{
    std::vector<Particle> particles = universe1(); // stable_orbit(); // universe1();
    
    // Create a CSV file for the particles and generate the header
    std::ofstream csvfile;
    csvfile.open("output.csv");

    // Run the simulation
    simulate(particles, 20, 0.1, &csvfile, false, 0.5);

    // Print the resulting values of all particles
    for (auto &p : particles)
    {
        std::cout << p << "\n";
    }
    csvfile.close();
}