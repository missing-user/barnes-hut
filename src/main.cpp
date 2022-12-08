#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

#include "Simulation/simulation.h"
#include "Simulation/distributions.h"

#include "QuadTree/QuadTree.h"

int main()
{
    std::vector<Particle> particles = universe4(); // stable_orbit(); // universe1();
    
    // Create a CSV file for the particles and generate the header
    std::ofstream csvfile;
    csvfile.open("output.csv");

    // Run the simulation
    simulate(particles, 10, 0.1, root, &csvfile);

    // Print the resulting values of all particles
    for (auto &p : particles)
    {
        std::cout << p << "\n";
    }
    csvfile.close();
}