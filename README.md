# Barnes-Hut Simulation
This project is a implementation of the Barnes-Hut algorithm for approximating the gravitational forces between objects in a system. It is commonly used in astrophysical simulations to model the motion of celestial bodies, such as planets and stars. The simulation can be initialized with the positions, masses, and velocities of the objects, and will iteratively update the positions and velocities based on the gravitational forces acting on them.

The basic idea behind the Barnes-Hut algorithm is to divide the system into a grid of cells, and to approximate the forces between objects in a cell using the center of mass and total mass of the cell. This allows the simulation to scale to systems with a large number of objects, as the calculation of forces between each individual object would become computationally infeasible. To improve the accuracy of the simulation, the size of the cells can be decreased, which will result in a more precise calculation of the forces between objects. However, this will also increase the computational cost of the simulation.

This project was made for the Advanced Programming course at the Technical University of Munich.
Authors: Phillip Jurašić & Islam Elgamal
Supervision: Gerasimos Chourdakis

References:
Barnes, J., & Hut, P. (1986). A hierarchical O(N log N) force-calculation algorithm. Nature, 324(6096), 446-449.

# Prerequisites 

You need cmake, gtest and glm installed for this project to work. The following commands will install them. The last two lines are for building gtest on your machine.

```sh
apt-get install g++ cmake 
apt-get install libgtest-dev
apt-get install libglm-dev

cd /usr/src/gtest
cmake CMakeLists.txt && make
```

Visualizing the results requires python3 and the following packages: plotly, numpy, pandas. Pandas and numpy are used for reading and transforming the data into a format that can be visualized. Plotly is used for creating the interactive 3D plot.

# Usage

To execute the project, run the following commands:

```sh
mkdir build && cd build && cmake ..

make barnes-hut
./src/barnes-hut

cd ..
pip install -r requirements.txt 
python plot.py
```

1. Creates a folder called build, navigates into it and creates the makefile 
2. Compiles and links the project
3. Runs the executable. This produces a .csv file called output.csv as an output of the simulation
4. Installs all the python dependencies for visualizing the results from requirements.txt
5. Runs the python script for visualizing the results

# Project: Barnes Hut galaxy Simulation 

Idea contributed by Philipp Jurašić and Islam Elgamal

## Motivation

Space is cool, and as admirers of science we want to know how it works. In this project we want to play god and create our own galaxy.

We want to approximate effect of gravity in the formation of gas clouds and clustering of stars over long periods of time. This will be done by simulating the attractive forces between thousands of particles (representing stars and other celestial bodies).

The [Barnes-Hut algorithm](https://en.m.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation#/media/File%3A2D_Quad-Tree_partitioning_of_100_bodies.png) is commonly used in computational physics, specifically to approximately compute n-body interactions.

![](https://upload.wikimedia.org/wikipedia/commons/9/93/2D_Quad-Tree_partitioning_of_100_bodies.png)

## Sprint 1 (basics)

In this sprint we will implement the n-body problem and a basic Barnes-Hut approximation. The program will randomly generate initial conditions for testing and the result will be outputted for visualization. We will compare the result of the brute force "reference" solution with the Barnes-Hut algorithm. 

### Sprint 1: Definition of "done"

- [x] Generate an array of initial masses, positions, and velocities of a system of bodies in 3D space.
- [x] Create a brute force n-body simulation O(n*n) (will be later used as a unit test for verification)
- [x] Create a function for space-dividing an array of coordinates into an octree data structure.
- [x] Implement the Barnes Hut algorithm for simulating the system of bodies over a specified duration of time.
- [x] Create a unit test, that compares the brute force reference solution with the Barnes-Hut approximation for a small test dataset
- [x] Output a timeseries of the resulting positions of all bodies into a file (e.g. .csv with timestamps) that can be visualized with external tools (e.g. Python script)

## Sprint 2 (OOP)

In this sprint, we will add visualization capabilities to the project, allow the user to configure simulation settings via a config file or command line parameters. The code will be restructured in an Object Oriented way such that the data structure will be based on classes. The functions will be abstracted with interfaces such that they can be switched to any arbitrary interaction function. 

### Sprint 2: Definition of "done"

- Add the ability to specify parameters like the initial conditions, timestep size, simulation duration etc. in either a configuration file or as command line parameters.
- Implement appropriate access control modifiers for functions and variables.
- Apply apropriate usage of references and pointers for optimal memory allocation.
- Add the ability to visualize the dynamics of the simulation in 3D space.
- Abstract the interaction function of the bodies (and add an example of how to use it, e.g. gravity potential and coloumb potential) 
- Abstract the space dividing function 
- Clean up and refactor the code

## Sprint 3 (performance and/or STL)

In this sprint, we will analyze and optimize the perfomance and computation time of the program. The focus will be to study how much impact each section of the code has on the total runtime and the effect of each optimization step taken to reduce computation time.

### Sprint 3: Definition of "done"

- Measure how much time is consumed during each section in the code 
- Utilize at least three different optimization techniques and study its impact on total runtime
- At least one function should utilize vectorized instructions 
- (Optional) Try Feedback-Directed Comipler Optimization (FDO) in g++ and MSVC https://learn.microsoft.com/en-us/cpp/build/profile-guided-optimizations?view=msvc-170


