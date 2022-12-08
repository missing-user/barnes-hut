# Usage

create a build folder at the top level, navigate into it.
From the build folder run: 
- cmake ..
- make
- make tests

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
- [ ] Implement the Barnes Hut algorithm for simulating the system of bodies over a specified duration of time.
- [ ] Create a unit test, that compares the brute force reference solution with the Barnes-Hut approximation for a small test dataset
- [x] Output a timeseries of the resulting positions of all bodies into a file (e.g. .csv with timestamps) that can be visualized with external tools (e.g. Python script)

## Sprint 2 (OOP)

In this sprint, we will add visualization capabilities to the project, allow the user to configure simulation settings via a config file or command line parameters. The code will be restructured in an Object Oriented way such that the data structure will be based on classes. The functions will be abstracted with interfaces such that they can be switched to any arbitrary interaction function. 

### Sprint 2: Definition of "done"

- Add the ability to specify parameters like the initial conditions, timestep size, simulation duration etc. in either a configuration file or as command line parameters.
- Create a class/datastructure that contains all the parameters (mass, position, velocity) of each body.
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


