[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitLab Release](https://img.shields.io/badge/Release-Sprint%201-red)](https://gitlab.lrz.de/advprog2022/13/barnes-hut/-/tree/version2)
![Test Result Badge](https://gitlab.lrz.de/advprog2022/13/barnes-hut/badges/version2/pipeline.svg)

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
5. Runs the python script for visualizing the behaviour of the system of particles over a set period of time in 3D space with dynamic POV

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
- (Optional) Try Feedback-Directed Compiler Optimization (FDO) in g++ and MSVC <https://learn.microsoft.com/en-us/cpp/build/profile-guided-optimizations?view=msvc-170>

- [x] Switched from shared_ptr to unique_ptr implementation for the tree (Improved tree build times by slightly)
- [x] Only leaf nodes store a vector of pointers, reducing the number of stored particle pointers from Nlog(N) to N
- [x] Sorting the particle array before the simulation to improve cache coherence (In single threaded benchmark UNIVERSE4 with 30k particles, 10s and 0.1s timestep the execution time was reduced from 1m13.605 to 55.010s. **That is a 30% improvement in execution time**) (Valgrind Cache miss report reports more cache misses tho, which is weird):
w/o sorting:
==8270== I   refs:      15,923,885,607
==8270== I1  misses:             2,789
==8270== LLi misses:             2,741
==8270== I1  miss rate:           0.00%
==8270== LLi miss rate:           0.00%
==8270==
==8270== D   refs:       7,372,560,122  (4,543,357,564 rd   + 2,829,202,558 wr)
==8270== D1  misses:       457,093,410  (  454,439,822 rd   +     2,653,588 wr)
==8270== LLd misses:         2,343,432  (    1,661,398 rd   +       682,034 wr)
==8270== D1  miss rate:            6.2% (         10.0%     +           0.1%  )
==8270== LLd miss rate:            0.0% (          0.0%     +           0.0%  )
==8270==
==8270== LL refs:          457,096,199  (  454,442,611 rd   +     2,653,588 wr)
==8270== LL misses:          2,346,173  (    1,664,139 rd   +       682,034 wr)
==8270== LL miss rate:             0.0% (          0.0%     +           0.0%  )

w sorting:
==8536==
==8536== I   refs:      16,025,914,742
==8536== I1  misses:             2,874
==8536== LLi misses:             2,839
==8536== I1  miss rate:           0.00%
==8536== LLi miss rate:           0.00%
==8536==
==8536== D   refs:       7,422,843,756  (4,573,795,704 rd   + 2,849,048,052 wr)
==8536== D1  misses:       452,273,941  (  449,099,081 rd   +     3,174,860 wr)
==8536== LLd misses:         3,035,603  (    1,981,384 rd   +     1,054,219 wr)
==8536== D1  miss rate:            6.1% (          9.8%     +           0.1%  )
==8536== LLd miss rate:            0.0% (          0.0%     +           0.0%  )
==8536==
==8536== LL refs:          452,276,815  (  449,101,955 rd   +     3,174,860 wr)
==8536== LL misses:          3,038,442  (    1,984,223 rd   +     1,054,219 wr)
==8536== LL miss rate:             0.0% (          0.0%     +           0.0%  )

with debug symbols:
w/o sorting:
==9184==
==9184== I   refs:      2,661,285,937
==9184== I1  misses:            2,800
==9184== LLi misses:            2,696
==9184== I1  miss rate:          0.00%
==9184== LLi miss rate:          0.00%
==9184==
==9184== D   refs:      1,441,670,368  (916,998,110 rd   + 524,672,258 wr)
==9184== D1  misses:       49,184,166  ( 48,825,583 rd   +     358,583 wr)
==9184== LLd misses:           81,154  (      8,037 rd   +      73,117 wr)
==9184== D1  miss rate:           3.4% (        5.3%     +         0.1%  )
==9184== LLd miss rate:           0.0% (        0.0%     +         0.0%  )
==9184==
==9184== LL refs:          49,186,966  ( 48,828,383 rd   +     358,583 wr)
==9184== LL misses:            83,850  (     10,733 rd   +      73,117 wr)
==9184== LL miss rate:            0.0% (        0.0%     +         0.0%  )

w sorting:
==8958==
==8958== I   refs:      2,680,896,810
==8958== I1  misses:            2,904
==8958== LLi misses:            2,734
==8958== I1  miss rate:          0.00%
==8958== LLi miss rate:          0.00%
==8958==
==8958== D   refs:      1,452,245,492  (923,634,238 rd   + 528,611,254 wr)
==8958== D1  misses:       48,611,530  ( 48,199,351 rd   +     412,179 wr)
==8958== LLd misses:          142,267  (      8,039 rd   +     134,228 wr)
==8958== D1  miss rate:           3.3% (        5.2%     +         0.1%  )
==8958== LLd miss rate:           0.0% (        0.0%     +         0.0%  )
==8958==
==8958== LL refs:          48,614,434  ( 48,202,255 rd   +     412,179 wr)
==8958== LL misses:           145,001  (     10,773 rd   +     134,228 wr)
==8958== LL miss rate:            0.0% (        0.0%     +         0.0%  )
)

- [x] Checking if the distance is zero in the force computation is more expensive than just computing and subtracting the softened force (always compute and subtract: 50k 52s 9.9s vs check if distance is zero: 50k 54s 10.5s) 5k: 6.8s 6.9s vs 6.7s 6.6s
- [x] Is passing the position vector by value faster than by reference in the computeAcceleration function? No. (17.3s by value, 10.0s by reference)
- [x] Also compare if by value and by reference make a difference in Tree.cpp selectOctant() and less_than_theta(). No measurable difference. (10s by value, 10s by reference)

- [ ] Currently the leading cause for L1 Cache misses is the `if(leaf)` statement in the Tree traversal. Could this be improved somehow?
- [ ] OpenMP parallel for loop for multithreading.
- [ ] OpenMP SIMD for vectorization.
