#ifndef TREE
#define TREE

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include "Cuboid.h"
#include "Particle.h"

class Tree
{
    int maxPnts = 8;
    Particle COM;

public:
    std::vector<Tree> branches = {};
    std::vector<std::shared_ptr<Particle>> particles = {};
    int level = 0;
    const Cuboid cuboid;
    bool leaf = true;

    Tree(const Cuboid &cuboidIn);

    void createBranches();

    void insertPnt(const Particle &p);

    Particle computeCOM();

    myvec3 computeAcc(const Particle &particle, myfloat theta) const;

    bool less_than_theta(const Particle &particle, double theta) const;

    std::string print() const;
};
#endif
