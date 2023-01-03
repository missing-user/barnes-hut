#ifndef TREE
#define TREE

#include <memory>
#include <string>
#include <vector>

#include "Cuboid.h"
#include "Particle.h"

class Tree {
  static int maxPointsPerNode;

private:
  Particle COM;
  bool less_than_theta(const Particle &particle, double theta) const;
  void createBranches();
  std::vector<Tree> branches = {};
  bool leaf = true;

public:
  std::vector<std::shared_ptr<Particle>> particles = {};
  const int level;
  const Cuboid cuboid;

  Tree(const Cuboid &cuboidIn);
  Tree(const Cuboid &cuboidIn, int levelIn);
  Tree(const std::vector<Particle> &particles);

  void insertPnt(const Particle &p);

  Particle computeCOM();

  myvec3 computeAcc(const Particle &particle, myfloat theta) const;

  std::string print() const;

  static void setMaxPointsPerNode(int maxPointsPerNodeIn) {
    maxPointsPerNode = maxPointsPerNodeIn;
  };
};
#endif
