#ifndef WEIGHT_TREE
#define WEIGHT_TREE

#include <iostream>
#include <memory>
#include <vector>

#include "Cuboid.h"
#include "Forces.h"
#include "Particle.h"

class Tree
{
protected:
  CenterOfMass COM;
  std::vector<Tree> branches;
  std::vector<std::unique_ptr<Particle>> particles;

  const Cuboid cuboid;
  const int level;  
  
  myvec3 divisor;
  bool leaf;

  bool lessThanTheta(const myvec3 &pos, double theta) const;
  void createBranches();
  void createBranches(const myvec3 &pos);
  int selectOctant(const myvec3 &pos) const;
  CenterOfMass computeCOM();

  void insert(std::unique_ptr<Particle> p);
  void subdivide();
  void print() const;

public:
  static int maxDepth;
  static int maxParticles;

  Tree(const Cuboid &cuboidIn, int levelIn);
  Tree(const std::vector<Particle> &particles);

  std::vector<std::size_t> DFS() const;
  std::pair<int, int> MaxDepthAndParticles() const;
  std::vector<DrawableCuboid> GetBoundingBoxes() const;

  myvec3 computeAccFromPos(const myvec3 &pos, myfloat theta) const;
  myvec3 computeAcc(const Particle &p1, myfloat theta) const;
};
#endif
