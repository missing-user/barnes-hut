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
  Tree* branches;
  std::vector<int> particles;

  Cuboid cuboid;
  int level;  

  bool lessThanTheta(const myvec3 &pos, double theta2) const;
  void createBranches();
  void createBranches(const myvec3 &pos);
  int selectOctant(const myvec3 &pos) const;
  CenterOfMass computeCOM();

  void insert(int p);
  void subdivide();
  void print() const;

public:
  static int maxDepth;
  static int maxParticles;
  static const Particle* firstP;

  Tree();
  Tree(const std::vector<Particle> &particles);
  ~Tree();
  Tree(const Cuboid &cuboidIn, int levelIn);

  std::vector<std::size_t> DFS() const;
  std::pair<int, int> MaxDepthAndParticles() const;
  std::vector<DrawableCuboid> GetBoundingBoxes() const;

  myvec3 computeAccFromPos(const myvec3 &pos, myfloat theta2) const;
  myvec3 computeAcc(const Particle &p1, myfloat theta2) const;
};
#endif
