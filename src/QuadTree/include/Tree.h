#ifndef TREE
#define TREE

#include <iostream>
#include <memory>
#include <vector>

#include "Cuboid.h"
#include "Forces.h"
#include "Particle.h"

class Tree
{
private:
  CenterOfMass COM;
  myvec3 divisor;
  std::vector<Tree> branches;
  std::vector<std::unique_ptr<Particle>> particles;
  bool leaf;

  const int level;
  const Cuboid cuboid;

  bool less_than_theta(const myvec3 &pos, double theta) const;
  void createBranches();
  void createBranchesAtP(const myvec3 &P);
  int selectOctant(const myvec3 &pos) const;
  CenterOfMass computeCOM();
  CenterOfMass computeCOM_NonRecursive();

public:
  static int maxDepth;
  static int maxParticles;

  Tree(const Cuboid &cuboidIn, int levelIn);
  Tree(const std::vector<Particle> &particles);

  void insert(const Particle &p);
  void insert(std::unique_ptr<Particle> p);
  void insertNonRecursive(const Particle &p);
  void insertNonRecursive(std::unique_ptr<Particle> p);
  void subdivideBreadth();
  void subdivideNonRecursive();
  void print() const;

  std::vector<DrawableCuboid> GetBoundingBoxes() const;
  std::pair<int, int> MaxDepthAndParticles() const;

  myvec3 computeAccFromPos(const myvec3 &pos, myfloat theta) const;
  myvec3 computeAcc(const Particle &p1, myfloat theta) const;
  std::vector<std::size_t> DFS() const;
};
#endif
