#ifndef WEIGHT_TREE
#define WEIGHT_TREE

#include <memory>
#include <vector>

#include "Cuboid.h"
#include "Forces.h"
#include "Particle.h"

class Tree
{
protected:
  Cuboid cuboid;
  int level;  
  CenterOfMass COM;
  
  bool leaf;
  myvec3 divisor;

  std::vector<Tree> branches;
  std::vector<std::shared_ptr<Particle>> particles;

  bool lessThanTheta(const myvec3 &pos, double theta) const;
  void createBranches();
  void createBranches(const myvec3 &pos);
  int selectOctant(const myvec3 &pos) const;
  CenterOfMass computeCOM();

  void insert(std::shared_ptr<Particle> p);
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
