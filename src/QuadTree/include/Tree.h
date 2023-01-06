#ifndef TREE
#define TREE

#include <memory>
#include <string>
#include <vector>

#include "Cuboid.h"
#include "Particle.h"

class Tree {
  static int maxDepth;
  static int maxParticles;

private:
  Particle COM;
  bool less_than_theta(const myvec3 pos, double theta) const;
  void createBranches();

  int selectOctant(const myvec3 pos) const;
  std::vector<Tree> branches;
  bool leaf = true;
  const int level;
  const Cuboid cuboid;

public:
  std::vector<std::shared_ptr<Particle>> particles;

  Tree(const Cuboid &cuboidIn);
  Tree(const Cuboid &cuboidIn, int levelIn);
  Tree(const std::vector<Particle> &particles);

  void insertPnt(const Particle &p);

  Particle computeCOM();

  myvec3 computeAcc(const Particle &particle, myfloat theta) const;

  std::string print() const;

  static void setMaxDepth(int maxDepthIn) { maxDepth = maxDepthIn; };
};
#endif
