#include "Tree.h"

Tree::Tree(const Cuboid &cuboidIn) : cuboid(cuboidIn) {}

void Tree::createBranches() // populates the branches array of this object with
                            // new trees from subdividing this node
{
  branches.reserve(
      8); // reserve space for 8 branches to save on resize operations
  for (Cuboid &r : cuboid.subdivide()) {
    Tree branch = Tree(r);
    branch.level = level + 1;
    branches.push_back(branch);
  }
}

void Tree::insertPnt(const Particle &p) // adds a point to the tree structure.
// Depending on the location of the point, new branches could be generated to
// accomodate the point
{
  if (!cuboid.contains(p)) // quits if point is outside of this node's domain
    return;
  particles.push_back(std::make_shared<Particle>(p));
  if (particles.size() >
      maxPointsPerNode) // if maximum number of points in a node is exceeded,
                        // add this point to one of the branches
  {
    if (leaf) {
      leaf = false;
      createBranches(); // If there aren't any branches, create them.
      for (Tree &b : branches) {
        for (auto &pnt : particles) {
          b.insertPnt(*pnt);
        }
      }
    } else {
      for (Tree &b : branches) {
        b.insertPnt(p);
      }
    }
  }
}

Particle Tree::computeCOM() // find the center of mass for this node and save it
                            // as the new COM
{
  Particle tempCOM;
  if (leaf) // if this node doesn't have branches, calculate the center of mass
            // of all contained particles
  {
    if (particles.size() != 0) {
      for (const std::shared_ptr<Particle> P : particles) {
        Particle tempP = *P;
        tempP.p *= tempP.m;
        tempCOM = tempCOM + tempP;
      }
      if (tempCOM.m != 0)
        tempCOM.p /= tempCOM.m;
    }
  } else // if there are branches, calculate the center of mass of the centers
         // of mass of the branches
  {
    for (Tree t : branches) {
      Particle tempP = t.computeCOM();
      tempP.p *= tempP.m;
      tempCOM = tempCOM + tempP;
    }
    if (tempCOM.m != 0)
      tempCOM.p /= tempCOM.m;
  }
  COM = tempCOM;

  std::cout << cuboid.print() << " level: " << level
            << " nPoints:" << particles.size() << " Center of Mass: " << COM
            << ", leaf:" << leaf << std::endl;
  return tempCOM;
}

myvec3 Tree::computeAcc(const Particle &particle, myfloat theta)
    const // compute the accelartion applied on a particle by this node
{
  const myfloat softening_param = 0.025;
  myvec3 acc{0, 0, 0};

  if (leaf) {
    for (const auto sp : particles) {
      if (sp->p !=
          particle.p) // if the queried particle isn't the input particle
      {
        auto diff = sp->p - particle.p;
        auto distance = glm::length(diff);
        acc += diff * sp->m /
               (distance * distance * distance +
                softening_param); // softened gravitational equation
      }
    }
  } else {
    if (less_than_theta(particle, theta)) { // Barnes-Hut threshold
      auto diff = COM.p - particle.p;
      auto distance = glm::length(diff);
      acc = diff * COM.m / (distance * distance * distance + softening_param);
    } else { // if threshold not met, compute the acceleration due to the
             // branches inside this node
      for (const Tree &b : branches) {
        acc += b.computeAcc(particle, theta);
      }
    }
  }

  return acc;

  // }
}

bool Tree::less_than_theta(const Particle &particle, double theta) const {
  // returns whether particle meets the threshold criteria for approximating the
  // acceleration due to this node
  myfloat distance = glm::length(particle.p - COM.p);
  // TODO: this should not be the length of the cuboid
  return glm::length(cuboid.max_extent - cuboid.min_extent) / distance < theta;
}

std::string Tree::print() const {
  std::ostringstream str;
  str << cuboid.print() << " level: " << level
      << " nPoints:" << particles.size() << " Center of Mass: " << COM
      << ", leaf:" << leaf << std::endl;
  for (const Tree b : branches) {
    str << b.print();
  }
  return str.str();
}