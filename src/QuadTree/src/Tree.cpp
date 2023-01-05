#include "Tree.h"

#include <iostream>
#include <sstream>

int Tree::maxDepth = 16;

Tree::Tree(const Cuboid &cuboidIn, int levelIn)
    : cuboid(cuboidIn), level(levelIn),
      dimension(glm::length2(cuboid.dimension())) {}

Tree::Tree(const std::vector<Particle> &particles)
    : cuboid(bounding_box(particles)), level(0),
      dimension(glm::length2(cuboid.dimension())) {
  for (const auto &p : particles) {
    insertPnt(p);
  }
  computeCOM();
}

void Tree::createBranches() // populates the branches array of this object with
                            // new trees from subdividing this node
{
  branches.reserve(
      8); // reserve space for 8 branches to save on resize operations
  for (Cuboid &r : cuboid.subdivide()) {
    Tree branch = Tree(r, level + 1);
    branches.push_back(branch);
  }
}

int Tree::selectOctant(const myvec3 pos) const {
  int octant = 0;
  octant += (pos.x > cuboid.center().x) << 0;
  octant += (pos.y > cuboid.center().y) << 1;
  octant += (pos.z > cuboid.center().z) << 2;
  return octant;
}

void Tree::insertPnt(const Particle &p) // adds a point to the tree structure.
// Depending on the location of the point, new branches could be generated to
// accomodate the point
{
  particles.push_back(std::make_shared<Particle>(p));
  if (particles.size() > 1 && level < maxDepth) {
    if (leaf) {
      leaf = false;
      createBranches(); // If there aren't any branches, create them.

      for (auto &pnt : particles) {
        branches[selectOctant(pnt->p)].insertPnt(*pnt);
      }
    } else { // If there are branches, add the point to the appropriate branch
      branches[selectOctant(p.p)].insertPnt(p);
    }
  }
}

Particle Tree::computeCOM() // find the center of mass for this node and save it
                            // as the new COM
{
  COM = {{0, 0, 0}, {0, 0, 0}, 0}; // reset COM to zero
  if (leaf) // if this node doesn't have branches, calculate the center of mass
            // of all contained particles
  {
    for (const auto &p : particles) {
      COM = COM + *p; // Adding two particles creates a new particle which
                      // represents the center of mass of the two particles
    }
  } else {
    for (Tree &t : branches) {
      COM = COM + t.computeCOM(); // This recursively calculates the center of
                                  // mass of each branch and updates them
    }
  }
  // std::cout << COM << " level: " << level << " leaf: " << leaf << "\n";
  return COM;
}

myvec3 Tree::computeAcc(const Particle &particle, myfloat theta)
    const // compute the accelartion applied on a particle by this node
{
  const myfloat softening_param = 0.025;
  myvec3 acc{0, 0, 0};

  if (leaf) {
    for (const auto &sp : particles) {
      if (sp->p != particle.p) // if the queried particle isn't the input
                               // particle (found by comparing positions)
      {
        auto diff =
            sp->p - particle.p; // Calculate the difference between the input
                                // particle and the particle in this node
        // sp->p dereferences the shared pointer to get the particle object
        auto distance = glm::length(diff);
        acc += diff * sp->m /
               (distance * distance * distance +
                softening_param); // softened gravitational equation
      }
    }
  } else {
    if (less_than_theta(particle.p, theta)) { // Barnes-Hut threshold
      // if the threshold is met, approximate the acceleration using the center
      // of mass instead of summing the individual particle contributions
      auto diff = COM.p - particle.p;
      auto distance = glm::length(diff);
      acc = diff * COM.m / (distance * distance * distance + softening_param);

      // std::cout << "Using center of mass: " << COM << "at node lvl: " <<
      // level<<"\n";

    } else { // if threshold not met, compute the acceleration due to the
      // branches inside this node
      for (const Tree &b : branches) {
        acc += b.computeAcc(particle, theta);
      }
    }
  }

  return acc;
}

bool Tree::less_than_theta(const myvec3 pos, double theta) const {
  // returns whether particle meets the threshold criteria for approximating the
  // acceleration due to this node
  myfloat distance = glm::length2(pos - COM.p);
  // TODO: this should not be the length of the cuboid
  return dimension < theta * distance;
}

std::string Tree::print() const {
  std::ostringstream str;
  str << cuboid.print() << " level: " << level
      << " nPoints:" << particles.size() << " Center of Mass: " << COM
      << ", leaf:" << leaf << "\n";
  for (const Tree &b : branches) {
    str << b.print();
  }
  return str.str();
}
