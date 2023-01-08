#include "Tree.h"

int Tree::maxDepth = 64;
int Tree::maxParticles = 1;

Tree::Tree(const Cuboid &cuboidIn, int levelIn)
    : cuboid(cuboidIn), level(levelIn), COM({{0, 0, 0}, 0}), leaf(true) {}

Tree::Tree(const std::vector<Particle> &particles)
    : cuboid(bounding_box(particles)), level(0), COM({{0, 0, 0}, 0}),
      leaf(true) {

  for (const auto &p : particles) {
    insert(p);
  }
  computeCOM();
}

void Tree::createBranches() { // populates the branches array of this object
                              // with new trees from subdividing this node
  branches.reserve(
      8); // reserve space for 8 branches to save on resize operations
  for (Cuboid &r : cuboid.subdivide()) {
    Tree branch = Tree(r, level + 1);
    branches.push_back(std::move(branch));
  }
}

int Tree::selectOctant(const myvec3 &pos) const {
  int octant = 0;
  octant += (pos.x > cuboid.center.x) << 0;
  octant += (pos.y > cuboid.center.y) << 1;
  octant += (pos.z > cuboid.center.z) << 2;
  return octant;
}

void Tree::insert(const Particle &p) { // adds a point to the tree structure.
  // Depending how full the node is, new branches may be generated
  insert(std::make_unique<Particle>(p));
}

void Tree::insert(std::unique_ptr<Particle> p) {
  if (level < maxDepth && (particles.size() >= maxParticles || !leaf)) {
    if (leaf) {
      leaf = false;
      createBranches(); // If there aren't any branches, create them.

      for (auto &pnt : particles) {
        branches[selectOctant(pnt->p)].insert(std::move(pnt));
      }
      particles.clear();
    }
    // Also add the new point to the appropriate branch
    branches[selectOctant(p->p)].insert(std::move(p));
  } else {
    particles.push_back(std::move(p));
  }
}

CenterOfMass Tree::computeCOM() { // calculate the center of mass for this node
                                  // and save it as the new COM
  // Assumtions:
  // This function is only called once per tree, and COM is initialized to zero

  if (leaf) {
    // if this node doesn't have branches, calculate the center of mass the
    // contained particles
    for (const auto &p : particles) {
      COM += *p; // Adding two particles creates a new particle which
                 // represents the center of mass of the two particles
    }
  } else {
    // This recursively calculates the center of mass of each branch sets it
    for (Tree &t : branches) {
      COM += t.computeCOM();
    }
  }

  return COM;
}

myvec3 Tree::computeAcc(const Particle &p1,
                        myfloat theta)
    const { // compute the accelartion applied on a particle by this node
  const auto pos = p1.p;
  return computeAccFromPos(pos, theta);
}

myvec3 Tree::computeAccFromPos(
    const myvec3 &pos,
    myfloat theta) const { // compute the total acceleration at this position
                           // due to all particles in this tree
  myvec3 acc{0, 0, 0};

  if (leaf) {
    for (const auto &sp : particles) {
      if (sp->p == pos)
        continue;

      acc += accelFunc(sp->p - pos, sp->m);
    }
  } else {
    if (less_than_theta(pos, theta)) { // Barnes-Hut threshold
      // if the threshold is met, approximate the acceleration using the center
      // of mass instead of summing the individual particle contributions
      acc = accelFunc(COM.p - pos, COM.m);

    } else { // if threshold not met, compute the acceleration due to the
      // branches inside this node
      for (const Tree &b : branches) {
        acc += b.computeAccFromPos(pos, theta);
      }
    }
  }

  return acc;
}

bool Tree::less_than_theta(const myvec3 &pos, double theta) const {
  return cuboid.diagonal2 < theta * glm::length2(pos - COM.p);
}

void Tree::print() const {
  std::cout << " level: " << level << " nPoints:" << particles.size()
            << " Center of Mass: " << COM << ", leaf:" << leaf << "\n";
  for (const auto &b : branches) {
    b.print();
  }
}

std::vector<std::size_t> Tree::DFS() const {
  std::vector<std::size_t> indices;
  if (leaf) {
    indices.reserve(particles.size());
    for (const auto &p : particles) {
      indices.push_back(p->id);
    }
  } else {
    for (const auto &b : branches) {
      const auto &v = b.DFS();
      indices.insert(indices.end(), v.begin(), v.end());
    }
  }
  return indices;
}

std::pair<int, int> Tree::MaxDepthAndParticles() const {
  if (leaf) {
    return {level, particles.size()};
  } else {
    std::pair<int, int> max = {0, 0};
    for (const auto &b : branches) {
      const auto &v = b.MaxDepthAndParticles();
      max.first = std::max(max.first, v.first);
      max.second = std::max(max.second, v.second);
    }
    return max;
  }
}