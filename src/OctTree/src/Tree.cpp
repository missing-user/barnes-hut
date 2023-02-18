#include "Tree.h"

// #define USE_CENTER_OF_MASS_FOR_SPLITTING

int Tree::maxDepth = 32;
int Tree::maxParticles = 4;

Tree::Tree(const Cuboid &cuboidIn, int levelIn)
    : cuboid(cuboidIn), level(levelIn), COM({{0, 0, 0}, 0}), leaf(true), divisor(cuboid.center) {}

Tree::Tree(const std::vector<Particle> &particles_in)
    : cuboid(bounding_box(particles_in)), level(0), COM({{0, 0, 0}, 0}),
      leaf(true), divisor(cuboid.center)
{
  for (const auto &p : particles_in)
  {
    particles.push_back(std::make_unique<Particle>(p));
  }
  subdivide();
}

void Tree::insert(std::unique_ptr<Particle> p){
  if (leaf){
    particles.push_back(std::move(p));
  }
  else{
    branches[selectOctant(p->p)].insert(std::move(p));
  }
}

void Tree::createBranches(const myvec3 &pos)
{ // populates the branches array of this object
  // with new trees from subdividing this node
  branches.reserve(
      8); // reserve space for 8 branches to save on resize operations
  for (Cuboid &r : cuboid.subdivideAtP(pos))
  {
    Tree branch = Tree(r, level + 1);
    branches.push_back(std::move(branch));
  }
}

int Tree::selectOctant(const myvec3 &pos) const
{
  int octant = 0;
  octant += (pos.x > divisor.x) << 0;
  octant += (pos.y > divisor.y) << 1;
  octant += (pos.z > divisor.z) << 2;
  return octant;
}

void Tree::subdivide()
{
  if (level < maxDepth && (particles.size() > maxParticles && leaf))
  {
    computeCOM(); // Compute COM while 
    createBranches(divisor); // If there aren't any branches, create them. 
    for (auto &pnt : particles) // This is parallelizable
    {
      branches[selectOctant(pnt->p)].insert(std::move(pnt));
    }
    particles.clear();
    for (auto &b : branches)
    {
      b.subdivide();
    }
    leaf = false;
  }
}

CenterOfMass Tree::computeCOM()
{ // calculate the center of mass for this node
  // and save it as the new COM
  COM.m = 0;
  COM.p *= 0;

  for (const auto &p : particles)
  {
    COM += *p; // Adding two particles creates a new particle which
                // represents the center of mass of the two particles
  }

  // If we are a leaf, then branches is empty and no loop will not execute
  for (Tree &t : branches)
  {
    COM += t.computeCOM();
  }

  #ifdef USE_CENTER_OF_MASS_FOR_SPLITTING
    divisor = COM.p;
  #endif

  return COM;
}

myvec3 Tree::computeAcc(const Particle &p1,
                        myfloat theta)
    const
{ // compute the accelartion applied on a particle by this node
  const auto pos = p1.p;
  return computeAccFromPos(pos, theta);
}

myvec3 Tree::computeAccFromPos(
    const myvec3 &pos,
    myfloat theta) const
{ // compute the total acceleration at this position
  // due to all particles in this Tree
  myvec3 acc{0, 0, 0};

  if (leaf)
  {
    for (const auto &sp : particles)
    {
      if (sp->p == pos)
        continue;

      acc += accelFunc(sp->p - pos, sp->m);
    }
  }
  else
  {
    if (lessThanTheta(pos, theta))
    { // Barnes-Hut threshold
      // if the threshold is met, approximate the acceleration using the center
      // of mass instead of summing the individual particle contributions
      acc = accelFunc(COM.p - pos, COM.m);
    }
    else
    { // if threshold not met, compute the acceleration due to the
      // branches inside this node
      for (const Tree &b : branches)
      {
        acc += b.computeAccFromPos(pos, theta);
      }
    }
  }

  return acc;
}

bool Tree::lessThanTheta(const myvec3 &pos, double theta) const
{
  return cuboid.diagonal2 < theta * glm::length2(pos - COM.p);
}

/***********************************************************************/
/* Helpers and Statistics */
/***********************************************************************/

void Tree::print() const
{
  std::cout << " level: " << level << " nPoints:" << particles.size()
            << " Center of Mass: " << COM << ", leaf:" << leaf << "\n";
  for (const auto &b : branches)
  {
    b.print();
  }
}

std::vector<std::size_t> Tree::DFS() const
{ // Depth first search, returns a vector with the order of the particles
  std::vector<std::size_t> indices;
  if (leaf)
  {
    indices.reserve(particles.size());
    for (const auto &p : particles)
    {
      indices.push_back(p->id);
    }
  }
  else
  {
    for (const auto &b : branches)
    {
      const auto &v = b.DFS();
      indices.insert(indices.end(), v.begin(), v.end());
    }
  }
  return indices;
}

std::pair<int, int> Tree::MaxDepthAndParticles() const
{ // returns the maximum depth and number of particles in the Tree recursively (depth, particles)
  if (leaf)
  {
    return {level, particles.size()};
  }
  else
  {
    std::pair<int, int> max = {0, 0};
    for (const auto &b : branches)
    {
      const auto &v = b.MaxDepthAndParticles();
      max.first = std::max(max.first, v.first);
      max.second = std::max(max.second, v.second);
    }
    return max;
  }
}

std::vector<DrawableCuboid> Tree::GetBoundingBoxes() const
{ // returns a vector of drawable cuboids for the bounding boxes of the Tree
  std::vector<DrawableCuboid> boxes;
  if (leaf)
  {
    boxes.push_back(std::move(DrawableCuboid(cuboid, level)));
  }
  else
  {
    for (const auto &b : branches)
    {
      const auto &v = b.GetBoundingBoxes();
      boxes.insert(boxes.end(),
                   std::make_move_iterator(v.begin()),
                   std::make_move_iterator(v.end()));
    }
  }
  return boxes;
}
