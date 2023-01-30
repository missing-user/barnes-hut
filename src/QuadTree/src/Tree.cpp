#include "Tree.h"

int Tree::maxDepth = 64;
int Tree::maxParticles = 1;

Tree::Tree(const Cuboid &cuboidIn, int levelIn)
    : cuboid(cuboidIn), level(levelIn), COM({{0, 0, 0}, 0}), leaf(true) {}

Tree::Tree(const std::vector<Particle> &particles_in)
    : cuboid(bounding_box(particles_in)), level(0), COM({{0, 0, 0}, 0}),
      leaf(true)
{
  for (const auto &p : particles_in)
  {
    insertNonRecursive(p);
  }
  subdivideBreadth();
  computeCOM();
}

void Tree::createBranches()
{ // populates the branches array of this object
  // with new trees from subdividing this node
  branches.reserve(
      8); // reserve space for 8 branches to save on resize operations
  for (Cuboid &r : cuboid.subdivide())
  {
    Tree branch = Tree(r, level + 1);
    branches.push_back(std::move(branch));
  }
}

void Tree::createBranchesAtP(myvec3 P)
{ // populates the branches array of this object
  // with new trees from subdividing this node
  branches.reserve(
      8); // reserve space for 8 branches to save on resize operations
  for (Cuboid &r : cuboid.subdivideAtP(P))
  {
    Tree branch = Tree(r, level + 1);
    branches.push_back(std::move(branch));
  }
}

int Tree::selectOctant(const myvec3 &pos) const
{
  int octant = 0;
  octant += (pos.x > cuboid.center.x) << 0;
  octant += (pos.y > cuboid.center.y) << 1;
  octant += (pos.z > cuboid.center.z) << 2;
  return octant;
}

void Tree::insertNonRecursive(const Particle &p)
{
  insertNonRecursive(std::make_unique<Particle>(p));
}
void Tree::insertNonRecursive(std::unique_ptr<Particle> p)
{
  particles.push_back(std::move(p));
}

void Tree::subdivideBreadth()
{
  if (level < maxDepth && (particles.size() > maxParticles && leaf))
  {
    leaf = false;
    computeCOM_NonRecursive();
    createBranchesAtP(COM.p);   // If there aren't any branches, create them.
    for (auto &pnt : particles) // This is parallelizable
    {
      branches[selectOctant(pnt->p)].insertNonRecursive(std::move(pnt));
    }
    particles.clear();
  }
  for (auto &b : branches)
  {
    b.subdivideBreadth();
  }
}

CenterOfMass Tree::computeCOM_NonRecursive()
{ // calculate the center of mass for this node
  // and save it as the new COM
  COM.m = 0;
  COM.p *= 0;
  for (const auto &p : particles)
  {
    COM += *p; // Adding two particles creates a new particle which
               // represents the center of mass of the two particles
  }
  return COM;
}

CenterOfMass Tree::computeCOM()
{ // calculate the center of mass for this node
  // and save it as the new COM
  // Assumtions:
  // This function is only called once per tree, and COM is initialized to zero
  COM.m = 0;
  COM.p *= 0;
  
  if (leaf)
  {
    // if this node doesn't have branches, calculate the center of mass the
    // contained particles
    for (const auto &p : particles) /// Shouldn't COM be reset first?
    {
      COM += *p; // Adding two particles creates a new particle which
                 // represents the center of mass of the two particles
    }
  }
  else
  {
    // This recursively calculates the center of mass of each branch sets it
    for (Tree &t : branches)
    {
      COM += t.computeCOM();
    }
  }

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
  // due to all particles in this tree
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
    if (less_than_theta(pos, theta))
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

bool Tree::less_than_theta(const myvec3 &pos, double theta) const
{
  return cuboid.diagonal2 < theta * glm::length2(pos - COM.p);
}

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
{
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
{
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
{
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
