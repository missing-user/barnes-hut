#include "Tree.h"

// #define USE_CENTER_OF_MASS_FOR_SPLITTING

int Tree::maxDepth = 64;
int Tree::maxParticles = 1;
const Particle* Tree::firstP = nullptr;

Tree::Tree()
    : cuboid(Cuboid({0,0,0},{0,0,0})), level(-1), 
    COM({{0, 0, 0}, 0}), branches(nullptr) {}

Tree::Tree(const Cuboid &cuboidIn, int levelIn)
: cuboid(cuboidIn), level(levelIn), COM({{0, 0, 0}, 0}), branches(nullptr) {}


Tree::Tree(const std::vector<Particle> &particles_in)
    : cuboid(bounding_box(particles_in)), level(0), COM({{0, 0, 0}, 0}),
      branches(nullptr)
{
  firstP = const_cast<Particle *>(&particles_in[0]);
  for (auto &p : particles_in)
  {
    particles.push_back(std::distance(&particles_in[0], &p));
  }
  subdivide();
}

Tree::~Tree(){
  delete [] branches;
}

void Tree::insert(int p){
  if (branches){
    branches[selectOctant((firstP+p)->p)].insert(p);
  }else{
    particles.push_back(p);
  }
}

void Tree::createBranches(const myvec3 &pos)
{ // populates the branches array of this object
  // with new trees from subdividing this node
  int i = 0;
  branches = new Tree[8];
  for (Cuboid &r : cuboid.subdivideAtP(pos))
  {
    branches[i] = Tree(r, level + 1);
    i++;
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

void Tree::subdivide()
{
  if (level < maxDepth && (particles.size() > maxParticles && !branches))
  {
    computeCOM(); // Compute COM while 
    createBranches(cuboid.center); // If there aren't any branches, create them. 
    for (auto &pnt : particles) // This is parallelizable
    {
      branches[selectOctant((firstP+pnt)->p)].insert(pnt);
    }
    particles.clear();
    for (int i = 0; i < 8; i++)
    {
      branches[i].subdivide();
    }
  }
}

CenterOfMass Tree::computeCOM()
{ // calculate the center of mass for this node
  // and save it as the new COM
  COM.m = 0;
  COM.p *= 0;

  for (const auto &p : particles)
  {
    COM += *(firstP+p); // Adding two particles creates a new particle which
                // represents the center of mass of the two particles
  }

  // If we are a leaf, then branches is empty and no loop will not execute
  if(branches)
    for (int i = 0; i < 8; i++)
    {
      COM += branches[i].computeCOM();
    }

  return COM;
}

myvec3 Tree::computeAcc(const Particle &p1,
                        myfloat theta2)
    const
{ // compute the accelartion applied on a particle by this node
  const auto pos = p1.p;
  return computeAccFromPos(pos, theta2);
}

myvec3 Tree::computeAccFromPos(
    const myvec3 &pos,
    myfloat theta2) const
{ // compute the total acceleration at this position
  // due to all particles in this Tree
  myvec3 acc{0, 0, 0};

  if (!branches)
  {
    for (const auto &sp : particles)
    {
      if ((firstP+sp)->p == pos)
        continue;

      acc += accelFunc((firstP+sp)->p - pos, (firstP+sp)->m);
    }
  }
  else
  {
    if (lessThanTheta(pos, theta2))
    { // Barnes-Hut threshold
      // if the threshold is met, approximate the acceleration using the center
      // of mass instead of summing the individual particle contributions
      acc = accelFunc(COM.p - pos, COM.m);
    }
    else
    { // if threshold not met, compute the acceleration due to the
      // branches inside this node
      for (int i = 0; i < 8; i++)
      {
        acc += branches[i].computeAccFromPos(pos, theta2);
      }
    }
  }

  return acc;
}

// Multipole acceptance criteria, uses theta^2 to avoid a sqrt() operation
bool Tree::lessThanTheta(const myvec3 &pos, double theta2) const
{
  return cuboid.diagonal2 < theta2 * glm::length2(pos - COM.p);
}

/***********************************************************************/
/* Helpers and Statistics */
/***********************************************************************/

void Tree::print() const
{
  std::cout << " level: " << level << " nPoints:" << particles.size()
            << " Center of Mass: " << COM << ", leaf:" << (!branches) << "\n";
  if(branches)
    for (int i = 0; i < 8; i++)
    {
      branches[i].print();
    }
}

std::vector<std::size_t> Tree::DFS() const
{ // Depth first search, returns a vector with the order of the particles
  std::vector<std::size_t> indices;
  if (branches)
  {
      for (int i = 0; i < 8; i++)
      {
        const auto &v = branches[i].DFS();
        indices.insert(indices.end(), v.begin(), v.end());
      }
  }else{
    indices.reserve(particles.size());
    for (auto &p : particles)
    {
      indices.push_back(p);
    }
  }
  return indices;
}

std::pair<int, int> Tree::MaxDepthAndParticles() const
{ // returns the maximum depth and number of particles in the Tree recursively (depth, particles)
  if (branches)
  {
    std::pair<int, int> max = {0, 0};
    for (int i = 0; i < 8; i++)
    {    
      const auto &v = branches[i].MaxDepthAndParticles();
      max.first = std::max(max.first, v.first);
      max.second = std::max(max.second, v.second);
    }
    return max;
  }else{
    return {level, particles.size()};
  }
}

std::vector<DrawableCuboid> Tree::GetBoundingBoxes() const
{ // returns a vector of drawable cuboids for the bounding boxes of the Tree
  std::vector<DrawableCuboid> boxes;
  if (branches)
  {
    for (int i = 0; i < 8; i++)
    {    
      const auto &v = branches[i].GetBoundingBoxes();
      boxes.insert(boxes.end(),
                   std::make_move_iterator(v.begin()),
                   std::make_move_iterator(v.end()));
    }
  }else{
    boxes.push_back(std::move(DrawableCuboid(cuboid, level)));
  }
  return boxes;
}
