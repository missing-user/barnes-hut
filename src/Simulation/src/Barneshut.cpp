#include "Barneshut.h"
#include "Order.h"
#include <queue>
#include <libmorton/morton.h>
#include <chrono>

#ifdef DEBUG_BUILD
#define DEBUG(x) do { std::cout << x; } while (0)
#else
#define DEBUG(x)
#endif

// Due to morton order we know that all current particles.pos are > node_pos, only check the upper bounds
#define IN_BOUNDS(i) ((particles.p.x[i] <= node_pos.x + node_dim.x) &&\
                      (particles.p.y[i] <= node_pos.y + node_dim.y) &&\
                      (particles.p.z[i] <= node_pos.z + node_dim.z))

// bool IN_BOUNDS(int i, const Vectors& p, const myvec3& node_pos, const myvec3& node_dim){
//   DEBUG("  IN_BOUNDS x: "<<p.x[i] << "<"<<(node_pos.x+node_dim.x)<<" "<<(p.x[i] <= node_pos.x + node_dim.x)<< "\n")
//   DEBUG("  IN_BOUNDS y: "<<p.y[i] << "<"<<(node_pos.y+node_dim.y)<<" "<<(p.y[i] <= node_pos.y + node_dim.y)<< "\n")
//   DEBUG("  IN_BOUNDS z: "<<p.z[i] << "<"<<(node_pos.z+node_dim.z)<<" "<<(p.z[i] <= node_pos.z + node_dim.z)<< "\n")
//   return ((p.x[i] <= node_pos.x + node_dim.x) && (p.y[i] <= node_pos.y + node_dim.y) && (p.z[i] <= node_pos.z + node_dim.z));
// }

struct Node{
  int start; // index of the first particle
  int count; // number of particles, since the particles are sorted, the last particle is start+count-1
  myvec3 center_of_mass;
  myfloat mass;
  // The immediate child of this node is always at index-1 if it exists.
  // Check if the child exists by checking if the index is negative
  // Check Figure 3. in https://www.tabellion.org/et/paper11/OutOfCorePBGI.pdf 
  int prev_sibling; // index of the previous sibling is unknown, since children may contain subnodes
};

void bh_superstep(Particles& particles, size_t count){
  auto boundingbox = bounding_box(particles.p, count);
  computeAndOrder(particles, boundingbox);

  auto time1 = std::chrono::high_resolution_clock::now();

  DEBUG("Particles: \n");
  for (size_t i = 0; i < count; i++)
  {
    DEBUG("("<<particles.p.x[i] << "," << particles.p.y[i] << "," << particles.p.z[i] <<"),"<<"\n");
  }
  
  // Create the tree
  std::vector<DrawableCuboid> drawcuboids;
  std::vector<Node> nodes;
  // build_tree();
  /* Since we know the particles are sorted in Morton order, we can use out of core
  * Tree construction as described in http://www.thesalmons.org/john/pubs/siam97/salmon.pdf
  * by only keeping track of the group currently being constructed and finalizing it 
  * once the first particle outside the group is found. (morton order ensures there will be 
  * no subsequent particles in the group)
  */
  const int depth_max = 21;
  const int leaf_max = 4; // maximum particles per node
  int i = 0;
  int depth = 0; // 0 = root
  myvec3 node_pos{boundingbox.min()}; // lower left corner
  myvec3 node_dim{boundingbox.dimension}; // length of each dimension
  char stack[depth_max+1];
  memset(stack, 0, sizeof(stack));

  std::queue<int> fifo; // indices of the particles
  while(i<particles.size()){
    // Add the first leaf_max+1 particles to the queue
    node_pos = boundingbox.min();
    for (int d = 1; d <= depth; d++)
    {
      node_pos.x += (stack[d] & 1) * boundingbox.dimension.x/(1<<d);
      node_pos.y += ((stack[d]>>1) & 1) * boundingbox.dimension.y/(1<<d);
      node_pos.z += ((stack[d]>>2) & 1) * boundingbox.dimension.z/(1<<d);
    }

    while (fifo.size() < leaf_max+1 && i<particles.size()-1)
    {
      fifo.push(i);
      i++;
    }
    DEBUG(i<<" pushed fifo: " << fifo.size() << std::endl);
    // Find the first depth level that does not contain all particles
    // i.e. the last particle (morton order)
    node_dim = boundingbox.dimension/static_cast<myfloat>(1<<depth);
    while (IN_BOUNDS(i) && depth < depth_max && fifo.size() > 0)
    { 
      depth++;
      node_dim = boundingbox.dimension/static_cast<myfloat>(1<<depth);
    }

    DEBUG("pos "<<node_pos<<" node_dim "<<node_dim <<" and stack: ");
    for (int d = 0; d <= depth; d++)
    {
      DEBUG((int)stack[d]<<" ");
    }
    DEBUG(std::endl);

    Node node;
    node.start = fifo.front();
    node.count = 0;
    // Add all particles that are in bounds to the node
    while(IN_BOUNDS(fifo.front()) && fifo.size() > 0){
      //if(fifo.size() == 0)
      //  throw std::runtime_error("Fifo is empty, but we are still in bounds. This should not happen, since the last particle must be out of bounds. Exceeded max depth?");
      
      DEBUG("popped "<<particles.p.x[fifo.front()] << " " << particles.p.y[fifo.front()] << " " << particles.p.z[fifo.front()] << std::endl);
      fifo.pop();
      node.count++;
    }
    nodes.push_back(node);
    DEBUG("Finalizing Leaf at depth "<<depth<<" with planets "<<node.count<<std::endl);
    drawcuboids.push_back(DrawableCuboid(minMaxCuboid(node_pos, node_pos+node_dim), depth));
    
    // Go to next node (at this depth if possible)
    if(stack[depth] < 7){
      stack[depth]++;
    }else{
      // Go to the next depth
      while(stack[depth] == 7){
        stack[depth] = 0;
        depth--;
        node_dim *= 2;
        DEBUG("Finalizing Node at depth "<<depth<<std::endl);
      }
      stack[depth]++;
      if(depth == 0){
        // We have finished the tree
        break;
      }
    }
  }

  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Tree building " << elapsed.count()<<"s"<< std::endl;


  // COM and mass calculation
  for (size_t i = 0; i < nodes.size(); i++)
  {
    myvec3 com{0,0,0};
    myfloat mass = 0;
    for (size_t j = nodes[i].start; j < nodes[i].start+nodes[i].count; j++)
    {
      com.x += particles.p.x[j] * particles.m[j];
      com.y += particles.p.y[j] * particles.m[j];
      com.z += particles.p.z[j] * particles.m[j];
      mass += particles.m[j];
    }
    com /= mass;
    nodes[i].center_of_mass = com;
    nodes[i].mass = mass;
  }
  
}

void stepSimulation(Particles& particles, myfloat dt, double theta) {
  // barnes hut optimized step
  // Buffer for the new state vector of particles
  Vectors acc{particles.size()};
  Vectors p2{particles.size()};

  bh_superstep(particles, particles.size());
  // We have constructed the tree, use it to efficiently compute the
  // accelerations. Far away particles get grouped and their contribution is
  // approximated using their center of mass (Barnes Hut algorithm)

  // Use velocity verlet algorithm to update the particle positions and
  // velocities
  // https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < particles.size(); i++) {
    // First update the positions
    //mytree.computeAccFromPos(acc[i], particles.p[i], theta);
    p2.x[i] = particles.p.x[i] + particles.v.x[i] * dt + acc.x[i] * dt * dt / 2.;
    p2.y[i] = particles.p.y[i] + particles.v.y[i] * dt + acc.y[i] * dt * dt / 2.;
    p2.z[i] = particles.p.z[i] + particles.v.z[i] * dt + acc.z[i] * dt * dt / 2.;
  }

  // New Octree using the updated particle positions
  // Since the following loop will not change the positions of particles_next,
  // we can safely reuse the vector for iteration and acceleration calculation
  // TODO: Reuse this octree for the next timestep
  
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < particles.size(); i++) {
    // Then update the velocities using v(t+1) = dt*(a(t) + a(t+dt))/2
    myfloat currentaccx, currentaccy, currentaccz;
    GlmView currentacc{&currentaccx, &currentaccy, &currentaccz};
    //mytree.computeAccFromPos(currentacc, p2[i], theta);
    particles.v.x[i] += (*currentacc.x + acc.x[i]) * dt / 2.;
    particles.v.y[i] += (*currentacc.y + acc.y[i]) * dt / 2.;
    particles.v.z[i] += (*currentacc.z + acc.z[i]) * dt / 2.;
  }
}