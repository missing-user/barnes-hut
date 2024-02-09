#include "Barneshut.h"
#include "Order.h"
#include <queue>
#include <libmorton/morton.h>
#include <chrono>

//#define DEBUG_BUILD
#ifdef DEBUG_BUILD
#define DEBUG(x) do { std::cout << x; } while (0)
#else
#define DEBUG(x)
#endif

// Due to morton order we know that all current particles.pos are > node_pos, only check the upper bounds
#define IN_BOUNDS(i) ((particles.p.x[i] <= node_pos.x + node_dim.x) &&\
                      (particles.p.y[i] <= node_pos.y + node_dim.y) &&\
                      (particles.p.z[i] <= node_pos.z + node_dim.z))

struct Node{
  // Check Figure 3. in https://www.tabellion.org/et/paper11/OutOfCorePBGI.pdf 
  int start; // index of the first particle
             // If internal node, start is the index of the first child. All Nodes up to start+8 are children 
  int count; // number of particles, since the particles are sorted, the last particle is start+count-1
             // Negative count indicates that this is an internal node
  myvec3 center_of_mass{0,0,0};
  myfloat mass{0};

  bool isLeaf() const { return count >= 0; }
  void setInternal(int first_child){
    start = first_child;
    count = -1;
  }
};

void bh_superstep(Particles& particles, size_t count){

  auto time1 = std::chrono::high_resolution_clock::now();
  auto boundingbox = bounding_box(particles.p, count);
  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Bounding Box calculation took " << elapsed.count()<<"s "<< std::endl;
  computeAndOrder(particles, boundingbox);
  time1 = std::chrono::high_resolution_clock::now();


  DEBUG("Particles: \n");
  for (size_t i = 0; i < count; i++)
  {
    DEBUG("("<<particles.p.x[i] << "," << particles.p.y[i] << "," << particles.p.z[i] <<"),"<<"\n");
  }
  
  // Create the tree
  const int depth_max = 15;
  const int leaf_max = 4; // maximum particles per node
  std::vector<DrawableCuboid> drawcuboids;
  std::array<std::vector<Node>, depth_max+1> tree;
  // build_tree();
  /* Since we know the particles are sorted in Morton order, we can use out of core
  * Tree construction as described in http://www.thesalmons.org/john/pubs/siam97/salmon.pdf
  * by only keeping track of the group currently being constructed and finalizing it 
  * once the first particle outside the group is found. (morton order ensures there will be 
  * no subsequent particles in the group)
  */
  size_t i = 0;
  int depth = 0; // 0 = root
  myvec3 node_pos{boundingbox.min()}; // lower left corner
  myvec3 node_dim{boundingbox.dimension}; // length of each dimension
  short stack[depth_max+1];
  std::memset(stack, 0, sizeof(stack));

  std::queue<size_t> fifo; // indices of the particles
  fifo.push(i);
  while(i<particles.size()){
    // Add the first leaf_max+1 particles to the queue
    node_pos = boundingbox.min();
    for (int d = 1; d <= depth; d++)
    {
      node_pos.x += (stack[d] & 1) * boundingbox.dimension.x/(1<<d);
      node_pos.y += ((stack[d]>>1) & 1) * boundingbox.dimension.y/(1<<d);
      node_pos.z += ((stack[d]>>2) & 1) * boundingbox.dimension.z/(1<<d);
    }
    for (int j = fifo.size(); j <= leaf_max; j++)
    {
      if(i >= particles.size()-1) [[unlikely]] { break; } // Early exit before we exceed the number of particles
      fifo.push(++i);
      DEBUG(i<<" pushed into fifo: " << fifo.size() <<"\n");
    }
    // Find the first depth level that does not contain all particles
    // i.e. the last particle (morton order)
    node_dim = boundingbox.dimension/static_cast<myfloat>(1<<depth);
    while (IN_BOUNDS(i) && depth < depth_max && fifo.size() > 1)
    { 
      depth++;
      node_dim = boundingbox.dimension/static_cast<myfloat>(1<<depth);
      DEBUG("subdivide, i="<<i<<" is still in range"<<node_dim<<(fifo.size())<<"\n");
    }
    if(depth >= depth_max) { 
      // Max depth reached, ignore the particle count limit and fill the node 
      // To avoid infinite deepening of the tree when particles are on the same position
      DEBUG("Max depth reached, ignoring particle count limit\n");
      while(IN_BOUNDS(i) && i < particles.size()-1){
        fifo.push(++i);
        DEBUG(i<<"[Override] pushed into fifo: " << fifo.size() <<"\n");
      }  
    }

    DEBUG("pos "<<node_pos<<" node_dim "<<node_dim <<" and stack: ");
    for (int d = 0; d <= depth; d++)
    {
      DEBUG((int)stack[d]<<" ");
    }
    DEBUG(std::endl);

    Node leaf;
    leaf.start = fifo.front(); // index of the first particle
    leaf.count = 0; // number of particles, since the particles are sorted, all particles up to start+count-1 are contained
    // Add all particles that are in bounds to the node
    while(IN_BOUNDS(fifo.front()) && fifo.size() > 0){
      DEBUG("popped "<<fifo.front()<<" ("<<particles.p.x[fifo.front()] << " " << particles.p.y[fifo.front()] << " " << particles.p.z[fifo.front()]<<")" << std::endl);
      fifo.pop();
      leaf.count++;
    }
    tree[depth].push_back(leaf);
    DEBUG("Finalizing Leaf at depth "<<depth<<" with num_planets="<<leaf.count<<std::endl);
    drawcuboids.push_back(DrawableCuboid(minMaxCuboid(node_pos, node_pos+node_dim), depth));
    
    // Go to next node (at this depth if possible)
    if(stack[depth] < 7){
      stack[depth]++;
    }else{
      // Go to the next depth
      while(stack[depth] == 7){
        stack[depth] = 0;

        Node node;
        node.setInternal(tree[depth].size()-8); // index of the first child (8 children per node)

        depth--;
        tree[depth].push_back(node);
        DEBUG("Finalizing Node at depth "<<depth<<std::endl);
        node_dim *= 2;

      }
      stack[depth]++;
      if(depth == 0){
        // We have finished the tree
        break;
      }
    }
  }

  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Tree building " << elapsed.count()<<"s"<< std::endl;

  time1 = std::chrono::high_resolution_clock::now();
  
  // COM and mass calculation
  for (int d = depth_max-1; d >= 0; d--) // unsigned int underflows to max value
  {
    DEBUG("Depth "<<d<<" has "<<tree[d].size()<<" nodes\n");
    for (size_t i = 0; i < tree[d].size(); i++)
    {
      auto& currentnode = tree[d][i];
      if(currentnode.isLeaf()){
        DEBUG("Leaf "<<i<<" at depth "<<d<<" has "<<currentnode.count<<" particles\n");
        for (size_t j = currentnode.start; j < currentnode.start+currentnode.count; j++)
        {
          currentnode.center_of_mass.x += particles.p.x[j] * particles.m[j];
          currentnode.center_of_mass.y += particles.p.y[j] * particles.m[j];
          currentnode.center_of_mass.z += particles.p.z[j] * particles.m[j];
          currentnode.mass += particles.m[j];
        }
        if(currentnode.mass != 0)
          currentnode.center_of_mass /= currentnode.mass;
      }else{
        // TODO: SIMD this, test unrolling
        DEBUG("Internal Node "<<i<<" with start "<<currentnode.start);
        for (size_t j = currentnode.start; j < currentnode.start + 8; j++)
        {
          currentnode.center_of_mass += tree[d+1][j].center_of_mass * tree[d+1][j].mass;
          currentnode.mass += tree[d+1][j].mass;

        }
        DEBUG(<<"=d COM "<<currentnode.center_of_mass<<" and mass "<<currentnode.mass<<"\n");
        if(currentnode.mass != 0)
          currentnode.center_of_mass /= currentnode.mass;
      }
    }
  }
  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "COM computation "<<tree[0][0].center_of_mass<<"took " << elapsed.count()<<"s"<< std::endl;
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