#include "Barneshut.h"
#include "Order.h"
#include <queue>
#include <bitset>
#include <libmorton/morton.h>
#include <chrono>

//#define DEBUG_BUILD
#ifdef DEBUG_BUILD
#define DEBUG_D(x, d) do { \
  for (int DEBUG_DEPTH = 0; DEBUG_DEPTH < d; DEBUG_DEPTH++){std::cout << " ";}\
  std::cout << x; } while (0)
#else
#define DEBUG_D(x, depth)
#endif

#define DEBUG(x) DEBUG_D(x, 0)

const int depth_max = 12;
const int leaf_max = 1; // maximum particles per node

void PRINT_BITS_BY3(uint_fast64_t x, int separator){
  for (size_t i = 1; i <= depth_max; i++)
  {
    auto prefix = std::bitset<3>(x>>(63-3*i));
    DEBUG(prefix);
    if(i == separator)
      DEBUG(":");
    else
      DEBUG(" ");
  }
}

void DEBUG_BITS(uint_fast64_t end,int depth){
  DEBUG_D("",depth);
  PRINT_BITS_BY3(end, depth);
  auto extsuf = std::bitset<63-3*depth_max>(end);
  assert(extsuf == std::bitset<63-3*depth_max>(std::numeric_limits<uint_fast64_t>::max())); //unused bottom part must be 1
  assert((end>>63) == 0); // Top bit must always be 0
}

// FIXME: accelFunc is a duplicate of the function in Forces.h, but for some reason the linker cant find it
//#include "Forces.h"
//#pragma omp declare simd linear(accx, accy, accz)
void accelFunc(myfloat*  accx, myfloat*  accy, myfloat*  accz, 
  myfloat diffx, myfloat diffy, myfloat diffz, myfloat mass) {
  constexpr myfloat softening_param = 0.025;
  auto r2 = length2(diffx, diffy, diffz);
  auto r = std::sqrt(r2);
  
  *accx += diffx * mass /(r2 * r + softening_param);
  *accy += diffy * mass /(r2 * r + softening_param);
  *accz += diffz * mass /(r2 * r + softening_param);
}

// Due to morton order we know that all current particles.pos are > node_pos, only check the upper bounds
// Also contains the array in bounds check (prevent segfault), i+1 to avoid overflow on unsigned int
#define IN_BOUNDS(i) ((i + 1 < particles.size()) && (mortoncodes[i] <= current_cell_end))

struct Node{
  // Check Figure 3. in https://www.tabellion.org/et/paper11/OutOfCorePBGI.pdf 
  int start; // index of the first particle
             // If internal node, start is the index of the first child. All Nodes up to start+8 are children 
  int count; // number of particles, since the particles are sorted, the last particle is start+count-1
             // Negative count indicates that this is an internal node

  bool isLeaf() const { return count >= 0; }
  void setInternal(int first_child){
    start = first_child;
    count = -1;
  }
};


inline bool isApproximationValid(myfloat dx,myfloat dy,myfloat dz, double theta, myfloat diagonal2)
{
  return diagonal2 < theta * length2(dx,dy,dz);
}

void recursive_force(
  const Particles& particles,
  const std::array<std::vector<Node>, depth_max+1>& tree,
  const myfloat* x, const myfloat* y, const myfloat* z, 
  const std::array<std::vector<myfloat>, depth_max+1>& commass, 
  const std::array<std::vector<myfloat>, depth_max+1>& comx,
  const std::array<std::vector<myfloat>, depth_max+1>& comy,
    const std::array<std::vector<myfloat>, depth_max+1>& comz,
     myfloat* accx,  myfloat* accy,  myfloat* accz, 
    const std::array<myfloat, depth_max+1>& diagonal2,
    int depth, int start){
  DEBUG_D("recursive_force enter with start "<<start<<"\n", depth);
  auto& node = tree[depth][start];
  if(node.isLeaf()){
    DEBUG_D("recursive_force leaf with "<<node.count<<" particles\n", depth);
    // Compute the force
    for (int i = node.start; i < node.start + node.count; i++)
    {
      if(!(i<particles.size())){
        std::cerr<<"At depth "<<depth<<" and start "<<start<<" with node "<<node.start<<" and count "<<node.count<<"\n";
        std::cerr<<"i="<<i<<" is out of bounds, particles.size()="<<particles.size()<<std::endl;
      }
      assert(node.count>0);
      assert(i < particles.size());
      accelFunc(accx, accy, accz,
                particles.p.x[i] - *x,
                particles.p.y[i] - *y,
                particles.p.z[i] - *z, particles.m[i]); // TODO replace 50 with the mass
    }
  }else{

    myfloat dx = comx.at(depth).at(start) - *x;
    myfloat dy = comy[depth][start] - *y;
    myfloat dz = comz[depth][start] - *z;
    if(isApproximationValid(dx,dy,dz, 1.5, diagonal2[depth])){
      DEBUG_D("recursive_force approximated with COM "<<start<<" = "<<myvec3(
        comx[depth][start], comy[depth][start], comz[depth][start]
      )<<"\n", depth);
      // Compute the COM force
      accelFunc(accx, accy, accz,dx,dy,dz,commass[depth][start]);
    }else{
      for (auto it = node.start; it < node.start+8; it++)
      {
        recursive_force(particles, tree, x,y,z, commass, comx, comy, comz,
        accx, accy, accz, diagonal2, depth+1, it);
      }
    }
  }

  DEBUG_D("recursive_force exit depth "<<depth<<"\n", depth);
}

std::vector<DrawableCuboid> bh_superstep(Particles& particles, size_t count, Vectors& acc){
  auto time1 = std::chrono::high_resolution_clock::now();
  auto boundingbox = bounding_box(particles.p, count);
  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Bounding Box calculation took " << elapsed.count()*1e3<<"ms "<< std::endl;
  
  time1 = std::chrono::high_resolution_clock::now();
  computeAndOrder(particles, boundingbox);
  auto mortoncodes = computeMortonCodes(particles, boundingbox);
  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Index reordering took " << elapsed.count()*1e3<<"ms"<< "\n";
  time1 = std::chrono::high_resolution_clock::now();


  DEBUG("Particles: \n");
  for (size_t i = 0; i < count; i++)
  {
    DEBUG("("<<particles.p.x[i] << "," << particles.p.y[i] << "," << particles.p.z[i] <<"), "<<std::bitset<63>(mortoncodes[i])<<"\n");
  }
  
  // Create the tree
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
  const myvec3 node_dim{boundingbox.dimension/std::pow(2.,21.)}; // length of each dimension
  uint_fast64_t current_cell_end = std::numeric_limits<uint_fast64_t>::max() - (1ul << 63); // First morton code that is not in the current cell
  short stack[depth_max+1];
  std::memset(stack, 0, sizeof(stack));

  std::queue<size_t> fifo; // indices of the particles
  fifo.push(i);
  std::vector<DrawableCuboid>  drawcuboids{};

  while(i<particles.size()){
    // Add the first leaf_max+1 particles to the queue
    for (int j = fifo.size(); j <= leaf_max; j++)
    {
      if(i+1 >= particles.size()) [[unlikely]] { break; } // Early exit before we exceed the number of particles
      fifo.push(++i);
      //DEBUG_D(i<<" pushed into fifo: " << fifo.size() <<"\n", depth);
    }
    // Find the first depth level that does not contain all particles
    // i.e. the last particle (morton order)
    while (IN_BOUNDS(i) && depth < depth_max)
    { 
      depth++;
      stack[depth] = 0;
      current_cell_end -= 7ul << (63-3*depth);
      DEBUG_BITS(current_cell_end, depth);
      DEBUG_D(" subdivide, i="<<i<<" is still in range"<<"\n", depth);

      if(depth >= depth_max){
        // Max depth reached, ignore the particle count limit and fill the node 
        // To avoid infinite deepening of the tree when particles are on the same position
        DEBUG("Max depth reached, ignoring particle count limit\n");
        while(IN_BOUNDS(i)){
          fifo.push(++i);
          DEBUG(i<<"[Override] pushed into fifo: " << fifo.size() <<"\n");
          std::cerr<<"Overriding fifo size to"<<fifo.size()<<std::endl;
        }
      }
    }

    Node leaf;
    leaf.start = fifo.front(); // index of the first particle
    leaf.count = 0; // number of particles, since the particles are sorted, all particles up to start+count-1 are contained
    // Add all particles that are in bounds to the node
    // fifo.size() > 0 must be checked BEFORE IN_BOUNDS, to avoid segfault
    while(fifo.size() > 0 && IN_BOUNDS(fifo.front())){
      //DEBUG_D("popped "<<fifo.front()<<" ("<<(std::bitset<63>(mortoncodes[fifo.front()]))<<")" << std::endl, depth);
      fifo.pop();
      leaf.count++;
    }
    // Add empty nodes when at the leaf level, but do not add empty nodes at the node level
    //if(leaf.count>0)
      tree[depth].push_back(leaf);
    //DEBUG_D("Finalizing Leaf with num_planets="<<leaf.count<<std::endl, depth);
    uint_fast32_t x,y,z;
    libmorton::morton3D_64_decode(current_cell_end, x,y,z);
    myvec3 localnode_dim = node_dim*(std::pow(2, 21-depth));
    myvec3 node_pos = myvec3(static_cast<double>(x)*node_dim.x,static_cast<double>(y)*node_dim.y,static_cast<double>(z)*node_dim.z) + boundingbox.min();
    //DEBUG_D(node_pos<<localnode_dim, depth);
    drawcuboids.push_back(DrawableCuboid(minMaxCuboid(node_pos-localnode_dim, node_pos), depth));
    
    // Go to next node (at this depth if possible)
    if(stack[depth] < 7){
      stack[depth]++;
      // Set the bits in position (3*(21-depth-1)) to stack[depth] in current_cell_end
      current_cell_end += 1ul << (63-3*depth);
      //DEBUG_D("Incremented stack="<<stack[depth]<<std::endl, depth);
      DEBUG_BITS(current_cell_end, depth);DEBUG("\n");
      //DEBUG_D("Morton "<<std::bitset<3>(current_cell_end>> (3*(21-depth)))<<" and stack "<<std::bitset<3>(stack[depth])<<std::endl, depth);
      assert(std::bitset<3>(current_cell_end >> (3*(21-depth)))==std::bitset<3>(stack[depth]));
    }else{
      // Go to the next depth
      while(stack[depth] == 7){
        Node node;
        node.setInternal(tree[depth].size()-8); // index of the first child (8 children per node)
        assert(tree[depth].size()%8 == 0);
        depth--;
        tree[depth].push_back(node);
        //assert(std::bitset<3>(current_cell_end >> (3*(21-depth)))==std::bitset<3>(stack[depth]));
        DEBUG_BITS(current_cell_end, depth);DEBUG("Finalizing Node "<<std::endl);
      }
      stack[depth]++;
      current_cell_end += 1ul << (63-3*depth);
    }

    if(depth <= 0){
      // We have finished the tree
      break;
    }
  }

  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Tree building " << elapsed.count()*1e3<<"ms\n";

  std::array<std::vector<myfloat>, depth_max+1> centers_of_massx;
  std::array<std::vector<myfloat>, depth_max+1> centers_of_massy;
  std::array<std::vector<myfloat>, depth_max+1> centers_of_massz;
  std::array<std::vector<myfloat>, depth_max+1> masses;
  for (int d = 0; d <= depth_max; d++)
  {
    //centers_of_mass.emplace_back(std::move(Vectors{tree[d].size()+1, 0.0}));
    centers_of_massx[d].resize(tree[d].size(), 0.0);
    centers_of_massy[d].resize(tree[d].size(), 0.0);
    centers_of_massz[d].resize(tree[d].size(), 0.0);
    masses[d].resize(tree[d].size(), 0.0);
  }
  time1 = std::chrono::high_resolution_clock::now();
  
  // COM and mass calculation
  for (int d = depth_max-1; d >= 0; d--) // unsigned int underflows to max value
  {
    DEBUG_D("COM computation, Depth "<<d<<" has "<<tree[d].size()<<" nodes"<<std::endl, d);
    for (size_t i = 0; i < tree[d].size(); i++)
    {
      auto& currentnode = tree[d][i];
      if(currentnode.isLeaf()){
        DEBUG_D("COM computation for Leaf "<<i<<" at depth "<<d<<" has "<<currentnode.count<<" particles\n", d);
        for (int j = currentnode.start; j < currentnode.start+currentnode.count; j++)
        {
          auto mass_child = particles.m[j];
          centers_of_massx[d][i] += particles.p.x[j] * mass_child;
          centers_of_massy[d][i] += particles.p.y[j] * mass_child;
          centers_of_massz[d][i] += particles.p.z[j] * mass_child;
          masses[d][i] += mass_child;
        }
      }else{
        // TODO: SIMD this, test unrolling
        DEBUG_D("COM computation for Internal Node "<<i<<" with start "<<currentnode.start<<"\n", d);
        for (int j = currentnode.start; j < currentnode.start + 8; j++)
        {
          auto mass_child = masses[d+1][j];
          centers_of_massx[d][i] += centers_of_massx[d+1][j] * mass_child;
          centers_of_massy[d][i] += centers_of_massy[d+1][j] * mass_child;
          centers_of_massz[d][i] += centers_of_massz[d+1][j] * mass_child;
          masses[d][i] += mass_child;
        }
      }
      if(masses[d][i] != 0){
        centers_of_massx[d][i] /= masses[d][i];
        centers_of_massy[d][i] /= masses[d][i];
        centers_of_massz[d][i] /= masses[d][i];

        DEBUG_D("COM for node "<<i<<" at depth "<<d<<" is "<<myvec3(
          centers_of_massx[d][i], centers_of_massy[d][i], centers_of_massz[d][i]
        )<<" with mass "<<masses[d][i]<<"\n", d);
      }
    }
  }
  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "COM computation took " << elapsed.count()*1e3<<"ms\n";

  // Force calculation
  std::array<myfloat, depth_max+1> diagonal2;
  for (int d = 0; d <= depth_max; d++)
  {
    diagonal2[d] = length2(boundingbox.dimension.x/(1<<d), 
                           boundingbox.dimension.y/(1<<d), 
                           boundingbox.dimension.z/(1<<d));
  }

  time1 = std::chrono::high_resolution_clock::now();
  myfloat* accx = acc.x, *accy=acc.y, *accz=acc.z;
  #pragma omp parallel for reduction(+:accx[:particles.size()],accy[:particles.size()],accz[:particles.size()])
  for (size_t i = 0; i < particles.size(); i++)
  {
    recursive_force(particles, tree, &particles.p.x[i], &particles.p.y[i], &particles.p.z[i], 
    masses, centers_of_massx, centers_of_massy, centers_of_massz, 
    &acc.x[i], &acc.y[i], &acc.z[i], diagonal2, 0, 0);
    DEBUG("Particle "<<i<<" has force "<<acc.x[i]<<" "<<acc.y[i]<<" "<<acc.z[i]<<"\n");
  }
  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Force calculation took " << elapsed.count()*1e3<<"ms\n";
    return drawcuboids;
}

std::vector<DrawableCuboid>  stepSimulation(Particles& particles, myfloat dt, double theta) {
  // barnes hut optimized step
  // Buffer for the new state vector of particles
  Vectors acc{particles.size(), 0.0};
  Vectors p2{particles.size()};

  auto drawcubes = bh_superstep(particles, particles.size(), acc);
  // We have constructed the tree, use it to efficiently compute the
  // accelerations. Far away particles get grouped and their contribution is
  // approximated using their center of mass (Barnes Hut algorithm)

  // Use velocity verlet algorithm to update the particle positions and
  // velocities
  // https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
#pragma omp parallel for //simd
  for (size_t i = 0; i < particles.size(); i++) {
    // First update the positions
    particles.p.x[i] = particles.p.x[i] + particles.v.x[i] * dt + acc.x[i] * dt * dt / 2.;
    particles.p.y[i] = particles.p.y[i] + particles.v.y[i] * dt + acc.y[i] * dt * dt / 2.;
    particles.p.z[i] = particles.p.z[i] + particles.v.z[i] * dt + acc.z[i] * dt * dt / 2.;
  }

#pragma omp parallel for
  for (size_t i = 0; i < particles.size(); i++) {
    // Then update the velocities using v(t+1) = dt*(a(t) + a(t+dt))/2
    particles.v.x[i] += (acc.x[i] + acc.x[i]) * dt / 2.;
    particles.v.y[i] += (acc.y[i] + acc.y[i]) * dt / 2.;
    particles.v.z[i] += (acc.z[i] + acc.z[i]) * dt / 2.;
  }
  return drawcubes;
}