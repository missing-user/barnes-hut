#include "Barneshut.h"
#include "Order.h"
#include <queue>
#include <bitset>
#include <libmorton/morton.h>
#include <chrono>

//#define LOG_TIME
//#define DEBUG_BUILD
#ifdef DEBUG_BUILD
#define DEBUG_D(x, d) do { \
  for (int DEBUG_DEPTH = 0; DEBUG_DEPTH < d; DEBUG_DEPTH++){std::cout << " ";}\
  std::cout << x; } while (0)
#else
#define DEBUG_D(x, depth)
#endif

#define DEBUG(x) DEBUG_D(x, 0)

const int depth_max = 14;
const int leaf_max = 8; // maximum particles per node
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
  //assert(extsuf == std::bitset<63-3*depth_max>(std::numeric_limits<uint_fast64_t>::max())); //unused bottom part must be 1
  //assert((end>>63) == 0); // Top bit must always be 0
}

// Due to morton order we know that all current particles.pos are > node_pos, only check the upper bounds
// Also contains the array in bounds check (prevent segfault), i+1 to avoid overflow on unsigned int
#define IN_BOUNDS(i) ((i < particles.size()) && (mortoncodes[i] <= current_cell_end))
#define LAST_ITERATION (i==particles.size()-1 && current_node_count<= 1)

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

void compute_centers_of_mass(const Particles& particles, const std::array<std::vector<Node>, depth_max+1>& tree, 
                             std::array<std::vector<myfloat>, depth_max+1>& centers_of_massx, 
                             std::array<std::vector<myfloat>, depth_max+1>& centers_of_massy, 
                             std::array<std::vector<myfloat>, depth_max+1>& centers_of_massz, 
                             std::array<std::vector<myfloat>, depth_max+1>& cumulative_mass)
{
  // Initialize the arrays
  for (int d = 0; d <= depth_max; d++)
  {
    //centers_of_mass.emplace_back(std::move(Vectors{tree[d].size()+1, 0.0}));
    centers_of_massx[d].resize(tree[d].size());
    centers_of_massy[d].resize(tree[d].size());
    centers_of_massz[d].resize(tree[d].size());
    cumulative_mass[d].resize(tree[d].size());
  }

  // No gain from parallelization
  for (int d = depth_max-1; d >= 0; d--) // unsigned int underflows to max value
  {
    DEBUG_D("COM computation, Depth "<<d<<" has "<<tree[d].size()<<" nodes"<<std::endl, d);
    for (size_t i = 0; i < tree[d].size(); i++)
    {
      auto& currentnode = tree[d][i];
      myfloat comx=0, comy=0, comz=0, comm=0;
      if(currentnode.isLeaf()){
        DEBUG_D("COM computation for Leaf "<<i<<" at depth "<<d<<" has "<<currentnode.count<<" particles\n", d);
        for (int j = currentnode.start; j < currentnode.start+currentnode.count; j++)
        {
          auto mass_child = particles.m[j];
          comx += particles.p.x[j] * mass_child;
          comy += particles.p.y[j] * mass_child;
          comz += particles.p.z[j] * mass_child;
          comm += mass_child;
        }
      }else{
        DEBUG_D("COM computation for Internal Node "<<i<<" with start "<<currentnode.start<<"\n", d);
        #pragma omp simd safelen(8)
        for (int j = currentnode.start; j < currentnode.start + 8; j++)
        {
          auto mass_child = cumulative_mass[d+1][j];
          comx += centers_of_massx[d+1][j] * mass_child;
          comy += centers_of_massy[d+1][j] * mass_child;
          comz += centers_of_massz[d+1][j] * mass_child;
          comm += mass_child;
        }
      }
      if(comm != 0){
        comx /= comm;
        comy /= comm;
        comz /= comm;

        DEBUG_D("COM for node "<<i<<" at depth "<<d<<" is "<<myvec3(comx, comy, comz)<<" with mass "<<comm<<"\n", d);
      }
      // Accumulating into temporary variables and then doing 
      // a single array access is much faster. Cache contention?
      centers_of_massx[d][i] = comx;
      centers_of_massy[d][i] = comy;
      centers_of_massz[d][i] = comz;
      cumulative_mass[d][i] = comm;
    }
  }
}
#pragma omp declare simd
inline bool isApproximationValid(myfloat dx,myfloat dy,myfloat dz, myfloat theta2, myfloat diagonal2)
{
  return diagonal2 < theta2 * length2(dx,dy,dz);
}

#pragma omp declare simd linear(accx, accy, accz) uniform(depth, diagonal2, start,x,y,z)
void recursive_force(
  const Particles& particles,
  const std::array<std::vector<Node>, depth_max+1>& tree,
  const myfloat x, const myfloat y, const myfloat z, 
  const std::array<std::vector<myfloat>, depth_max+1>& commass, 
  const std::array<std::vector<myfloat>, depth_max+1>& comx,
  const std::array<std::vector<myfloat>, depth_max+1>& comy,
    const std::array<std::vector<myfloat>, depth_max+1>& comz,
     myfloat* __restrict accx,  myfloat* __restrict accy,  myfloat* __restrict accz, 
    const std::array<myfloat, depth_max+1>& diagonal2,
    int depth, int start, myfloat theta2){
  DEBUG_D("recursive_force enter with start "<<start<<"\n", depth);
  auto& node = tree[depth][start];
  if(node.isLeaf()){
    DEBUG_D("recursive_force leaf with "<<node.count<<" particles\n", depth);
    // Compute the force
    myfloat dvx = 0, dvy = 0, dvz = 0;
    #pragma omp simd
    for (int i = node.start; i < node.start + node.count; i++)
    {
      myfloat diffx = particles.p.x[i] - x;
      myfloat diffy = particles.p.y[i] - y;
      myfloat diffz = particles.p.z[i] - z;

      constexpr myfloat softening_param = 0.025;
      myfloat r2 = length2(diffx, diffy, diffz)+softening_param;
      myfloat mOverDist3 = particles.m[i] / (r2 * std::sqrt(r2));
      
      dvx += diffx * mOverDist3;
      dvy += diffy * mOverDist3;
      dvz += diffz * mOverDist3;
    }

    *accx += dvx;
    *accy += dvy;
    *accz += dvz;
  }else{
    myfloat dx = comx[depth][start] - x;
    myfloat dy = comy[depth][start] - y;
    myfloat dz = comz[depth][start] - z;
    if(isApproximationValid(dx,dy,dz, theta2, diagonal2[depth])){
      DEBUG_D("recursive_force approximated with COM "<<start<<" = "<<myvec3(
        comx[depth][start], comy[depth][start], comz[depth][start]
      )<<"\n", depth);
      // Compute the COM force
      constexpr myfloat softening_param = 0.025;
      myfloat r2 = length2(dx, dy, dz)+softening_param;
      myfloat mOverDist3 = commass[depth][start] / (r2 * std::sqrt(r2));
      
      *accx += dx * mOverDist3;
      *accy += dy * mOverDist3;
      *accz += dz * mOverDist3;
    }else{
      #pragma omp simd
      for (auto it = node.start; it < node.start+8; it++)
      {
        recursive_force(particles, tree, x,y,z, commass, comx, comy, comz,
        accx, accy, accz, diagonal2, depth+1, it, theta2);
      }
    }
  }

  DEBUG_D("recursive_force exit depth "<<depth<<"\n", depth);
}

void compute_accelerations(const Particles& particles, std::array<std::vector<Node>, depth_max+1> tree,
                          std::array<std::vector<myfloat>, depth_max+1>& centers_of_massx, 
                          std::array<std::vector<myfloat>, depth_max+1>& centers_of_massy, 
                          std::array<std::vector<myfloat>, depth_max+1>& centers_of_massz, 
                          std::array<std::vector<myfloat>, depth_max+1>& cumulative_mass,
                          myfloat dt, myfloat theta2, const Cuboid& boundingbox){
  std::array<myfloat, depth_max+1> diagonal2;
  for (int d = 0; d <= depth_max; d++)
  {
    // The diagonal the bounding boxes at a given depth are half of the previous level.
    // Since we are using the squared distance, we need to divide by 4!s
    diagonal2[d] = boundingbox.diagonal2/(1<<(d*2));
  }

  // Force calculation
#pragma omp parallel for schedule(dynamic,64) 
  for (size_t i = 0; i < particles.size(); i++) 
  {
    myfloat dvx = 0, dvy = 0, dvz = 0;
    recursive_force(particles, tree, particles.p.x[i], particles.p.y[i], particles.p.z[i], 
    cumulative_mass, centers_of_massx, centers_of_massy, centers_of_massz, 
    &dvx, &dvy, &dvz, diagonal2, 0, 0, theta2);
    particles.v.x[i] += dvx*dt;
    particles.v.y[i] += dvy*dt;
    particles.v.z[i] += dvz*dt;
  }
}

std::array<std::vector<Node>, depth_max+1> build_tree(Particles& particles, const Cuboid& boundingbox){
  auto mortoncodes = computeMortonCodes(particles, boundingbox);
  assert(std::is_sorted(mortoncodes.begin(), mortoncodes.end()));
  std::array<std::vector<Node>, depth_max+1> tree;
  /* Since we know the particles are sorted in Morton order, we can use out of core
  * Tree construction as described in http://www.thesalmons.org/john/pubs/siam97/salmon.pdf
  * by only keeping track of the group currently being constructed and finalizing it 
  * once the first particle outside the group is found. (morton order ensures there will be 
  * no subsequent particles in the group)
  */
  size_t i = 0;
  int depth = 0; // 0 = root
  uint_fast64_t current_cell_end = std::numeric_limits<uint_fast64_t>::max() - (1ul << 63); // First morton code that is not in the current cell
  short stack[depth_max+1];
  std::memset(stack, 0, sizeof(stack));
  // FIFO queue: all particle indices between start and start+count are designated for the current node
  size_t current_node_start = 0, current_node_count = 1;
  std::vector<DrawableCuboid>  drawcuboids{};

  while(i<particles.size()){
    // Add the first leaf_max+1 particles to the queue
    auto prev_node_count = current_node_count;
    current_node_count = std::min(particles.size()-current_node_start, static_cast<size_t>(leaf_max)+1);
    i += current_node_count-prev_node_count;

    // Find the first depth level that does not contain all particles in the fifo
    // due to morton order we know that the last particle (i) will be the first to be out of range
    // We have to handle the special case when only one particle is left, to avoid infinite deepening of the tree 
    while (IN_BOUNDS(i) && depth < depth_max && !LAST_ITERATION)
    {
      depth++;
      stack[depth] = 0;
      current_cell_end -= 7ul << (63-3*depth);
      DEBUG_BITS(current_cell_end, depth);
      DEBUG_D(" Subdivide, i="<<i<<" is still in range"<<"\n", depth);

      if(depth >= depth_max){
        // Max depth reached, ignore the particle count limit and fill the node 
        // To avoid infinite deepening of the tree when particles are on the same position
        DEBUG("Max depth reached, ignoring particle count limit\n");
        while(IN_BOUNDS(i)){
          i++;
          current_node_count++;
          DEBUG(i<<"[Override] pushed into fifo: " << current_node_count <<"\n");
          std::cerr<<"Overriding fifo size to"<<current_node_count<<std::endl;
        }
      }
    }
    Node leaf;
    leaf.start = current_node_start; // index of the first particle
    leaf.count = 0; // number of particles, since the particles are sorted, all particles up to start+count-1 are contained
    // Add all particles that are in bounds to the node
    // current_node_count > 0 must be checked BEFORE IN_BOUNDS, to avoid segfault
    while(current_node_count > 0 && IN_BOUNDS(current_node_start)){
      DEBUG_D("popped "<<current_node_start<<" ("<<(std::bitset<63>(mortoncodes[current_node_start]))<<")" << std::endl, depth);
      current_node_start++;
      current_node_count--;
      leaf.count++;
    }
    tree[depth].push_back(leaf);
    //DEBUG_D("Finalizing Leaf with num_planets="<<leaf.count<<std::endl, depth);
    uint_fast32_t x,y,z;
    libmorton::morton3D_64_decode(current_cell_end, x,y,z);
    
    // Go to next node (at this depth if possible)
    if(stack[depth] < 7){
      stack[depth]++;
      // Set the bits in position (3*(21-depth-1)) to stack[depth] in current_cell_end
      current_cell_end += 1ul << (63-3*depth);
      DEBUG_BITS(current_cell_end, depth);DEBUG("\n");
      //assert(std::bitset<3>(current_cell_end >> (3*(21-depth)))==std::bitset<3>(stack[depth]));
    }else{
      // Go to the next depth
      while(stack[depth] == 7){
        Node node;
        node.setInternal(tree[depth].size()-8); // index of the first child (8 children per node)
        //assert(tree[depth].size()%8 == 0);
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
  return tree;
}

std::vector<DrawableCuboid> bh_superstep(Particles& particles, size_t count, myfloat dt, myfloat theta2){
  auto time1 = std::chrono::high_resolution_clock::now();
  auto boundingbox = bounding_box(particles.p, count);
  
#ifdef LOG_TIME
  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Bounding Box calculation took " << elapsed.count()*1e3<<"ms "<< std::endl;
  time1 = std::chrono::high_resolution_clock::now();
#endif
  computeAndOrder(particles, boundingbox);
#ifdef LOG_TIME
elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Index reordering took " << elapsed.count()*1e3<<"ms"<< "\n";
  time1 = std::chrono::high_resolution_clock::now();
#endif
#ifdef LOG_TIME
  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Recomputing morton codes took " << elapsed.count()*1e3<<"ms "<< std::endl;
  time1 = std::chrono::high_resolution_clock::now();
#endif
  
  DEBUG("Particles: \n");
  for (size_t i = 0; i < count; i++)
  {
    DEBUG("("<<particles.p.x[i] << "," << particles.p.y[i] << "," << particles.p.z[i] <<"), "<<std::bitset<63>(mortoncodes[i])<<"\n");
  }
  
  // Create the tree
  auto tree = build_tree(particles, boundingbox);
  

#ifdef LOG_TIME
  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Tree building " << elapsed.count()*1e3<<"ms\n";
  time1 = std::chrono::high_resolution_clock::now();
#endif
  std::array<std::vector<myfloat>, depth_max+1> centers_of_massx;
  std::array<std::vector<myfloat>, depth_max+1> centers_of_massy;
  std::array<std::vector<myfloat>, depth_max+1> centers_of_massz;
  std::array<std::vector<myfloat>, depth_max+1> masses;
  
  // COM and mass calculation
  compute_centers_of_mass(particles, tree, centers_of_massx, centers_of_massy, centers_of_massz, masses);
#ifdef LOG_TIME
  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "COM computation took " << elapsed.count()*1e3<<"ms\n";
  time1 = std::chrono::high_resolution_clock::now();
#endif
  compute_accelerations(particles, tree, centers_of_massx, centers_of_massy, centers_of_massz, masses, 
                        dt, theta2, boundingbox);

#ifdef LOG_TIME
  elapsed = std::chrono::high_resolution_clock::now() - time1;
  std::cout << "Force calculation took " << elapsed.count()*1e3<<"ms\n";
#endif
  return {};
}

std::vector<DrawableCuboid>  stepSimulation(Particles& particles, myfloat dt, myfloat theta2) {
  // barnes hut optimized step
  

  auto drawcubes = bh_superstep(particles, particles.size(), dt, theta2);
  // We have constructed the tree, use it to efficiently compute the
  // accelerations. Far away particles get grouped and their contribution is
  // approximated using their center of mass (Barnes Hut algorithm)

  // Use velocity verlet algorithm to update the particle positions and
  // velocities
  // https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
#pragma omp simd
  for (size_t i = 0; i < particles.size(); i++)
  {
    particles.p.x[i] += particles.v.x[i] * dt;
    particles.p.y[i] += particles.v.y[i] * dt;
    particles.p.z[i] += particles.v.z[i] * dt;
  }

  return drawcubes;
}