#include "Barneshut.h"
#include "Order.h"
#include <queue>
#include <bitset>
#include <libmorton/morton.h>
#include <xsimd/xsimd.hpp>
#include <chrono>

namespace xs = xsimd;
using b_type = xs::batch<myfloat>;
using b_bool_type = xs::batch_bool<myfloat>;

// #define DEBUG_BUILD
#ifdef DEBUG_BUILD
#define DEBUG_D(x, d)                                         \
  do                                                          \
  {                                                           \
    for (int DEBUG_DEPTH = 0; DEBUG_DEPTH < d; DEBUG_DEPTH++) \
    {                                                         \
      std::cout << " ";                                       \
    }                                                         \
    std::cout << x;                                           \
  } while (0)
#else
#define DEBUG_D(x, depth)
#endif

#define DEBUG(x) DEBUG_D(x, 0)

const int depth_max = 18;
const int leaf_max = 8; // maximum particles per node
void PRINT_BITS_BY3(uint_fast64_t x, int separator)
{
  for (size_t i = 1; i <= depth_max - 1; i++)
  {
    auto prefix = std::bitset<3>(x >> (63 - 3 * i));
    DEBUG(prefix);
    if (i == separator)
      DEBUG(":");
    else
      DEBUG(" ");
  }
}

void DEBUG_BITS(uint_fast64_t end, int depth)
{
  DEBUG_D("", depth);
  PRINT_BITS_BY3(end, depth);
  auto extsuf = std::bitset<63 - 3 * depth_max - 1>(end);
  // assert(extsuf == std::bitset<63-3*depth_max>(std::numeric_limits<uint_fast64_t>::max())); //unused bottom part must be 1
  // assert((end>>63) == 0); // Top bit must always be 0
}

// Due to morton order we know that all current particles.pos are > node_pos, only check the upper bounds
// Also contains the array in bounds check (prevent segfault), i+1 to avoid overflow on unsigned int
#define IN_BOUNDS(i) ((i < particles.size()) && (mortoncodes[i] <= current_max_morton))
#define LAST_ITERATION (i == particles.size() - 1 && current_node_count <= 1)

struct Node
{
  // Check Figure 3. in https://www.tabellion.org/et/paper11/OutOfCorePBGI.pdf
  int start; // index of the first particle
             // If internal node, start is the index of the first child. All Nodes up to start+8 are children
  int count; // number of particles, since the particles are sorted, the last particle is start+count-1
             // Negative count indicates that this is an internal node

  bool isLeaf() const { return count >= 0; }
  void setInternal(int first_child)
  {
    start = first_child;
    count = -1;
  }
};

typedef std::array<std::vector<Node>, depth_max> Tree;

struct CentersOfMass
{
  std::array<vector_type, depth_max> x;
  std::array<vector_type, depth_max> y;
  std::array<vector_type, depth_max> z;
  std::array<vector_type, depth_max> m;
};

CentersOfMass compute_centers_of_mass(const Particles &particles, const Tree &tree)
{
  CentersOfMass com;
  // Initialize the arrays
  for (int d = 0; d < depth_max; d++)
  {
    com.x[d].resize(tree[d].size());
    com.y[d].resize(tree[d].size());
    com.z[d].resize(tree[d].size());
    com.m[d].resize(tree[d].size());
  }

  // No gain from parallelization, start one above the leaf level
  for (int d = depth_max - 2; d >= 0; d--) // unsigned int underflows to max value
  {
    DEBUG_D("COM computation, Depth " << d << " has " << tree[d].size() << " nodes" << std::endl, d);
    for (size_t i = 0; i < tree[d].size(); i++)
    {
      auto &currentnode = tree[d][i];
      myfloat comx = 0, comy = 0, comz = 0, comm = 0;
      if (currentnode.isLeaf())
      {
        DEBUG_D("COM computation for Leaf " << i << " at depth " << d << " has " << currentnode.count << " particles\n", d);
#pragma omp simd
        for (int j = currentnode.start; j < currentnode.start + currentnode.count; j++)
        {
          auto mass_child = particles.m[j];
          comx += particles.p.x[j] * mass_child;
          comy += particles.p.y[j] * mass_child;
          comz += particles.p.z[j] * mass_child;
          comm += mass_child;
        }
      }
      else
      {
        DEBUG_D("COM computation for Internal Node " << i << " with start " << currentnode.start << "\n", d);
#pragma omp simd safelen(8)
        for (int j = currentnode.start; j < currentnode.start + 8; j++)
        {
          auto mass_child = com.m[d + 1][j];
          comx += com.x[d + 1][j] * mass_child;
          comy += com.y[d + 1][j] * mass_child;
          comz += com.z[d + 1][j] * mass_child;
          comm += mass_child;
        }
      }
      if (comm != 0)
      {
        comx /= comm;
        comy /= comm;
        comz /= comm;

        DEBUG_D("COM for node " << i << " at depth " << d << " is " << myvec3(comx, comy, comz) << " with mass " << comm << "\n", d);
      }
      // Accumulating into temporary variables and then doing
      // a single array access is much faster. Cache contention?
      com.x[d][i] = comx;
      com.y[d][i] = comy;
      com.z[d][i] = comz;
      com.m[d][i] = comm;
    }
  }
  return com;
}

template <typename T>
inline auto isApproximationValid(T dx, T dy, T dz, myfloat theta2, myfloat diagonal2)
{
  return static_cast<T>(diagonal2) < theta2 * length2(dx, dy, dz);
}

// Template implementation works with float, double and xsimd::batch of float or double
template <typename T>
void recursive_force(
    T *__restrict accx, T *__restrict accy, T *__restrict accz,
    const T x, const T y, const T z,
    const Particles &__restrict particles, const Tree &__restrict tree, const CentersOfMass &__restrict com,
    const std::array<myfloat, depth_max> &diagonal2, int depth, int start, myfloat theta2)
{
  auto &node = tree[depth][start];
  T dvx = 0, dvy = 0, dvz = 0;
  if (node.isLeaf())
  {
    bruteForceAcc(&dvx, &dvy, &dvz,
                  particles.p.x.data() + node.start, 
                  particles.p.y.data() + node.start, 
                  particles.p.z.data() + node.start,
                  x, y, z, particles.m.data() + node.start, node.count);
  }
  else
  {
    auto dx = com.x[depth][start] - x;
    auto dy = com.y[depth][start] - y;
    auto dz = com.z[depth][start] - z;

    // Multipole acceptance criterion
    auto mask = isApproximationValid<T>(dx, dy, dz, theta2, diagonal2[depth]);
    bool approximate;
    if constexpr (xs::is_batch<T>::value)
      approximate = xs::all(mask); // All elements of the batch must be true
    else
      approximate = mask; // Scalar case, the mask is a single bool

    // If approximating, we can use the center of mass of the node, otherwise we need to go deeper
    if (approximate)
    {
      accelFunc<T>(accx, accy, accz, dx, dy, dz, com.m[depth][start]);
    }
    else
    {
      for (int it = node.start; it < node.start + 8; it++)
      {
        recursive_force<T>(&dvx, &dvy, &dvz, x, y, z, particles, tree, com, diagonal2, depth + 1, it, theta2);
      }
    }
  }
  *accx += dvx;
  *accy += dvy;
  *accz += dvz;
}

std::array<myfloat, depth_max> precompute_diagonals(const myfloat top_level_diagonal2)
{
  std::array<myfloat, depth_max> diagonal2;
  for (int d = 0; d < depth_max; d++)
  {
    // The diagonal the bounding boxes at a given depth are half of the previous level.
    // Since we are using the squared distance, we need to divide by 4
    diagonal2[d] = top_level_diagonal2 / (1 << (d * 2));
  }
  return diagonal2;
}

void compute_accelerations(Vectors &acc, const Particles &particles, const Tree &tree, const CentersOfMass &com,
                           myfloat dt, myfloat theta2, const Cuboid &boundingbox)
{
  // Precompute the squared diagonals
  auto diagonal2 = precompute_diagonals(boundingbox.diagonal2);

  // Force calculation
  auto vectorized_size = particles.size() - (particles.size() % b_type::size);
#pragma omp parallel
  {
#pragma omp for schedule(dynamic, 64) nowait
    for (int i = 0; i < vectorized_size; i += b_type::size)
    {
      b_type dvx = 0, dvy = 0, dvz = 0;
      auto x = b_type::load_unaligned(&particles.p.x[i]);
      auto y = b_type::load_unaligned(&particles.p.y[i]);
      auto z = b_type::load_unaligned(&particles.p.z[i]);
      recursive_force(&dvx, &dvy, &dvz, x, y, z, particles, tree, com, diagonal2, 0, 0, theta2);
      dvx.store_unaligned(&acc.x[i]);
      dvy.store_unaligned(&acc.y[i]);
      dvz.store_unaligned(&acc.z[i]);
    }

// Scalar remainder loop
#pragma omp single
    for (size_t i = vectorized_size; i < particles.size(); i++)
    {
      myfloat dvx = 0, dvy = 0, dvz = 0;
      recursive_force(&dvx, &dvy, &dvz, particles.p.x[i], particles.p.y[i], particles.p.z[i],
                      particles, tree, com, diagonal2, 0, 0, theta2);
      acc.x[i] = dvx;
      acc.y[i] = dvy;
      acc.z[i] = dvz;
    }
  }
}

Tree build_tree(Particles &particles, const Cuboid &boundingbox)
{
  auto mortoncodes = computeMortonCodes(particles, boundingbox);
  // assert(std::is_sorted(mortoncodes.begin(), mortoncodes.end()));
  Tree tree;
  const int internal_depth_max = depth_max - 1;
  /* Since we know the particles are sorted in Morton order, we can use out of core
   * Tree construction as described in http://www.thesalmons.org/john/pubs/siam97/salmon.pdf
   * by only keeping track of the group currently being constructed and finalizing it
   * once the first particle outside the group is found. (morton order ensures there will be
   * no subsequent particles in the group)
   */
  size_t i = 0;
  int depth = 0;                                                                              // 0 = root
  uint64_t current_max_morton = 0x7FFFFFFFFFFFFFFF; // Last morton code that is still in the current cell
  short stack[depth_max];
  std::memset(stack, 0, sizeof(stack));
  // FIFO queue: all particle indices between start and start+count are designated for the current node
  int current_node_start = 0, current_node_count = 1;
  std::vector<DrawableCuboid> drawcuboids{};

  while (i < particles.size())
  {
    // Add the first leaf_max+1 particles to the queue
    auto prev_node_count = current_node_count;
    current_node_count = std::min(static_cast<int>(particles.size()) - current_node_start, leaf_max + 1);
    i += current_node_count - prev_node_count;

    // Find the first depth level that does not contain all particles in the fifo
    // due to morton order we know that the last particle (i) will be the first to be out of range
    // We have to handle the special case when only one particle is left, to avoid infinite deepening of the tree
    while (IN_BOUNDS(i) && depth < internal_depth_max && !LAST_ITERATION)
    {
      depth++;
      stack[depth] = 0;
      current_max_morton -= static_cast<uint64_t>(7) << static_cast<uint64_t>(63 - 3 * depth);
      DEBUG_BITS(current_max_morton, depth);
      DEBUG_D(" Subdivide, i=" << i << " is still in range"
                               << "\n",
              depth);

      if (depth >= internal_depth_max)
      {
        // Max depth reached, ignore the particle count limit and fill the node
        // To avoid infinite deepening of the tree when particles are on the same position
        DEBUG("Max depth reached, ignoring particle count limit\n");
        while (IN_BOUNDS(i))
        {
          i++;
          current_node_count++;
          DEBUG(i << "[Override] pushed into fifo: " << current_node_count << "\n");
          std::cerr << "Overriding fifo size to" << current_node_count << "\n";
        }
      }
    }
    Node leaf;
    leaf.start = current_node_start; // index of the first particle
    leaf.count = 0;                  // number of particles, since the particles are sorted, all particles up to start+count-1 are contained
    // Add all particles that are in bounds to the node
    // current_node_count > 0 must be checked BEFORE IN_BOUNDS, to avoid segfault
    while (current_node_count > 0 && IN_BOUNDS(current_node_start))
    {
      DEBUG_D("popped " << current_node_start << " (" << (std::bitset<63>(mortoncodes[current_node_start])) << ")" << std::endl, depth);
      current_node_start++;
      current_node_count--;
      leaf.count++;
    }
    tree[depth].push_back(leaf);
    // DEBUG_D("Finalizing Leaf with num_planets="<<leaf.count<<std::endl, depth);
    uint_fast32_t x, y, z;
    libmorton::morton3D_64_decode(current_max_morton, x, y, z);

    // Go to next node (at this depth if possible)
    if (stack[depth] < 7)
    {
      stack[depth]++;
      // Set the bits in position (3*(21-depth-1)) to stack[depth] in current_max_morton
      current_max_morton += static_cast<uint64_t>(1) << static_cast<uint64_t>(63 - 3 * depth);
      DEBUG_BITS(current_max_morton, depth);
      DEBUG("\n");
      // assert(std::bitset<3>(current_max_morton >> (3*(21-depth)))==std::bitset<3>(stack[depth]));
    }
    else
    {
      // Go to the next depth
      while (stack[depth] == 7)
      {
        Node node;
        node.setInternal(tree[depth].size() - 8); // index of the first child (8 children per node)
        // assert(tree[depth].size()%8 == 0);
        depth--;
        tree[depth].push_back(node);
        // assert(std::bitset<3>(current_max_morton >> (3*(21-depth)))==std::bitset<3>(stack[depth]));
        DEBUG_BITS(current_max_morton, depth);
        DEBUG("Finalizing Node " << std::endl);
      }
      stack[depth]++;
      current_max_morton += static_cast<uint64_t>(1) << static_cast<uint64_t>(63 - 3 * depth);
    }

    if (depth <= 0)
    {
      // We have finished the tree
      break;
    }
  }
  return tree;
}

void bh_superstep(Vectors &acc, Particles &particles, size_t count, myfloat dt, myfloat theta2)
{
    auto time0 = std::chrono::high_resolution_clock::now();
  auto boundingbox = bounding_box(particles.p, count);
    auto time1 = std::chrono::high_resolution_clock::now();
  computeAndOrder(particles, boundingbox);
    auto time2 = std::chrono::high_resolution_clock::now();
  auto tree = build_tree(particles, boundingbox);
    auto time3 = std::chrono::high_resolution_clock::now();
  auto com = compute_centers_of_mass(particles, tree);
    auto time4 = std::chrono::high_resolution_clock::now();
  compute_accelerations(acc, particles, tree, com, dt, theta2, boundingbox);
    auto time5 = std::chrono::high_resolution_clock::now();
  std::cout<< "Times: "<<
  std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count()<<" ms, "<<
  std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count()<<" ms, "<<
  std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count()<<" ms, "<<
  std::chrono::duration_cast<std::chrono::milliseconds>(time4 - time3).count()<<" ms, "<<
  std::chrono::duration_cast<std::chrono::milliseconds>(time5 - time4).count()<<" ms\n";
}

std::vector<DrawableCuboid> draw_approximations(
    const myfloat x, const myfloat y, const myfloat z,
    const Particles &__restrict particles,
    const Tree &__restrict tree,
    const CentersOfMass &__restrict com,
    const std::array<myfloat, depth_max> &diagonal2,
    const Cuboid &boundingbox,
    int depth, int start, myfloat theta2)
{
  std::vector<DrawableCuboid> draw;
  auto &node = tree[depth][start];
  if (!node.isLeaf())
  {
    myfloat dx = com.x[depth][start] - x;
    myfloat dy = com.y[depth][start] - y;
    myfloat dz = com.z[depth][start] - z;
    if (isApproximationValid(dx, dy, dz, theta2, diagonal2[depth]))
    {
      // Compute the COM force
      draw.push_back(DrawableCuboid{{myvec3(com.x[depth][start], com.y[depth][start], com.z[depth][start]),
                                     boundingbox.dimension / static_cast<myfloat>(1 << depth)},
                                    depth});
    }
    else
    {
      for (int it = node.start; it < node.start + 8; it++)
      {
        auto tmp = draw_approximations(x, y, z,
                                       particles, tree, com, diagonal2, boundingbox, depth + 1, it, theta2);
        draw.insert(draw.end(), tmp.begin(), tmp.end());
      }
    }
  }
  return draw;
}

debug_information bh_superstep_debug(Particles &particles, size_t count, myfloat theta2, myvec3 position)
{
  debug_information info;
  auto boundingbox = bounding_box(particles.p, count);
  computeAndOrder(particles, boundingbox);
  auto tree = build_tree(particles, boundingbox);

  info.depth = 0;
  for (auto &depth : tree)
  {
    if (depth.size() > 0)
    {
      info.depth++;
    }
    for (auto &node : depth)
    {
      if (node.isLeaf())
      {
        info.max_particles_in_leaf = std::max(info.max_particles_in_leaf, node.count);
      }
    }
  }

  auto com = compute_centers_of_mass(particles, tree);
  auto diagonal2 = precompute_diagonals(boundingbox.diagonal2);

  info.debug_boxes = draw_approximations(position.x, position.y, position.z, particles, tree, com, diagonal2, boundingbox, 0, 0, theta2);
  return info;
}

void stepSimulation(Particles &particles, myfloat dt, myfloat theta2)
{
  // Far away particles get grouped and their contribution is
  // approximated using their center of mass (Barnes Hut algorithm)
  Vectors acc{particles.size()};
  bh_superstep(acc, particles, particles.size(), dt, theta2);

  // Use velocity verlet algorithm to update the particle positions and
  // velocities
  // https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
  // Completely memory bound, no use in parallelizing
#pragma omp simd
  for (size_t i = 0; i < particles.size(); i++)
  {
    particles.v.x[i] += acc.x[i] * dt;
    particles.v.y[i] += acc.y[i] * dt;
    particles.v.z[i] += acc.z[i] * dt;

    particles.p.x[i] += particles.v.x[i] * dt;
    particles.p.y[i] += particles.v.y[i] * dt;
    particles.p.z[i] += particles.v.z[i] * dt;
  }
}