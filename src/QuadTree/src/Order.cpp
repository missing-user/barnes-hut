#include "Order.h"
#include <chrono>

void reorder(std::vector<Particle> &data,
             std::vector<std::size_t> const &order) {
  // Reorder function from
  // https://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices
  std::vector<Particle> tmp; // create an empty vector
  tmp.reserve(data.size());  // ensure memory and avoid moves in the vector
  for (std::size_t i = 0; i < order.size(); ++i) {
    tmp.push_back(data[order[i]]);
  }
  data.swap(tmp); // swap vector contents
}

void computeAndOrder(std::vector<Particle> &particles) {
  std::cout << "Starting reorder\n";
  auto begin = std::chrono::steady_clock::now();

  // Assign ids to particles for reordering
  for (size_t i = 0; i < particles.size(); i++) {
    particles[i].id = i;
  }

  Tree mytree(particles);
  std::cout << "Reorder tree constructed\n";
  reorder(particles, mytree.DFS());

  auto end = std::chrono::steady_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
  std::cout << "Reorder finished " << elapsed.count() << " ms\n";
}
