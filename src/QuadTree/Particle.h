#ifndef PARTICLE
#define PARTICLE

#include <sstream>
#include <vector>

#include <glm/glm.hpp>    // pvec3
#include <glm/gtx/io.hpp> // Allows us to easily std::cout << pvec3;

typedef glm::dvec3 myvec3; // We can easily switch the entire implementation to float precision by adjusting these two variables
typedef double myfloat;

struct Particle
{
public:
    myvec3 p{0, 0, 0};
    myvec3 v{0, 0, 0};
    myfloat m{0};
};

inline std::ostream &operator<<(std::ostream &out, const Particle &p)
{
    // Helper function for the particle class, so we can print and debug it easily
    return out << p.p << "\tv=" << p.v << "\tm=" << p.m;
}

inline std::ostream &operator<<(std::ostream &out, const std::vector<Particle> &particles)
{
    // Print the resulting values of all particles in a .csv format
    for (auto &p : particles)
    {
        out << p.p[0] << "," << p.p[1] << "," << p.p[2] << ",";
    }
    return out;
}

std::pair<myvec3, myvec3> bounding_box(const std::vector<Particle> &);
// std::pair<myvec3, myvec3> bounding_box(std::vector<Particle, std::allocator<Particle>> const &);

Particle operator+(const Particle &P1, const Particle &P2);
#endif