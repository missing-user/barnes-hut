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

std::pair<myvec3, myvec3> bounding_box(const std::vector<Particle> &particles)
{
    myvec3 bmin = particles.at(0).p;
    myvec3 bmax = particles.at(0).p;

    for (const auto &p : particles)
    {
        bmin = glm::min(bmin, p.p);
        bmax = glm::max(bmax, p.p);
    }
    // Bounding box containment is defined as a half open interval: [min, max)
    // To actually make the particles contained within bmax, we need to add a small offset
    bmax += myvec3{1e-6,1e-6,1e-6};

    return {bmin, bmax};
}

Particle operator+(const Particle &P1, const Particle &P2)
{
    Particle p;
    p.p = P1.p + P2.p;
    // p.v = P1.v + P2.v; // Only used for center of mass calc, no v
    p.m = P1.m + P2.m;
    return p;
}
#endif