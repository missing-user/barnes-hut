#include <vector>

// #include "Node.h"
#include "Particle.h"
#include "Rectangle.h"
#include "Cuboid.h"
#include "Tree.h"

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
} // https : // stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c

int main()
{
    // Cuboid base(0., 0., 0. 1., 1., 1.);
    Cuboid base2(0., 0., 0., 1., 1., 1.);
    Tree root(base2);
    root.level = 0;
    std::vector<Particle> points(30);
    for (Particle &p : points)
    {
        p.p.x = fRand(0, 1);
        p.p.y = fRand(0, 1);
        p.p.z = fRand(0, 1);
        p.m = fRand(0, 1);
        root.insertPnt(p);
    }
    root.computeCOM();

    std::cout << root.print();
    return 0;
}