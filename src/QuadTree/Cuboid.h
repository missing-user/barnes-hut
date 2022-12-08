#ifndef CUBOID
#define CUBOID

#include "Particle.h"

class Cuboid
{
public:
    myvec3 min_extent;
    myvec3 max_extent;
    Cuboid() {}
    Cuboid(myvec3 min, myvec3 max)
    {
        min_extent = min;
        max_extent = max;
    }
    bool contains(const Particle &P)
    {        
        return glm::all(glm::greaterThanEqual(P.p, min_extent)) && glm::all(glm::lessThan(P.p, max_extent));
    }
    std::vector<Cuboid> subdivide()
    {
        std::vector<Cuboid> subcuboids(8);
        myfloat x1=min_extent.x,x2=max_extent.x,
                y1=min_extent.y,y2=max_extent.y,
                z1=min_extent.z,z2=max_extent.z;

        auto middle = (min_extent+max_extent)/2.0;

        subcuboids[0] = Cuboid(min_extent, middle);
        subcuboids[1] = Cuboid(myvec3((x1 + x2) / 2, y1, z1),myvec3(x2, (y1 + y2) / 2, (z1 + z2) / 2));
        subcuboids[2] = Cuboid(myvec3(x1, (y1 + y2) / 2, z1),myvec3((x1 + x2) / 2, y2, (z1 + z2) / 2));
        subcuboids[3] = Cuboid(myvec3((x1 + x2) / 2, (y1 + y2) / 2, z1),myvec3(x2, y2, (z1 + z2) / 2));
        subcuboids[4] = Cuboid(myvec3(x1, y1, (z1 + z2) / 2),myvec3((x1 + x2) / 2, (y1 + y2) / 2, z2));
        subcuboids[5] = Cuboid(myvec3((x1 + x2) / 2, y1, (z1 + z2) / 2),myvec3(x2, (y1 + y2) / 2, z2));
        subcuboids[6] = Cuboid(myvec3(x1, (y1 + y2) / 2, (z1 + z2) / 2),myvec3((x1 + x2) / 2, y2, z2));
        subcuboids[7] = Cuboid((min_extent+max_extent)/2.0, max_extent);
        return subcuboids;
    }
    std::string print() const
    {
        std::ostringstream str;
        str << min_extent<< " "<<max_extent;
        std::string s = str.str();
        return s;
    }
};

#endif

// (Particle P1,Particle P2)
//     {
//         x1 = P1.x;
//         y1 = P1.y;
//         z1 = P1.z;
//         x2 = P2.x;
//         y2 = P2.y;
//         z2 = P2.z;
//     }