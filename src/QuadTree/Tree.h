#ifndef TREE
#define TREE

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include "Rectangle.h"
#include "Cuboid.h"
#include "Particle.h"
class Tree
{
    int maxPnts = 1;
    Particle COM;

public:
    std::vector<Tree> branches = {};
    //  p = std::make_unique<D>();
    std::vector<std::shared_ptr<Particle>> particles = {};
    int level = 0;
    Cuboid cuboid;
    bool leaf = true;
    Tree(Cuboid cuboidIn)
    {
        cuboid = cuboidIn;
    }

    void createBranches()
    {
        for (Cuboid r : cuboid.subdivide())
        {
            Tree branch = Tree(r);
            branch.level = level + 1;
            branches.push_back(branch);
        }
    }

    void insertPnt(const Particle &p) // make this boolean so that if you do insert it to stop inserting to other branches
    {
        if (!cuboid.contains(p))
            return;
        particles.push_back(std::make_shared<Particle>(p));
        if (particles.size() > maxPnts)
        {
            if (leaf)
            {
                leaf = false;
                createBranches();
                for (Tree &b : branches)
                {
                    for (auto &pnt : particles)
                    {
                        b.insertPnt(*pnt);
                    }
                }
            }
            else
            {
                for (Tree &b : branches)
                {
                    b.insertPnt(p);
                }
            }
        }
    }

    Particle computeCOM()
    {
        Particle tempCOM;
        if (leaf)
        {
            if (particles.size() != 0)
            {
                for (const std::shared_ptr<Particle> P : particles)
                {
                    Particle tempP = *P;
                    tempP.p *= tempP.m;
                    tempCOM = tempCOM + tempP;
                }
                if(tempCOM.m != 0)
                    tempCOM.p /= tempCOM.m;
            }
        }
        else
        {
            for (Tree t : branches)
            {
                Particle tempP = t.computeCOM();
                tempP.p *= tempP.m;
                tempCOM = tempCOM + tempP;
            }
            if(tempCOM.m != 0)
                tempCOM.p /= tempCOM.m;
        }
        COM = tempCOM;
        std::cout << cuboid.print() << " level: " << level << " nPoints:" << particles.size() << " Center of Mass: " << COM << ", leaf:" << leaf << std::endl;
        return tempCOM;
    }

    myvec3 computeAcc(const Particle &particle, myfloat theta) const
    {
        myfloat softening_param = 0.025;
        myvec3 acc{0,0,0};

        if (leaf)
        {
            for (auto sp : particles)
            {
                if ((*sp).p != particle.p)
                {
                    auto diff = (*sp).p - particle.p;
                    auto distance = glm::length(diff);
                    acc += diff * (*sp).m / (distance * distance * distance + softening_param);
                }
            }
        }
        else
        {
            if(less_than_theta(particle, theta)){
                auto diff = COM.p - particle.p;
                auto distance = glm::length(diff);
                acc = diff * COM.m / (distance * distance * distance + softening_param);
            }
            else{
                for (const Tree &b : branches)
                {
                    acc += b.computeAcc(particle, theta);
                }
            }
        }

        return acc;

        // }
    }

    bool less_than_theta(const Particle &particle, double theta) const {
        myfloat distance = glm::length(particle.p - COM.p);
        //TODO: this should not be the length of the cuboid
        return glm::length(cuboid.max_extent - cuboid.min_extent)/distance < theta;
    }

    std::string print() const
    {
        std::ostringstream str;
        str << cuboid.print() << " level: " << level << " nPoints:" << particles.size() << " Center of Mass: " << COM << ", leaf:" << leaf << std::endl;
        for (const Tree b : branches)
        {
            str << b.print();
        }
        return str.str();
    }
};
#endif
