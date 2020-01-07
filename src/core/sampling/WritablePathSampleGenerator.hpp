#ifndef WRITABLEPATHSAMPLEGENERATOR_HPP_
#define WRITABLEPATHSAMPLEGENERATOR_HPP_

#include "PathSampleGenerator.hpp"

namespace Tungsten {

class WritablePathSampleGenerator : public PathSampleGenerator
{
public:
    virtual ~WritablePathSampleGenerator() {}

    virtual void seek(int vertex) = 0;
    virtual void putBoolean(Float pTrue, bool choice) = 0;
    virtual void putDiscrete(int numChoices, int choice) = 0;
    virtual void put1D(Float value) = 0;
    virtual void put2D(Vec2f value) = 0;

    virtual Float untracked1D() = 0;
    inline Vec2f untracked2D()
    {
        Float a = untracked1D();
        Float b = untracked1D();
        return Vec2f(a, b);
    }
    inline bool untrackedBoolean(Float pTrue)
    {
        return untracked1D() < pTrue;
    }
    inline int untrackedDiscrete(int numChoices)
    {
        return int(untracked1D()*numChoices);
    }
};

}

#endif /* WRITABLEPATHSAMPLEGENERATOR_HPP_ */
