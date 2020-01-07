#ifndef MEDIUMSAMPLE_HPP_
#define MEDIUMSAMPLE_HPP_

#include "math/Vec.hpp"

namespace Tungsten {

class PhaseFunction;

struct MediumSample
{
    PhaseFunction *phase;
    Vec3f p;
    Float continuedT;
    Vec3f continuedWeight;
    Float t;
    Vec3f weight;
    Vec3f emission;
    Float pdf;
    bool exited;
};

}

#endif /* MEDIUMSAMPLE_HPP_ */
