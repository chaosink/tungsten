#ifndef LIGHTSAMPLE_HPP_
#define LIGHTSAMPLE_HPP_

#include "math/Vec.hpp"

namespace Tungsten {

class Medium;

struct LightSample
{
    Vec3f d;
    Float dist;
    Float pdf;
    const Medium *medium;
};

}

#endif /* LIGHTSAMPLE_HPP_ */
