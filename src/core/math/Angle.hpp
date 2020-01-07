#ifndef ANGLE_HPP_
#define ANGLE_HPP_

#include "IntTypes.hpp"

#include <cmath>

namespace Tungsten {

CONSTEXPR Float PI          = 3.1415926536f;
CONSTEXPR Float PI_HALF     = PI*0.5f;
CONSTEXPR Float TWO_PI      = PI*2.0f;
CONSTEXPR Float FOUR_PI     = PI*4.0f;
CONSTEXPR Float INV_PI      = 1.0f/PI;
CONSTEXPR Float INV_TWO_PI  = 0.5f*INV_PI;
CONSTEXPR Float INV_FOUR_PI = 0.25f*INV_PI;
CONSTEXPR Float SQRT_PI     = 1.77245385091f;
CONSTEXPR Float INV_SQRT_PI = 1.0f/SQRT_PI;

class Angle
{
public:
    static CONSTEXPR Float radToDeg(Float a)
    {
        return a*(180.0f/PI);
    }

    static CONSTEXPR Float degToRad(Float a)
    {
        return a*(PI/180.0f);
    }

    static Float angleRepeat(Float a)
    {
        return std::fmod(std::fmod(a, TWO_PI) + TWO_PI, TWO_PI);
    }
};

}

#endif /* ANGLE_HPP_ */
