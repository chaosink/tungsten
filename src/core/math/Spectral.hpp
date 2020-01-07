#ifndef SPECTRAL_HPP_
#define SPECTRAL_HPP_

#include "math/MathUtil.hpp"
#include "math/Vec.hpp"

namespace Tungsten {

namespace Spectral {

static CONSTEXPR int CIE_samples = 471;
static CONSTEXPR Float CIE_Min = 360.0f;
static CONSTEXPR Float CIE_Max = 830.0f;

extern const Float CIE_X_entries[];
extern const Float CIE_Y_entries[];
extern const Float CIE_Z_entries[];

void spectralXyzWeights(int samples, Float lambdas[], Vec3f weights[]);

inline Vec3f xyzToRgb(Vec3f xyz) {
    return Vec3f(
        3.240479f*xyz.x() + -1.537150f*xyz.y() + -0.498535f*xyz.z(),
       -0.969256f*xyz.x() +  1.875991f*xyz.y() +  0.041556f*xyz.z(),
        0.055648f*xyz.x() + -0.204043f*xyz.y() +  1.057311f*xyz.z()
    );
}

static inline Vec3f wavelengthToXyz(Float lambda)
{
    Float x = CIE_samples*(lambda - CIE_Min)/(CIE_Max - CIE_Min);
    int i = clamp(int(x), 0, CIE_samples - 2);
    Float u = x - i;
    Vec3f xyz0 = Vec3f(CIE_X_entries[i + 0], CIE_Y_entries[i + 0], CIE_Z_entries[i + 0]);
    Vec3f xyz1 = Vec3f(CIE_X_entries[i + 1], CIE_Y_entries[i + 1], CIE_Z_entries[i + 1]);

    return xyz0*(1.0f - u) + xyz1*u;
}

static inline Vec3f wavelengthToRgb(Float lambda)
{
    return xyzToRgb(wavelengthToXyz(lambda));
}

}

}

#endif /* SPECTRAL_HPP_ */
