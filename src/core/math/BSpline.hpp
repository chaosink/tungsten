#ifndef BSPLINE_HPP_
#define BSPLINE_HPP_

namespace Tungsten {

namespace BSpline {

// http://www.answers.com/topic/b-spline
template<typename T>
inline T quadratic(T p0, T p1, T p2, Float t)
{
    return (p0*0.5f - p1 + p2*0.5f)*t*t + (p1 - p0)*t + (p0 + p1)*0.5f;
}

template<typename T>
inline T quadraticDeriv(T p0, T p1, T p2, Float t)
{
    return (p0 - p1*2.0f + p2)*t + (p1 - p0);
}

inline Vec2f quadraticMinMax(Float p0, Float p1, Float p2)
{
    Float xMin = (p0 + p1)*0.5f;
    Float xMax = (p1 + p2)*0.5f;
    if (xMin > xMax)
        std::swap(xMin, xMax);

    Float tFlat = (p0 - p1)/(p0 - 2.0f*p1 + p2);
    if (tFlat > 0.0f && tFlat < 1.0f) {
        Float xFlat = quadratic(p0, p1, p2, tFlat);
        xMin = min(xMin, xFlat);
        xMax = max(xMax, xFlat);
    }
    return Vec2f(xMin, xMax);
}

}

}

#endif /* BSPLINE_HPP_ */
