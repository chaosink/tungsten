#ifndef RAY_HPP_
#define RAY_HPP_

#include "MathUtil.hpp"
#include "Angle.hpp"
#include "Vec.hpp"

#include <limits>

namespace Tungsten {

class Ray
{
    Vec3f _pos;
    Vec3f _dir;
    Float _nearT;
    Float _farT;
    Float _time;
    bool _primaryRay;

public:
    Ray() = default;

    Ray(const Vec3f &pos, const Vec3f &dir, Float nearT = 1e-4f, Float farT = infinity(), Float time = 0.0f)
    : _pos(pos), _dir(dir), _nearT(nearT), _farT(farT), _time(time), _primaryRay(false)
    {
    }

    Ray scatter(const Vec3f &newPos, const Vec3f &newDir, Float newNearT, Float newFarT = infinity()) const
    {
        Ray ray(*this);
        ray._pos = newPos;
        ray._dir = newDir;
        ray._nearT = newNearT;
        ray._farT = newFarT;
        return ray;
    }

    Vec3f hitpoint() const
    {
        return _pos + _dir*_farT;
    }

    const Vec3f& dir() const
    {
        return _dir;
    }

    void setDir(const Vec3f& dir)
    {
        _dir = dir;
    }

    const Vec3f& pos() const
    {
        return _pos;
    }

    void setPos(const Vec3f& pos)
    {
        _pos = pos;
    }

    Float farT() const
    {
        return _farT;
    }

    void setFarT(Float farT)
    {
        _farT = farT;
    }

    Float nearT() const
    {
        return _nearT;
    }

    void setNearT(Float nearT)
    {
        _nearT = nearT;
    }

    Float time() const
    {
        return _time;
    }

    void setTime(Float time)
    {
        _time = time;
    }

    bool isPrimaryRay() const
    {
        return _primaryRay;
    }

    void setPrimaryRay(bool value)
    {
        _primaryRay = value;
    }

    static inline Float infinity()
    {
        return std::numeric_limits<Float>::infinity();
    }
};

}

#endif /* RAY_HPP_ */
