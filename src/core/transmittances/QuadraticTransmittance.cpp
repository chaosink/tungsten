#include "QuadraticTransmittance.hpp"

#include "sampling/UniformPathSampler.hpp"

#include "math/MathUtil.hpp"

#include "io/JsonObject.hpp"

#include "Memory.hpp"

namespace Tungsten {

QuadraticTransmittance::QuadraticTransmittance()
: _maxT(0.75f)
{
}

void QuadraticTransmittance::fromJson(JsonPtr value, const Scene &scene)
{
    Transmittance::fromJson(value, scene);
    value.getField("max_t", _maxT);
}

rapidjson::Value QuadraticTransmittance::toJson(Allocator &allocator) const
{
    return JsonObject{Transmittance::toJson(allocator), allocator,
        "type", "quadratic",
        "max_t", _maxT
    };
}

Vec3f QuadraticTransmittance::surfaceSurface(const Vec3f &tau) const
{
    Vec3f t = min(tau/_maxT, Vec3f(1.0f));
    return Float(1.0f) - t*2.0f + t*t;
}
Vec3f QuadraticTransmittance::surfaceMedium(const Vec3f &tau) const
{
    return (2.0f/_maxT)*(Float(1.0f) - min(tau/_maxT, Vec3f(1.0f)));
}
Vec3f QuadraticTransmittance::mediumSurface(const Vec3f &tau) const
{
    return Float(1.0f) - min(tau/_maxT, Vec3f(1.0f));
}
Vec3f QuadraticTransmittance::mediumMedium(const Vec3f &tau) const
{
    Vec3f result(1.0f/_maxT);
    for (int i = 0; i < 3; ++i)
        if (tau[i] > _maxT)
            result[i] = 0.0f;
    return result;
}

Float QuadraticTransmittance::sigmaBar() const
{
    return 2.0f/_maxT;
}

Float QuadraticTransmittance::sampleSurface(PathSampleGenerator &sampler) const
{
    return _maxT*(1.0f - std::sqrt(1.0f - sampler.next1D()));
}
Float QuadraticTransmittance::sampleMedium(PathSampleGenerator &sampler) const
{
    return _maxT*sampler.next1D();
}

}
