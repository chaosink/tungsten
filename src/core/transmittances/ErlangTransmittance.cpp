#include "ErlangTransmittance.hpp"

#include "sampling/UniformPathSampler.hpp"

#include "math/MathUtil.hpp"

#include "io/JsonObject.hpp"

#include "Memory.hpp"

namespace Tungsten {

ErlangTransmittance::ErlangTransmittance()
: _lambda(5.0f)
{
}

void ErlangTransmittance::fromJson(JsonPtr value, const Scene &scene)
{
    Transmittance::fromJson(value, scene);
    value.getField("rate", _lambda);
}

rapidjson::Value ErlangTransmittance::toJson(Allocator &allocator) const
{
    return JsonObject{Transmittance::toJson(allocator), allocator,
        "type", "erlang",
        "rate", _lambda,
    };
}

Vec3f ErlangTransmittance::surfaceSurface(const Vec3f &tau) const
{
    return Float(0.5f)*std::exp(-_lambda*tau)*(Float(2.0f) + _lambda*tau);
}
Vec3f ErlangTransmittance::surfaceMedium(const Vec3f &tau) const
{
    return mediumSurface(tau)*_lambda*0.5f;
}
Vec3f ErlangTransmittance::mediumSurface(const Vec3f &tau) const
{
    return std::exp(-_lambda*tau)*(Float(1.0f) + _lambda*tau);
}
Vec3f ErlangTransmittance::mediumMedium(const Vec3f &tau) const
{
    return sqr(_lambda)*tau*std::exp(-_lambda*tau);
}

Float ErlangTransmittance::sigmaBar() const
{
    return _lambda*0.5f;
}

Float ErlangTransmittance::sampleSurface(PathSampleGenerator &sampler) const
{
    Float xi = sampler.next1D();
    Float x = 0.5f;
    for (int i = 0; i < 10; ++i) {
        x += (xi - (1.0f - surfaceSurface(Vec3f(x))[0]))/surfaceMedium(Vec3f(x))[0];
        x = max(x, Float(0.0f));
    }
    return x;
}
Float ErlangTransmittance::sampleMedium(PathSampleGenerator &sampler) const
{
    return -1.0f/_lambda*std::log(sampler.next1D()*sampler.next1D());
}

}
