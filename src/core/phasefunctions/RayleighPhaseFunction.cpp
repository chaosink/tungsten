#include "RayleighPhaseFunction.hpp"

#include "sampling/PathSampleGenerator.hpp"
#include "sampling/SampleWarp.hpp"

#include "math/TangentFrame.hpp"
#include "math/MathUtil.hpp"
#include "math/Angle.hpp"

#include "io/JsonObject.hpp"

namespace Tungsten {

inline Float RayleighPhaseFunction::rayleigh(Float cosTheta)
{
    return (3.0f/(16.0f*PI))*(1.0f + cosTheta*cosTheta);
}

rapidjson::Value RayleighPhaseFunction::toJson(Allocator &allocator) const
{
    return JsonObject{PhaseFunction::toJson(allocator), allocator,
        "type", "rayleigh"
    };
}

Vec3f RayleighPhaseFunction::eval(const Vec3f &wi, const Vec3f &wo) const
{
    return Vec3f(rayleigh(wi.dot(wo)));
}

bool RayleighPhaseFunction::sample(PathSampleGenerator &sampler, const Vec3f &wi, PhaseSample &sample) const
{
    Vec2f xi = sampler.next2D();
    Float phi = xi.x()*TWO_PI;
    Float z = xi.y()*4.0f - 2.0f;
    Float invZ = std::sqrt(z*z + 1.0f);
    Float u = std::cbrt(z + invZ);
    Float cosTheta = u - 1.0f/u;

    Float sinTheta = std::sqrt(max(1.0f - cosTheta*cosTheta, Float(0.0f)));
    sample.w = TangentFrame(wi).toGlobal(Vec3f(
        std::cos(phi)*sinTheta,
        std::sin(phi)*sinTheta,
        cosTheta
    ));
    sample.weight = Vec3f(1.0f);
    sample.pdf = rayleigh(cosTheta);
    return true;
}

bool RayleighPhaseFunction::invert(WritablePathSampleGenerator &sampler, const Vec3f &wi, const Vec3f &wo) const
{
    Vec3f w = TangentFrame(wi).toLocal(wo);
    Float cosTheta = w.z();
    Float u = 0.5f*(cosTheta + std::sqrt(4.0f + cosTheta*cosTheta));
    Float u3 = u*u*u;
    Float z = (u3*u3 - 1.0f)/(2.0f*u3);
    Float xi2 = (z + 2.0f)*0.25f;
    Float xi1 = SampleWarp::invertPhi(w, sampler.untracked1D());

    sampler.put2D(Vec2f(xi1, xi2));

    return true;
}

Float RayleighPhaseFunction::pdf(const Vec3f &wi, const Vec3f &wo) const
{
    return rayleigh(wi.dot(wo));
}

}
