#include "OrenNayarBsdf.hpp"

#include "samplerecords/SurfaceScatterEvent.hpp"

#include "textures/ConstantTexture.hpp"

#include "sampling/PathSampleGenerator.hpp"
#include "sampling/SampleWarp.hpp"

#include "math/Angle.hpp"
#include "math/Vec.hpp"

#include "io/JsonObject.hpp"
#include "io/Scene.hpp"

namespace Tungsten {

OrenNayarBsdf::OrenNayarBsdf()
: _roughness(std::make_shared<ConstantTexture>(1.0f))
{
    _lobes = BsdfLobes(BsdfLobes::DiffuseReflectionLobe);
}

void OrenNayarBsdf::fromJson(JsonPtr value, const Scene &scene)
{
    Bsdf::fromJson(value, scene);

    if (auto roughness = value["roughness"])
        _roughness = scene.fetchTexture(roughness, TexelConversion::REQUEST_AVERAGE);
}

rapidjson::Value OrenNayarBsdf::toJson(Allocator &allocator) const
{
    return JsonObject{Bsdf::toJson(allocator), allocator,
        "type", "oren_nayar",
        "roughness", *_roughness
    };
}


bool OrenNayarBsdf::sample(SurfaceScatterEvent &event) const
{
    if (!event.requestedLobe.test(BsdfLobes::DiffuseReflectionLobe))
        return false;
    if (event.wi.z() <= 0.0f)
        return false;

    Float roughness = (*_roughness)[*event.info].x();
    Float ratio = clamp(roughness, Float(0.01f), Float(1.0f));
    if (event.sampler->nextBoolean(ratio))
        event.wo  = SampleWarp::uniformHemisphere(event.sampler->next2D());
    else
        event.wo  = SampleWarp::cosineHemisphere(event.sampler->next2D());

    event.pdf = SampleWarp::uniformHemispherePdf(event.wo)*ratio + SampleWarp::cosineHemispherePdf(event.wo)*(1.0f - ratio);
    event.weight = eval(event)/event.pdf;
    event.sampledLobe = BsdfLobes::DiffuseReflectionLobe;
    return event.wo.z() > 0.0f;
}

Vec3f OrenNayarBsdf::eval(const SurfaceScatterEvent &event) const
{
    if (!event.requestedLobe.test(BsdfLobes::DiffuseReflectionLobe))
        return Vec3f(0.0f);
    if (event.wi.z() <= 0.0f || event.wo.z() <= 0.0f)
        return Vec3f(0.0f);

    const Vec3f &wi = event.wi;
    const Vec3f &wo = event.wo;

    Float thetaR = std::acos(event.wo.z());
    Float thetaI = std::acos(event.wi.z());
    Float alpha = max(thetaR, thetaI);
    Float beta  = min(thetaR, thetaI);
    Float sinAlpha = std::sin(alpha);
    Float denom = (wi.x()*wi.x() + wi.y()*wi.y())*(wo.x()*wo.x() + wo.y()*wo.y());
    Float cosDeltaPhi;
    if (denom == 0.0f)
        cosDeltaPhi = 1.0f;
    else
        cosDeltaPhi = (wi.x()*wo.x() + wi.y()*wo.y())/std::sqrt(denom);

    const Float RoughnessToSigma = 1.0f/std::sqrt(2.0f);
    Float sigma = RoughnessToSigma*(*_roughness)[*event.info].x();
    Float sigmaSq = sigma*sigma;

    Float C1 = 1.0f - 0.5f*sigmaSq/(sigmaSq + 0.33f);
    Float C2 = 0.45f*sigmaSq/(sigmaSq + 0.09f);
    if (cosDeltaPhi >= 0.0f)
        C2 *= sinAlpha;
    else
        C2 *= sinAlpha - cube((2.0f*INV_PI)*beta);
    Float C3 = 0.125f*(sigmaSq/(sigmaSq + 0.09f))*sqr((4.0f*INV_PI*INV_PI)*alpha*beta);

    Float fr1 = (C1 + cosDeltaPhi*C2*std::tan(beta) + (1.0f - std::abs(cosDeltaPhi))*C3*std::tan(0.5f*(alpha + beta)));
    Float fr2 = 0.17f*sigmaSq/(sigmaSq + 0.13f)*(1.0f - cosDeltaPhi*sqr((2.0f*INV_PI)*beta));

    Vec3f diffuseAlbedo = albedo(event.info);
    return (diffuseAlbedo*fr1 + diffuseAlbedo*diffuseAlbedo*fr2)*wo.z()*INV_PI;
}

bool OrenNayarBsdf::invert(WritablePathSampleGenerator &sampler, const SurfaceScatterEvent &event) const
{
    if (!event.requestedLobe.test(BsdfLobes::DiffuseReflectionLobe))
        return 0.0f;
    if (event.wi.z() <= 0.0f || event.wo.z() <= 0.0f)
        return 0.0f;

    Float roughness = (*_roughness)[*event.info].x();
    Float ratio = clamp(roughness, Float(0.01f), Float(1.0f));

    Float pdf0 = SampleWarp::uniformHemispherePdf(event.wo)*ratio;
    Float pdf1 = SampleWarp::cosineHemispherePdf(event.wo)*(1.0f - ratio);

    if (sampler.untrackedBoolean(pdf0/(pdf0 + pdf1))) {
        sampler.putBoolean(ratio, true);
        sampler.put2D(SampleWarp::invertUniformHemisphere(event.wo, sampler.untracked1D()));
    } else {
        sampler.putBoolean(ratio, false);
        sampler.put2D(SampleWarp::invertCosineHemisphere(event.wo, sampler.untracked1D()));
    }
    return true;
}

Float OrenNayarBsdf::pdf(const SurfaceScatterEvent &event) const
{
    if (!event.requestedLobe.test(BsdfLobes::DiffuseReflectionLobe))
        return 0.0f;
    if (event.wi.z() <= 0.0f || event.wo.z() <= 0.0f)
        return 0.0f;

    Float roughness = (*_roughness)[*event.info].x();
    Float ratio = clamp(roughness, Float(0.01f), Float(1.0f));
    return SampleWarp::uniformHemispherePdf(event.wo)*ratio + SampleWarp::cosineHemispherePdf(event.wo)*(1.0f - ratio);
}

}
