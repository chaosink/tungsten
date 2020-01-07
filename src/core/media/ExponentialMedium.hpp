#ifndef EXPONENTIALMEDIUM_HPP_
#define EXPONENTIALMEDIUM_HPP_

#include "Medium.hpp"

namespace Tungsten {

class ExponentialMedium : public Medium
{
    Vec3f _materialSigmaA, _materialSigmaS;
    Float _density;
    Float _falloffScale;
    Vec3f _unitPoint;
    Vec3f _falloffDirection;

    Vec3f _unitFalloffDirection;
    Vec3f _sigmaA, _sigmaS;
    Vec3f _sigmaT;
    bool _absorptionOnly;

    inline Float density(Vec3f p) const;
    inline Float density(Float x, Float dx, Float t) const;
    inline Float densityIntegral(Float x, Float dx, Float tMax) const;
    inline Float inverseOpticalDepth(Float x, Float dx, Float tau) const;

public:
    ExponentialMedium();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    virtual bool isHomogeneous() const override;

    virtual void prepareForRender() override;

    virtual Vec3f sigmaA(Vec3f p) const override;
    virtual Vec3f sigmaS(Vec3f p) const override;
    virtual Vec3f sigmaT(Vec3f p) const override;

    virtual bool sampleDistance(PathSampleGenerator &sampler, const Ray &ray,
            MediumState &state, MediumSample &sample) const override;
    virtual Vec3f transmittance(PathSampleGenerator &sampler, const Ray &ray, bool startOnSurface,
            bool endOnSurface) const override;
    virtual Float pdf(PathSampleGenerator &sampler, const Ray &ray, bool startOnSurface, bool endOnSurface) const override;
};

}



#endif /* EXPONENTIALMEDIUM_HPP_ */
