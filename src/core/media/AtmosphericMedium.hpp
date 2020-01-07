#ifndef ATMOSPHERICMEDIUM_HPP_
#define ATMOSPHERICMEDIUM_HPP_

#include "Medium.hpp"

#include "textures/Texture.hpp"

#include <memory>

namespace Tungsten {

class AtmosphericMedium : public Medium
{
    const Scene *_scene;
    std::string _primName;

    Vec3f _materialSigmaA, _materialSigmaS;
    Float _density;
    Float _falloffScale;
    Float _radius;
    Vec3f _center;

    Float _effectiveFalloffScale;
    Vec3f _sigmaA, _sigmaS;
    Vec3f _sigmaT;
    bool _absorptionOnly;

    inline Float density(Vec3f p) const;
    inline Float density(Float h, Float t0) const;
    inline Float densityIntegral(Float h, Float t0, Float t1) const;
    inline Float inverseOpticalDepth(double h, double t0, double tau) const;

public:
    AtmosphericMedium();

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

#endif /* ATMOSPHERICMEDIUM_HPP_ */
