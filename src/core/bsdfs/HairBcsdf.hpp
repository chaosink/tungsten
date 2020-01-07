#ifndef DEONHAIR_HPP_
#define DEONHAIR_HPP_

#include "PrecomputedAzimuthalLobe.hpp"

#include "bsdfs/Bsdf.hpp"

#include "math/Angle.hpp"

#include <memory>

namespace Tungsten {

// An implementation of the papers "An Energy-Conserving Hair Reflectance Model"
// and "Importance Sampling for Physically-Based Hair Fiber Models"
// using precomputed azimuthal scattering functions
class HairBcsdf : public Bsdf
{
    const Float Eta = 1.55f;

    Float _scaleAngleDeg;
    Float _melaninRatio;
    Float _melaninConcentration;
    bool _overridesSigmaA;
    Vec3f _sigmaA;
    Float _roughness;

    Float _scaleAngleRad;
    std::unique_ptr<PrecomputedAzimuthalLobe> _nR, _nTT, _nTRT;
    Float _betaR, _betaTT, _betaTRT;
    Float _vR, _vTT, _vTRT;

    static Float I0(Float x);
    static Float logI0(Float x);

    static Float g(Float beta, Float theta);
    static Float D(Float beta, Float phi);

    static Float Phi(Float gammaI, Float gammaT, int p);

    Float M(Float v, Float sinThetaI, Float sinThetaO, Float cosThetaI, Float cosThetaO) const;

    Float NrIntegrand(Float beta, Float wiDotWo, Float phi, Float h) const;
    Vec3f NpIntegrand(Float beta, Float cosThetaD, Float phi, int p, Float h) const;

    Float sampleM(Float v, Float sinThetaI, Float cosThetaI, Float xi1, Float xi2) const;

    void precomputeAzimuthalDistributions();

public:
    HairBcsdf();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    rapidjson::Value toJson(Allocator &allocator) const override;

    virtual Vec3f eval(const SurfaceScatterEvent &event) const override;
    virtual bool sample(SurfaceScatterEvent &event) const override;
    virtual Float pdf(const SurfaceScatterEvent &event) const override;

    virtual void prepareForRender() override;
};

}

#endif /* DEONHAIR_HPP_ */
