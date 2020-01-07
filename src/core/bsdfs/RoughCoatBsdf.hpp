#ifndef ROUGHCOATBSDF_HPP_
#define ROUGHCOATBSDF_HPP_

#include "Bsdf.hpp"
#include "Microfacet.hpp"

namespace Tungsten {

class RoughCoatBsdf : public Bsdf
{
    Float _ior;
    Float _thickness;
    Vec3f _sigmaA;
    std::shared_ptr<Bsdf> _substrate;
    Microfacet::Distribution _distribution;
    std::shared_ptr<Texture> _roughness;

    Float _avgTransmittance;
    Vec3f _scaledSigmaA;

public:
    RoughCoatBsdf();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    void substrateEvalAndPdf(const SurfaceScatterEvent &event, Float eta,
            Float Fi, Float cosThetaTi, Float &pdf, Vec3f &brdf) const;

    //TODO transmissive substrate
    virtual bool sample(SurfaceScatterEvent &event) const override;
    virtual Vec3f eval(const SurfaceScatterEvent &event) const override;
    virtual bool invert(WritablePathSampleGenerator &sampler, const SurfaceScatterEvent &event) const override;
    virtual Float pdf(const SurfaceScatterEvent &event) const override;

    virtual void prepareForRender() override;

    const char *distributionName() const
    {
        return _distribution.toString();
    }

    Float ior() const
    {
        return _ior;
    }

    const std::shared_ptr<Texture> &roughness() const
    {
        return _roughness;
    }

    Vec3f sigmaA() const
    {
        return _sigmaA;
    }

    const std::shared_ptr<Bsdf> &substrate() const
    {
        return _substrate;
    }

    Float thickness() const
    {
        return _thickness;
    }

    void setDistributionName(const std::string &distributionName)
    {
        _distribution = distributionName;
    }

    void setIor(Float ior)
    {
        _ior = ior;
    }

    void setRoughness(const std::shared_ptr<Texture> &roughness)
    {
        _roughness = roughness;
    }

    void setSigmaA(Vec3f sigmaA)
    {
        _sigmaA = sigmaA;
    }

    void setSubstrate(const std::shared_ptr<Bsdf> &substrate)
    {
        _substrate = substrate;
    }

    void setThickness(Float thickness)
    {
        _thickness = thickness;
    }
};

}

#endif /* ROUGHCOATBSDF_HPP_ */
