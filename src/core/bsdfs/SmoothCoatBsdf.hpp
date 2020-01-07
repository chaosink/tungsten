#ifndef SMOOTHCOATBSDF_HPP_
#define SMOOTHCOATBSDF_HPP_

#include "Bsdf.hpp"

namespace Tungsten {

class SmoothCoatBsdf : public Bsdf
{
    Float _ior;
    Float _thickness;
    Vec3f _sigmaA;
    std::shared_ptr<Bsdf> _substrate;

    Float _avgTransmittance;
    Vec3f _scaledSigmaA;

public:
    SmoothCoatBsdf();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    //TODO transmissive substrate
    virtual bool sample(SurfaceScatterEvent &event) const override;
    virtual bool invert(WritablePathSampleGenerator &sampler, const SurfaceScatterEvent &event) const override;
    virtual Vec3f eval(const SurfaceScatterEvent &event) const override;
    virtual Float pdf(const SurfaceScatterEvent &event) const override;

    virtual void prepareForRender() override;

    Float ior() const
    {
        return _ior;
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

    void setIor(Float ior)
    {
        _ior = ior;
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

#endif /* SMOOTHCOATBSDF_HPP_ */
