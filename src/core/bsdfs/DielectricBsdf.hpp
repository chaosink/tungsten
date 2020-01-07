#ifndef DIELECTRICBSDF_HPP_
#define DIELECTRICBSDF_HPP_

#include "Bsdf.hpp"

namespace Tungsten {

class Scene;

class DielectricBsdf : public Bsdf
{
    Float _ior;
    Float _invIor;
    bool _enableT;

public:
    DielectricBsdf();
    DielectricBsdf(Float ior);

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    virtual bool sample(SurfaceScatterEvent &event) const override;
    virtual Vec3f eval(const SurfaceScatterEvent &event) const override;
    virtual bool invert(WritablePathSampleGenerator &sampler, const SurfaceScatterEvent &event) const override;
    virtual Float pdf(const SurfaceScatterEvent &event) const override;
    virtual Float eta(const SurfaceScatterEvent &event) const override;

    virtual void prepareForRender() override;


    bool enableTransmission() const
    {
        return _enableT;
    }

    Float ior() const
    {
        return _ior;
    }

    void setEnableTransmission(bool enableTransmission)
    {
        _enableT = enableTransmission;
    }

    void setIor(Float ior)
    {
        _ior = ior;
    }
};

}


#endif /* DIELECTRICBSDF_HPP_ */
