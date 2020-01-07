#ifndef PHONGBSDF_HPP_
#define PHONGBSDF_HPP_

#include "Bsdf.hpp"

namespace Tungsten {

class Scene;

class PhongBsdf : public Bsdf
{
    Float _exponent;
    Float _invExponent;
    Float _pdfFactor;
    Float _brdfFactor;
    Float _diffuseRatio;

public:
    PhongBsdf(Float exponent = 64.0f, Float diffuseRatio = 0.2f);

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    virtual bool sample(SurfaceScatterEvent &event) const override;
    virtual Vec3f eval(const SurfaceScatterEvent &event) const override;
    virtual Float pdf(const SurfaceScatterEvent &event) const override;

    virtual void prepareForRender() override;

    Float exponent() const
    {
        return _exponent;
    }

    Float diffuseRatio() const
    {
        return _diffuseRatio;
    }

    void setDiffuseRatio(Float diffuseRatio)
    {
        _diffuseRatio = diffuseRatio;
    }

    void setExponent(Float exponent)
    {
        _exponent = exponent;
    }
};

}


#endif /* PHONGBSDF_HPP_ */
