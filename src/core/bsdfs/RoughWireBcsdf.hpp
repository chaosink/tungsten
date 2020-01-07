#ifndef ROUGHWIREBCSDF_HPP_
#define ROUGHWIREBCSDF_HPP_

#include "bsdfs/Bsdf.hpp"

namespace Tungsten {

class RoughWireBcsdf : public Bsdf
{
    std::string _materialName;
    Vec3f _eta;
    Vec3f _k;
    Float _roughness;

    Float _v;

    bool lookupMaterial();

    static Float I0(Float x);
    static Float logI0(Float x);

    Float N(Float cosPhi) const;
    Float M(Float v, Float sinThetaI, Float sinThetaO, Float cosThetaI, Float cosThetaO) const;

    Float sampleN(Float xi) const;
    Float sampleM(Float v, Float sinThetaI, Float cosThetaI, Float xi1, Float xi2) const;

public:
    RoughWireBcsdf();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    rapidjson::Value toJson(Allocator &allocator) const override;

    virtual Vec3f eval(const SurfaceScatterEvent &event) const override;
    virtual bool sample(SurfaceScatterEvent &event) const override;
    virtual Float pdf(const SurfaceScatterEvent &event) const override;

    virtual void prepareForRender() override;
};

}

#endif /* ROUGHWIREBCSDF_HPP_ */
