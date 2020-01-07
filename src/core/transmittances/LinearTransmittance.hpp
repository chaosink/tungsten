#ifndef LINEARTRANSMITTANCE_HPP_
#define LINEARTRANSMITTANCE_HPP_

#include "Transmittance.hpp"

namespace Tungsten {

class LinearTransmittance : public Transmittance
{
    Float _maxT;

public:
    LinearTransmittance();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    virtual Vec3f surfaceSurface(const Vec3f &tau) const override final;
    virtual Vec3f surfaceMedium(const Vec3f &tau) const override final;
    virtual Vec3f mediumSurface(const Vec3f &tau) const override final;
    virtual Vec3f mediumMedium(const Vec3f &tau) const override final;

    virtual bool isDirac() const override final;

    virtual Float sigmaBar() const override final;

    virtual Float sampleSurface(PathSampleGenerator &sampler) const override final;
    virtual Float sampleMedium(PathSampleGenerator &sampler) const override final;
};

}

#endif /* LINEARTRANSMITTANCE_HPP_ */
