#ifndef DAVISWEINSTEINTRANSMITTANCE_HPP_
#define DAVISWEINSTEINTRANSMITTANCE_HPP_

#include "Transmittance.hpp"

namespace Tungsten {

class DavisWeinsteinTransmittance : public Transmittance
{
    Float _h;
    Float _c;

    Float computeAlpha(Float tau) const;

public:
    DavisWeinsteinTransmittance();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    virtual Vec3f surfaceSurface(const Vec3f &tau) const override final;
    virtual Vec3f surfaceMedium(const Vec3f &tau) const override final;
    virtual Vec3f mediumSurface(const Vec3f &tau) const override final;
    virtual Vec3f mediumMedium(const Vec3f &tau) const override final;

    virtual Float sigmaBar() const override final;

    virtual Float sampleSurface(PathSampleGenerator &sampler) const override final;
    virtual Float sampleMedium(PathSampleGenerator &sampler) const override final;
};

}

#endif /* DAVISWEINSTEINTRANSMITTANCE_HPP_ */
