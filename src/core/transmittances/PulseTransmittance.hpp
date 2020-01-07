#ifndef PULSETRANSMITTANCE_HPP_
#define PULSETRANSMITTANCE_HPP_

#include "Transmittance.hpp"

namespace Tungsten {

class PulseTransmittance : public Transmittance
{
    Float _a, _b;
    int _numPulses;

public:
    PulseTransmittance();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    virtual bool isDirac() const override final;

    virtual Vec3f surfaceSurface(const Vec3f &tau) const override final;
    virtual Vec3f surfaceMedium(const Vec3f &tau) const override final;
    virtual Vec3f mediumSurface(const Vec3f &tau) const override final;
    virtual Vec3f mediumMedium(const Vec3f &tau) const override final;

    virtual Float sigmaBar() const override final;

    virtual Float sampleSurface(PathSampleGenerator &sampler) const override final;
    virtual Float sampleMedium(PathSampleGenerator &sampler) const override final;
};


}



#endif /* PULSETRANSMITTANCE_HPP_ */
