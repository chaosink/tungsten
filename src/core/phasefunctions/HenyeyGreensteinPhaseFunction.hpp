#ifndef HENYEYGREENSTEINPHASEFUNCTION_HPP_
#define HENYEYGREENSTEINPHASEFUNCTION_HPP_

#include "PhaseFunction.hpp"

namespace Tungsten {

class HenyeyGreensteinPhaseFunction : public PhaseFunction
{
    Float _g;

    inline Float henyeyGreenstein(Float cosTheta) const;

public:
    HenyeyGreensteinPhaseFunction();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    virtual Vec3f eval(const Vec3f &wi, const Vec3f &wo) const override;
    virtual bool sample(PathSampleGenerator &sampler, const Vec3f &wi, PhaseSample &sample) const override;
    virtual bool invert(WritablePathSampleGenerator &sampler, const Vec3f &wi, const Vec3f &wo) const;
    virtual Float pdf(const Vec3f &wi, const Vec3f &wo) const override;

    Float g() const
    {
        return _g;
    }
};

}

#endif /* HENYEYGREENSTEINPHASEFUNCTION_HPP_ */
