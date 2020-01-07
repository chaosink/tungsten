#ifndef MICROFACET_HPP_
#define MICROFACET_HPP_

#include "sampling/SampleWarp.hpp"

#include "math/MathUtil.hpp"
#include "math/Angle.hpp"
#include "math/Vec.hpp"

#include "StringableEnum.hpp"

namespace Tungsten {

class Microfacet
{
    enum DistributionEnum
    {
        Beckmann,
        Phong,
        GGX
    };

public:
    typedef StringableEnum<DistributionEnum> Distribution;
    friend Distribution;

    static Float roughnessToAlpha(DistributionEnum dist, Float roughness)
    {
        CONSTEXPR Float MinAlpha = 1e-3f;
        roughness = max(roughness, MinAlpha);

        if (dist == Phong)
            return 2.0f/(roughness*roughness) - 2.0f;
        else
            return roughness;
    }

    static Float D(DistributionEnum dist, Float alpha, const Vec3f &m)
    {
        if (m.z() <= 0.0f)
            return 0.0f;

        switch (dist) {
        case Beckmann: {
            Float alphaSq = alpha*alpha;
            Float cosThetaSq = m.z()*m.z();
            Float tanThetaSq = max(1.0f - cosThetaSq, Float(0.0f))/cosThetaSq;
            Float cosThetaQu = cosThetaSq*cosThetaSq;
            return INV_PI*std::exp(-tanThetaSq/alphaSq)/(alphaSq*cosThetaQu);
        }
        case Phong:
            return (alpha + 2.0f)*INV_TWO_PI*Float(std::pow(double(m.z()), double(alpha)));
        case GGX: {
            Float alphaSq = alpha*alpha;
            Float cosThetaSq = m.z()*m.z();
            Float tanThetaSq = max(1.0f - cosThetaSq, Float(0.0f))/cosThetaSq;
            Float cosThetaQu = cosThetaSq*cosThetaSq;
            return alphaSq*INV_PI/(cosThetaQu*sqr(alphaSq + tanThetaSq));
        }
        }

        return 0.0f;
    }

    static Float G1(DistributionEnum dist, Float alpha, const Vec3f &v, const Vec3f &m)
    {
        if (v.dot(m)*v.z() <= 0.0f)
            return 0.0f;

        switch (dist) {
        case Beckmann: {
            Float cosThetaSq = v.z()*v.z();
            Float tanTheta = std::abs(std::sqrt(max(1.0f - cosThetaSq, Float(0.0f)))/v.z());
            Float a = 1.0f/(alpha*tanTheta);
            if (a < 1.6f)
                return (3.535f*a + 2.181f*a*a)/(1.0f + 2.276f*a + 2.577f*a*a);
            else
                return 1.0f;
        } case Phong: {
            Float cosThetaSq = v.z()*v.z();
            Float tanTheta = std::abs(std::sqrt(max(1.0f - cosThetaSq, Float(0.0f)))/v.z());
            Float a = std::sqrt(0.5f*alpha + 1.0f)/tanTheta;
            if (a < 1.6f)
                return (3.535f*a + 2.181f*a*a)/(1.0f + 2.276f*a + 2.577f*a*a);
            else
                return 1.0f;
        } case GGX: {
            Float alphaSq = alpha*alpha;
            Float cosThetaSq = v.z()*v.z();
            Float tanThetaSq = max(1.0f - cosThetaSq, Float(0.0f))/cosThetaSq;
            return 2.0f/(1.0f + std::sqrt(1.0f + alphaSq*tanThetaSq));
        }
        }

        return 0.0f;
    }

    static Float G(DistributionEnum dist, Float alpha, const Vec3f &i, const Vec3f &o, const Vec3f &m)
    {
        return G1(dist, alpha, i, m)*G1(dist, alpha, o, m);
    }

    static Float pdf(DistributionEnum dist, Float alpha, const Vec3f &m)
    {
        return D(dist, alpha, m)*m.z();
    }

    static Vec3f sample(DistributionEnum dist, Float alpha, Vec2f xi)
    {
        Float phi = xi.y()*TWO_PI;
        Float cosTheta = 0.0f;

        switch (dist) {
        case Beckmann: {
            Float tanThetaSq = -alpha*alpha*std::log(1.0f - xi.x());
            cosTheta = 1.0f/std::sqrt(1.0f + tanThetaSq);
            break;
        } case Phong:
            cosTheta = Float(std::pow(double(xi.x()), 1.0/(double(alpha) + 2.0)));
            break;
        case GGX: {
            Float tanThetaSq = alpha*alpha*xi.x()/(1.0f - xi.x());
            cosTheta = 1.0f/std::sqrt(1.0f + tanThetaSq);
            break;
        }
        }

        Float r = std::sqrt(max(1.0f - cosTheta*cosTheta, Float(0.0f)));
        return Vec3f(std::cos(phi)*r, std::sin(phi)*r, cosTheta);
    }

    static inline Vec2f invert(Distribution dist, Float alpha, const Vec3d &m, Float mu)
    {
        Vec2d xi;
        xi.y() = SampleWarp::invertPhi(m, mu);

        double cosTheta = m.z();
        switch (dist) {
        case Beckmann: {
            double tanThetaSq = 1.0/(cosTheta*cosTheta) - 1.0;
            xi.x() = 1.0 - std::exp(-tanThetaSq/(alpha*alpha));
            break;
        } case Phong: {
            xi.x() = std::pow(cosTheta, double(alpha) + 2.0);
            break;
        } case GGX: {
            double tanThetaSq = 1.0/(cosTheta*cosTheta) - 1.0;
            double gamma = tanThetaSq/(alpha*alpha);
            xi.x() = gamma/(1.0 + gamma);
            break;
        } default:
            xi.x() = 0.0;
        }

        return Vec2f(xi);
    }
};

}

#endif /* MICROFACET_HPP_ */
