#ifndef FRESNEL_HPP_
#define FRESNEL_HPP_

#include "math/MathUtil.hpp"
#include "math/Angle.hpp"

#include <cmath>

namespace Tungsten {

namespace Fresnel {

// Computes total reflectance from an infinitesimally thin film with refraction index eta
// from all internal reflection/refraction events
static inline Float thinFilmReflectance(Float eta, Float cosThetaI, Float &cosThetaT)
{
    Float sinThetaTSq = eta*eta*(1.0f - cosThetaI*cosThetaI);
    if (sinThetaTSq > 1.0f) {
        cosThetaT = 0.0f;
        return 1.0f;
    }
    cosThetaT = std::sqrt(max(1.0f - sinThetaTSq, Float(0.0f)));

    Float Rs = sqr((eta*cosThetaI - cosThetaT)/(eta*cosThetaI + cosThetaT));
    Float Rp = sqr((eta*cosThetaT - cosThetaI)/(eta*cosThetaT + cosThetaI));

    return 1.0f - ((1.0f - Rs)/(1.0f + Rs) + (1.0f - Rp)/(1.0f + Rp))*0.5f;
}

static inline Float thinFilmReflectance(Float eta, Float cosThetaI)
{
    Float cosThetaT;
    return thinFilmReflectance(eta, cosThetaI, cosThetaT);
}

// Computes total reflectance including spectral interference from a thin film
// that is `thickness' nanometers thick and has refraction index eta
// See http://www.gamedev.net/page/resources/_/technical/graphics-programming-and-theory/thin-film-interference-for-computer-graphics-r2962
static inline Vec3f thinFilmReflectanceInterference(Float eta, Float cosThetaI, Float thickness, Float &cosThetaT)
{
    const Vec3f invLambdas = Float(1.0f)/Vec3f(650.0f, 510.0f, 475.0f);

    Float cosThetaISq = cosThetaI*cosThetaI;
    Float sinThetaISq = 1.0f - cosThetaISq;
    Float invEta = 1.0f/eta;

    Float sinThetaTSq = eta*eta*sinThetaISq;
    if (sinThetaTSq > 1.0f) {
        cosThetaT = 0.0f;
        return Vec3f(1.0f);
    }
    cosThetaT = std::sqrt(1.0f - sinThetaTSq);

    Float Ts = 4.0f*eta*cosThetaI*cosThetaT/sqr(eta*cosThetaI + cosThetaT);
    Float Tp = 4.0f*eta*cosThetaI*cosThetaT/sqr(eta*cosThetaT + cosThetaI);

    Float Rs = 1.0f - Ts;
    Float Rp = 1.0f - Tp;

    Vec3f phi = (thickness*cosThetaT*FOUR_PI*invEta)*invLambdas;
    Vec3f cosPhi(std::cos(phi.x()), std::cos(phi.y()), std::cos(phi.z()));

    Vec3f tS = sqr(Ts)/((sqr(Rs) + 1.0f) - 2.0f*Rs*cosPhi);
    Vec3f tP = sqr(Tp)/((sqr(Rp) + 1.0f) - 2.0f*Rp*cosPhi);

    return Float(1.0f) - (tS + tP)*0.5f;
}

static inline Vec3f thinFilmReflectanceInterference(Float eta, Float cosThetaI, Float thickness)
{
    Float cosThetaT;
    return thinFilmReflectanceInterference(eta, cosThetaI, thickness, cosThetaT);
}

static inline Float dielectricReflectance(Float eta, Float cosThetaI, Float &cosThetaT)
{
    if (cosThetaI < 0.0f) {
        eta = 1.0f/eta;
        cosThetaI = -cosThetaI;
    }
    Float sinThetaTSq = eta*eta*(1.0f - cosThetaI*cosThetaI);
    if (sinThetaTSq > 1.0f) {
        cosThetaT = 0.0f;
        return 1.0f;
    }
    cosThetaT = std::sqrt(max(1.0f - sinThetaTSq, Float(0.0f)));

    Float Rs = (eta*cosThetaI - cosThetaT)/(eta*cosThetaI + cosThetaT);
    Float Rp = (eta*cosThetaT - cosThetaI)/(eta*cosThetaT + cosThetaI);

    return (Rs*Rs + Rp*Rp)*0.5f;
}

static inline Float dielectricReflectance(Float eta, Float cosThetaI)
{
    Float cosThetaT;
    return dielectricReflectance(eta, cosThetaI, cosThetaT);
}

// From "PHYSICALLY BASED LIGHTING CALCULATIONS FOR COMPUTER GRAPHICS" by Peter Shirley
// http://www.cs.virginia.edu/~jdl/bib/globillum/shirley_thesis.pdf
static inline Float conductorReflectance(Float eta, Float k, Float cosThetaI)
{
    Float cosThetaISq = cosThetaI*cosThetaI;
    Float sinThetaISq = max(1.0f - cosThetaISq, Float(0.0f));
    Float sinThetaIQu = sinThetaISq*sinThetaISq;

    Float innerTerm = eta*eta - k*k - sinThetaISq;
    Float aSqPlusBSq = std::sqrt(max(innerTerm*innerTerm + 4.0f*eta*eta*k*k, Float(0.0f)));
    Float a = std::sqrt(max((aSqPlusBSq + innerTerm)*0.5f, Float(0.0f)));

    Float Rs = ((aSqPlusBSq + cosThetaISq) - (2.0f*a*cosThetaI))/
               ((aSqPlusBSq + cosThetaISq) + (2.0f*a*cosThetaI));
    Float Rp = ((cosThetaISq*aSqPlusBSq + sinThetaIQu) - (2.0f*a*cosThetaI*sinThetaISq))/
               ((cosThetaISq*aSqPlusBSq + sinThetaIQu) + (2.0f*a*cosThetaI*sinThetaISq));

    return 0.5f*(Rs + Rs*Rp);
}

static inline Float conductorReflectanceApprox(Float eta, Float k, Float cosThetaI)
{
    Float cosThetaISq = cosThetaI*cosThetaI;
    Float ekSq = eta*eta* + k*k;
    Float cosThetaEta2 = cosThetaI*2.0f*eta;

    Float Rp = (ekSq*cosThetaISq - cosThetaEta2 + 1.0f)/(ekSq*cosThetaISq + cosThetaEta2 + 1.0f);
    Float Rs = (ekSq - cosThetaEta2 + cosThetaISq)/(ekSq + cosThetaEta2 + cosThetaISq);
    return (Rs + Rp)*0.5f;
}

static inline Vec3f conductorReflectance(const Vec3f &eta, const Vec3f &k, Float cosThetaI)
{
    return Vec3f(
        conductorReflectance(eta.x(), k.x(), cosThetaI),
        conductorReflectance(eta.y(), k.y(), cosThetaI),
        conductorReflectance(eta.z(), k.z(), cosThetaI)
    );
}

// Computes hemispherical integral of dielectricReflectance(ior, cos(theta))*cos(theta)
static inline Float computeDiffuseFresnel(Float ior, const int sampleCount)
{
    double diffuseFresnel = 0.0;
    Float fb = Fresnel::dielectricReflectance(ior, 0.0f);
    for (int i = 1; i <= sampleCount; ++i) {
        Float cosThetaSq = Float(i)/sampleCount;
        Float fa = Fresnel::dielectricReflectance(ior, min(std::sqrt(cosThetaSq), Float(1.0f)));
        diffuseFresnel += double(fa + fb)*(0.5/sampleCount);
        fb = fa;
    }

    return Float(diffuseFresnel);
}

}

}

#endif /* FRESNEL_HPP_ */
