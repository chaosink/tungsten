#ifndef SAMPLE_HPP_
#define SAMPLE_HPP_

#include "math/MathUtil.hpp"
#include "math/Angle.hpp"
#include "math/Vec.hpp"
#include "math/Box.hpp"

#include <algorithm>
#include <cmath>

namespace Tungsten {

namespace SampleWarp {

template<typename T>
static inline Float invertPhi(Vec<T, 3> w, Float mu)
{
    Float result = (w.x() == 0.0f && w.y() == 0.0f) ? mu*INV_TWO_PI : std::atan2(w.y(), w.x())*INV_TWO_PI;
    if (result < T(0))
        result += T(1);
    return result;
}

static inline Vec3f uniformHemisphere(const Vec2f &xi)
{
    Float phi  = TWO_PI*xi.x();
    Float r = std::sqrt(max(1.0f - xi.y()*xi.y(), Float(0.0f)));
    return Vec3f(std::cos(phi)*r, std::sin(phi)*r, xi.y());
}

static inline Float uniformHemispherePdf(const Vec3f &/*p*/)
{
    return INV_TWO_PI;
}

static inline Vec2f invertUniformHemisphere(const Vec3f &w, Float mu)
{
    return Vec2f(invertPhi(w, mu), w.z());
}

static inline Vec3f cosineHemisphere(const Vec2f &xi)
{
    Float phi = xi.x()*TWO_PI;
    Float r = std::sqrt(xi.y());

    return Vec3f(
            std::cos(phi)*r,
            std::sin(phi)*r,
            std::sqrt(max(1.0f - xi.y(), Float(0.0f)))
    );
}

static inline Float cosineHemispherePdf(const Vec3f &p)
{
    return std::abs(p.z())*INV_PI;
}

static inline Vec2f invertCosineHemisphere(const Vec3f &w, Float mu)
{
    return Vec2f(invertPhi(w, mu), max(1.0f - w.z()*w.z(), Float(0.0f)));
}

static inline Vec3f uniformDisk(const Vec2f &xi)
{
    Float phi = xi.x()*TWO_PI;
    Float r = std::sqrt(xi.y());
    return Vec3f(std::cos(phi)*r, std::sin(phi)*r, 0.0f);
}

static inline Float uniformDiskPdf()
{
    return INV_PI;
}

static inline Vec2f invertUniformDisk(const Vec3f &p, Float mu)
{
    return Vec2f(invertPhi(p, mu), p.xy().lengthSq());
}

static inline Vec3f uniformCylinder(Vec2f &xi)
{
    Float phi = xi.x()*TWO_PI;
    return Vec3f(
        std::cos(phi),
        std::sin(phi),
        xi.y()*2.0f - 1.0f
    );
}

static inline Float uniformCylinderPdf()
{
    return INV_PI;
}

static inline Vec3f uniformSphere(const Vec2f &xi)
{
    Float phi = xi.x()*TWO_PI;
    Float z = xi.y()*2.0f - 1.0f;
    Float r = std::sqrt(max(1.0f - z*z, Float(0.0f)));

    return Vec3f(
        std::cos(phi)*r,
        std::sin(phi)*r,
        z
    );
}

static inline Float uniformSpherePdf()
{
    return INV_FOUR_PI;
}

static inline Vec2f invertUniformSphere(const Vec3f &w, Float mu)
{
    return Vec2f(invertPhi(w, mu), (w.z() + 1.0f)*0.5f);
}

static inline Vec3f uniformSphericalCap(const Vec2f &xi, Float cosThetaMax)
{
    Float phi = xi.x()*TWO_PI;
    Float z = xi.y()*(1.0f - cosThetaMax) + cosThetaMax;
    Float r = std::sqrt(max(1.0f - z*z, Float(0.0f)));
    return Vec3f(
        std::cos(phi)*r,
        std::sin(phi)*r,
        z
    );
}

static inline Float uniformSphericalCapPdf(Float cosThetaMax)
{
    return INV_TWO_PI/(1.0f - cosThetaMax);
}

static inline bool invertUniformSphericalCap(const Vec3f &w, Float cosThetaMax, Vec2f &xi, Float mu)
{
    Float xiY = (w.z() - cosThetaMax)/(1.0f - cosThetaMax);
    if (xiY >= 1.0f || xiY < 0.0f)
        return false;

    xi = Vec2f(invertPhi(w, mu), xiY);
    return true;
}

static inline Vec3f phongHemisphere(const Vec2f &xi, Float n)
{
    Float phi = xi.x()*TWO_PI;
    Float cosTheta = std::pow(xi.y(), 1.0f/(n + 1.0f));
    Float r = std::sqrt(std::max(1.0f - cosTheta*cosTheta, Float(0.0f)));
    return Vec3f(std::cos(phi)*r, std::sin(phi)*r, cosTheta);
}

static inline Float phongHemispherePdf(const Vec3f &v, Float n)
{
    return INV_TWO_PI*(n + 1.0f)*std::pow(v.z(), n);
}

static inline Vec2f invertPhongHemisphere(const Vec3f &w, Float n, Float mu)
{
    return Vec2f(invertPhi(w, mu), std::pow(w.z(), n + 1.0f));
}

static inline Vec2f uniformTriangleUv(const Vec2f &xi)
{
    Float uSqrt = std::sqrt(xi.x());
    Float alpha = 1.0f - uSqrt;
    Float beta = (1.0f - xi.y())*uSqrt;

    return Vec2f(alpha, beta);
}

static inline Vec2f invertUniformTriangleUv(const Vec2f &uv)
{
    return Vec2f(sqr(1.0f - uv.x()), 1.0f - uv.y()/(1.0f - uv.x()));
}

static inline Vec3f uniformTriangle(const Vec2f &xi, const Vec3f& a, const Vec3f& b, const Vec3f& c)
{
    Vec2f uv = uniformTriangleUv(xi);
    return a*uv.x() + b*uv.y() + c*(1.0f - uv.x() - uv.y());
}

static inline Float uniformTrianglePdf(const Vec3f &a, const Vec3f &b, const Vec3f &c)
{
    return 2.0f/((b - a).cross(c - a).length());
}

static inline Float powerHeuristic(Float pdf0, Float pdf1)
{
    return (pdf0*pdf0)/(pdf0*pdf0 + pdf1*pdf1);
}

static inline Vec3f projectedBox(const Box3f &box, const Vec3f &direction, Float faceXi, const Vec2f &xi)
{
    Vec3f diag = box.diagonal();
    Float areaX = diag.y()*diag.z()*std::abs(direction.x());
    Float areaY = diag.z()*diag.x()*std::abs(direction.y());
    Float areaZ = diag.x()*diag.y()*std::abs(direction.z());

    Float u = faceXi*(areaX + areaY + areaZ);

    Vec3f result;
    if (u < areaX) {
        result.x() = direction.x() < 0.0f ? box.max().x() : box.min().x();
        result.y() = box.min().y() + diag.y()*xi.x();
        result.z() = box.min().z() + diag.z()*xi.y();
    } else if (u < areaX + areaY) {
        result.y() = direction.y() < 0.0f ? box.max().y() : box.min().y();
        result.z() = box.min().z() + diag.z()*xi.x();
        result.x() = box.min().x() + diag.x()*xi.y();
    } else {
        result.z() = direction.z() < 0.0f ? box.max().z() : box.min().z();
        result.x() = box.min().x() + diag.x()*xi.x();
        result.y() = box.min().y() + diag.y()*xi.y();
    }

    return result;
}

static inline Float projectedBoxPdf(const Box3f &box, const Vec3f &direction)
{
    Vec3f diag = box.diagonal();
    Float areaX = diag.y()*diag.z()*std::abs(direction.x());
    Float areaY = diag.z()*diag.x()*std::abs(direction.y());
    Float areaZ = diag.x()*diag.y()*std::abs(direction.z());

    return 1.0f/(areaX + areaY + areaZ);
}

static inline bool invertProjectedBox(const Box3f &box, const Vec3f &o, const Vec3f &d, Float &faceXi, Vec2f &xi, Float mu)
{
    Vec3f invD = Float(1.0f)/d;
    Vec3f relMin((box.min() - o));
    Vec3f relMax((box.max() - o));

    Float ttMin = 0, ttMax = 1e30f;
    int dim = -1;
    for (int i = 0; i < 3; ++i) {
        if (invD[i] >= 0.0f) {
            ttMin = max(ttMin, relMin[i]*invD[i]);
            Float x = relMax[i]*invD[i];
            if (x < ttMax) {
                ttMax = x;
                dim = i;
            }
        } else {
            Float x = relMin[i]*invD[i];
            if (x < ttMax) {
                ttMax = x;
                dim = i;
            }
            ttMin = max(ttMin, relMax[i]*invD[i]);
        }
    }

    if (ttMin <= ttMax) {
        Vec3f diag = box.diagonal();
        int dim1 = (dim + 1) % 3;
        int dim2 = (dim + 2) % 3;

        xi = Vec2f(
            (o[dim1] + d[dim1]*ttMax - box.min()[dim1])/diag[dim1],
            (o[dim2] + d[dim2]*ttMax - box.min()[dim2])/diag[dim2]
        );

        Float areaX = diag.y()*diag.z()*std::abs(d.x());
        Float areaY = diag.z()*diag.x()*std::abs(d.y());
        Float areaZ = diag.x()*diag.y()*std::abs(d.z());

        if (dim == 0)
            faceXi = mu*areaX/(areaX + areaY + areaZ);
        else if (dim == 1)
            faceXi = (areaX + mu*areaY)/(areaX + areaY + areaZ);
        else
            faceXi = (areaX + areaY + mu*areaZ)/(areaX + areaY + areaZ);


        return true;
    }
    return false;
}

}

}

#endif /* SAMPLE_HPP_ */
