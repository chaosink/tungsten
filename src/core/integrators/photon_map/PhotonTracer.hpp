#ifndef PHOTONTRACER_HPP_
#define PHOTONTRACER_HPP_

#include "PhotonMapSettings.hpp"
#include "FrustumBinner.hpp"
#include "PhotonRange.hpp"
#include "KdTree.hpp"
#include "Photon.hpp"

#include "integrators/TraceBase.hpp"

#include <unordered_map>

namespace Tungsten {

namespace Bvh {
class BinaryBvh;
}
class GridAccel;

class HashedShadowCache
{
    std::unordered_map<uint64, Float> _cache;

public:
    HashedShadowCache(uint64 initialSize)
    {
        _cache.reserve(initialSize);
    }

    void clear()
    {
        _cache.clear();
    }

    template<typename Tracer>
    inline Float hitDistance(uint32 photon, uint32 bin, Tracer tracer)
    {
        uint64 key = (uint64(photon) << 32ull) | uint64(bin);
        auto iter = _cache.find(key);
        if (iter == _cache.end()) {
            Float dist = tracer();
            _cache.insert(std::make_pair(key, dist));
            return dist;
        } else {
            return iter->second;
        }
    }
};

class LinearShadowCache
{
    const uint32 MaxCacheBins = 1024*1024;

    std::unique_ptr<uint32[]> _photonIndices;
    std::unique_ptr<Float[]> _distances;

public:
    LinearShadowCache()
    : _photonIndices(new uint32[MaxCacheBins]),
      _distances(new Float[MaxCacheBins])
    {
        clear();
    }

    void clear()
    {
        std::memset(_photonIndices.get(), 0, MaxCacheBins*sizeof(_photonIndices[0]));
    }

    template<typename Tracer>
    inline Float hitDistance(uint32 photon, uint32 bin, Tracer tracer)
    {
        if (bin < MaxCacheBins && _photonIndices[bin] == photon) {
            return _distances[bin];
        } else {
            Float dist = tracer();
            _photonIndices[bin] = photon;
            _distances[bin] = dist;
            return dist;
        }
    }
};

class PhotonTracer : public TraceBase
{
    PhotonMapSettings _settings;
    uint32 _mailIdx;
    std::unique_ptr<const Photon *[]> _photonQuery;
    std::unique_ptr<Float[]> _distanceQuery;
    std::unique_ptr<uint32[]> _mailboxes;

    LinearShadowCache _directCache;
    HashedShadowCache _indirectCache;

    FrustumBinner _frustumGrid;

    void clearCache();

public:
    PhotonTracer(TraceableScene *scene, const PhotonMapSettings &settings, uint32 threadId);

    void evalPrimaryRays(const PhotonBeam *beams, const PhotonPlane0D *planes0D, const PhotonPlane1D *planes1D,
            uint32 start, uint32 end, Float radius, const Ray *depthBuffer, PathSampleGenerator &sampler, Float scale);

    Vec3f traceSensorPath(Vec2u pixel, const KdTree<Photon> &surfaceTree,
            const KdTree<VolumePhoton> *mediumTree, const Bvh::BinaryBvh *mediumBvh, const GridAccel *mediumGrid,
            const PhotonBeam *beams, const PhotonPlane0D *planes0D, const PhotonPlane1D *planes1D, PathSampleGenerator &sampler,
            Float gatherRadius, Float volumeGatherRadius,
            PhotonMapSettings::VolumePhotonType photonType, Ray &depthRay, bool useFrustumGrid);

    void tracePhotonPath(SurfacePhotonRange &surfaceRange, VolumePhotonRange &volumeRange,
            PathPhotonRange &pathRange, PathSampleGenerator &sampler);
};

}

#endif /* PHOTONTRACER_HPP_ */
