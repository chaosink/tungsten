#ifndef DISTRIBUTION1D_HPP_
#define DISTRIBUTION1D_HPP_

#include "math/MathUtil.hpp"

#include <algorithm>
#include <vector>

namespace Tungsten {

class Distribution1D
{
    std::vector<Float> _pdf;
    std::vector<Float> _cdf;
public:
    Distribution1D(std::vector<Float> weights)
    : _pdf(std::move(weights))
    {
        _cdf.resize(_pdf.size() + 1);
        _cdf[0] = 0.0f;
        for (size_t i = 0; i < _pdf.size(); ++i)
            _cdf[i + 1] = _cdf[i] + _pdf[i];
        Float totalWeight = _cdf.back();
        for (Float &p : _pdf)
            p /= totalWeight;
        for (Float &c : _cdf)
            c /= totalWeight;
        _cdf.back() = 1.0f;
    }

    void warp(Float &u, int &idx) const
    {
        idx = int(std::distance(_cdf.begin(), std::upper_bound(_cdf.begin(), _cdf.end(), u)) - 1);
        u = clamp((u - _cdf[idx])/_pdf[idx], Float(0.0f), Float(1.0f));
    }

    Float pdf(int idx) const
    {
        return _pdf[idx];
    }

    Float cdf(int idx) const
    {
        return _cdf[idx];
    }
};

}

#endif /* DISTRIBUTION1D_HPP_ */
