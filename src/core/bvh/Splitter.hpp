#ifndef SPLITTER_HPP_
#define SPLITTER_HPP_

#include "math/Box.hpp"

namespace Tungsten {

namespace Bvh {

struct SplitInfo
{
    Box3fp lBox, rBox;
    Box3fp lCentroidBox, rCentroidBox;
    int dim;
    uint32 idx;
    Float cost;
};

namespace Splitter {
    static CONSTEXPR Float IntersectionCost = 1.0f;
    static CONSTEXPR Float TraversalCost = 1.0f;
};

}

}

#endif /* SPLITTER_HPP_ */
