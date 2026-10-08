/*! \file gridgho.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "gridgho.h"

namespace DRL {

namespace GRID {

GridGho::GridGho(const GridSpace& gsp, const va_idx_t& sz, size_t ghoSz)
    : Grid(gsp, sz)
    , FG(Grid(gsp, sz + ghoSz * 2))
    , GhoSz(ghoSz)
    , i2fGs(0), i2fMap(0)

{

}

GridGho::GridGho(const va_idx_t &sz, size_t q, size_t ghoSz)
    : Grid(sz, q)
    , FG(Grid(S, sz + ghoSz * 2))
    , GhoSz(ghoSz)
    , i2fGs(0), i2fMap(0)
{

}

GridGho::~GridGho()
{
    delete i2fGs;
    delete i2fMap;
}

} // namespace GRID

} // namespace DRL

