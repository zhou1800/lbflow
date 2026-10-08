/*! \file gridmc.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "gridmc.h"
#include "fn.h"
#include <vector>
#include <algorithm>

namespace DRL {

namespace GRID {

using namespace std;

GridMc::GridMc(const GridSpace &gsp, const va_idx_t &sz, size_t ghoSz)
    : GridGho(gsp, sz, ghoSz)
    , CG(gsp, initCoreDistribution())
    , LCp(CalcI2P(mpi().getId(), CG.Sz))
    , AdjCi(initAdjCi())
    , LGs(initLclGs(LCp))
    , LGsWgl(addGho2Gs(LGs, *this))
    , LStart(CalcI2P(LGs.start(), Sz))
    , LG(gsp, LGs.size(), ghoSz)
    , subCp(), subGs(), subGsWgl(), subStart() , subG()
{

}

GridMc::GridMc(const va_idx_t &sz, size_t q, size_t ghoSz)
    : GridGho(sz, q, ghoSz)
    , CG(S, initCoreDistribution())
    , LCp(CalcI2P(mpi().getId(), CG.Sz))
    , AdjCi(initAdjCi())
    , LGs(initLclGs(LCp))
    , LGsWgl(addGho2Gs(LGs, *this))
    , LStart(CalcI2P(LGs.start(), Sz))
    , LG(S, LGs.size(), ghoSz)
    , subCp(), subGs(), subGsWgl(), subStart() , subG()
{

}

GridMc::~GridMc()
{
}

va_idx_t GridMc::initAdjCi() const
{
    va_idx_t a(S.Q);
    for (size_t q = 0; q < S.Q; ++q)
        a[q] = calcAdjI(mpi().getId(), q, CG);
    return a;
}

va_idx_t GridMc::initCoreDistribution() const
{
    vector<double> sz(S.D);  // use vector for iterator.
    for (size_t i = 0; i < S.D; ++i)
        sz[i] = static_cast<double>(Sz[i]);

    va_idx_t procDist(1, S.D);

    list<size_t> procPrime {ALG::Prime(mpi().getCnt())};
    for (auto prim : procPrime) {
        size_t iMax
        {static_cast<size_t>(max_element(sz.begin(), sz.end()) - sz.begin())};
        if (sz[iMax] < prim)
            throw std::invalid_argument("too many processors");
        procDist[iMax] *= prim;
        sz[iMax] /= prim;
    }
    return procDist;
}

gslice GridMc::initLclGs(const va_idx_t &cP) const
{
    // invert the mode to ensure first core always get the smallest size.
    va_idx_t mode {CG.Sz - Sz % CG.Sz};
    va_idx_t length {Sz / CG.Sz};      // floor length.

    // For consistency, caculation of start is based on floor length.
    va_idx_t start((size_t)(0), S.D);
    for (size_t i = 0; i < S.D; ++i) {
        size_t ci {cP[i]};
        size_t nCell {ci * length[i]};
        nCell += mode[i] >= ci ? 0 : ci - mode[i];
        start[i] += nCell;
    }
    start *= Stride;

    for (size_t i = 0; i < S.D; ++i)     // caculate actual length
        length[i] += mode[i] > cP[i] ? 0 : 1;

    return gslice(accumulate(begin(start), end(start), static_cast<size_t>(0)),
                  move(length), Stride);
}

} // namespace GRID

} //namespace DRL
