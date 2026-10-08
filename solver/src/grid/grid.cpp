/*! \file grid.cpp
 * \author Yu Zhou
 */
#include "grid.h"
#include <stdexcept>

namespace DRL {

namespace GRID {

using std::invalid_argument;

Grid::Grid(const GridSpace &s, const va_idx_t &sz)
    : pS(0), S(s)
    , Sz(sz), Stride(calcStride(sz)), N(Stride[0] * sz[0])
    , adjI(0)
{
    if (Sz.size() != S.D) throw invalid_argument(
                "Grid::Grid() - Number of dimensions do not match grid size.");
}

Grid::Grid(const va_idx_t &sz, size_t q)
    : pS(new GridSpace(sz.size(), q))
    , S(*pS)
    , Sz(sz), Stride(calcStride(sz)), N(Stride[0] * sz[0])
    , adjI(0)
{

}

Grid::~Grid()
{
    delete pS;
    delete adjI;
}

} // namespace GRID

} // namespace DRL
