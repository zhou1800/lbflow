/*! \file pipt/markedge.cpp
 * \author Yu Zhou
 */

#include "markedge.h"

namespace DRL {

namespace PIPT {

using namespace DRL::GRID;

MarkEdge::MarkEdge(const GridMc &g)
    : ImgProc(g)
{
}

MarkEdge::MarkEdge(const va_idx_t &sz, size_t q, size_t ghoSz)
    : ImgProc(sz, q, ghoSz)
{

}

} // namespace PIPT

} // namespace DRL
