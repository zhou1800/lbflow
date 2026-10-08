/*! \file imgproc.cpp
 * \author Yu Zhou
 */

#include "imgproc.h"

namespace DRL {

namespace PIPT {

using namespace std;

ImgProc::ImgProc(const GridMc &g)
    : pG(0)
    , G(g)
{
}

ImgProc::ImgProc(const va_idx_t &sz, size_t q, size_t ghoSz)
    : pG(new GridMc(sz, q, ghoSz))
    , G(*pG)
{

}

ImgProc::~ImgProc()
{
    delete pG;
}

} // namespace PIPT

} // namespace DRL
