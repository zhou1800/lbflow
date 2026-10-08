/*! \file roi.cpp
 * \author Yu Zhou
 */
#include "roi.h"
#include <algorithm>
#include <stdexcept>

namespace DRL {

namespace GRID {

using namespace std;

ROI::ROI()
{

}

va_idx_t ROI::Intersect(va_idx_t &idx, const gslice &gs, size_t gsN,
                         va_idx_t *pos)
{
    va_idx_t gsIdx = gs2idx(gs, gsN);
    return Intersect(gsIdx, idx, pos);
}

va_idx_t ROI::Intersect(const gslice &gs1, size_t gsN1,
                         const gslice &gs2, size_t gsN2, va_idx_t *pos)
{
    va_idx_t gsIdx1 = gs2idx(gs1, gsN1);
    va_idx_t gsIdx2 = gs2idx(gs2, gsN2);
    return Intersect(gsIdx1, gsIdx2, pos);
}

va_idx_t ROI::Intersect(va_idx_t &idx1, va_idx_t &idx2, va_idx_t *pos)
{
    va_idx_t idxF(min(idx1.size(), idx2.size()));
    sort(begin(idx1), end(idx1));
    sort(begin(idx2), end(idx2));

    va_idx_t idx;
    if (pos == 0){
        auto it = set_intersection(begin(idx1), end(idx1),
                                   begin(idx2), end(idx2), begin(idxF));
        size_t len = it - begin(idxF);
        idx = va_idx_t(&(idxF[0]), len);
    }
    else {
        pos->resize(idxF.size());

        size_t p1 = 0, p2 = 0, pp = 0;
        while (p1 < idx1.size() && p2 < idx2.size()){
            if (idx1[p1] < idx2[p2]) ++p1;
            else if (idx1[p1] > idx2[p2]) ++p2;
            else{
                (*pos)[pp] = p1;
                idxF[pp] = idx1[p1];
                ++pp; ++p1; ++p2;
            }
        }
        size_t len = pp;
        idx = va_idx_t(&(idxF[0]), len);
        *pos = va_idx_t(&((*pos)[0]), len);
    }

    return idx;
}

va_idx_t ROI::Merge(va_idx_t &idx1, va_idx_t &idx2)
{
    sort(begin(idx1), end(idx1));
    sort(begin(idx2), end(idx2));

    va_idx_t result(idx1.size() + idx2.size());
    merge(begin(idx1), end(idx1), begin(idx2), end(idx2), begin(result));
    return result;
}

va_idx_t ROI::gs2idx(const gslice &gs, size_t gsN)
{
    va_idx_t gd(gsN);
    for (size_t i = 0; i < gsN; ++i) gd[i] = i;
    return gd[gs];
}

va_gslice ROI::getIBd(const GridGho &g)
{
    vava_idx_t marker(g.S.D);  // marker for inner boundary
    for (size_t d = 0; d < g.S.D; ++d)
        marker[d] = va_idx_t({g.GhoSz,g.GhoSz, g.FG.Sz[d] - g.GhoSz* 2});

    va_gslice iBd(g.S.Q);
    for (size_t q = 0; q < g.S.Q; ++q)
        iBd[q] = gslice(calcStart(g, q, marker), calcSize(g, q), g.FG.Stride);

    return iBd;
}

va_gslice ROI::getOBd(const GridGho &g)
{
    vava_idx_t marker(g.S.D);  // marker for outterBD
    for (size_t d = 0; d < g.S.D; ++d)
        marker[d] = va_idx_t({0, g.GhoSz, g.FG.Sz[d] - g.GhoSz});

    va_gslice oBd(g.S.Q);
    for (size_t q = 0; q < g.S.Q; ++q)
        oBd[q] = gslice(calcStart(g, q, marker), calcSize(g, q), g.FG.Stride);

    return oBd;
}

va_gslice ROI::getSft(const GridGho &g)
{
    vava_idx_t marker({0, g.GhoSz, g.GhoSz * 2}, g.S.D); // marker for shift

    va_gslice sft(g.S.Q);
    for (size_t q = 0; q < g.S.Q; ++q)
        sft[q] = gslice(calcStart(g, q, marker), g.Sz, g.FG.Stride);

    return sft;
}

edge_gs_t ROI::getEdge(const Grid &g, const va_idx_t &sftSz)
{
    // parameter check.
    if (g.Sz.size() != sftSz.size())
        throw invalid_argument("ROI::getEdge - dimension of shiftSz not equal"
                               " to grid dimension.");
    for (size_t d = 0; d < g.S.D; ++d)
        if (g.Sz[d] < sftSz[d])
            throw invalid_argument("ROI::getEdge - sftSz greater than grid"
                                   " size.");

    valarray<map<size_t, gslice>> edge(g.S.Q);
    for (size_t qM = 0; qM < g.S.Q; ++qM) {
        va_idx_t div = calcDivider(g, sftSz, qM);
        edge[qM][0] = calcSlice(g, div, qM, 0);
        for (size_t qS : g.S.Qadj[qM]) {
            edge[qM][qS] = calcSlice(g, div, qM, qS);
        }
    }
    return edge;
}

size_t ROI::calcStart(const GridGho& g, size_t q, const vava_idx_t &marker)
{
    va_idx_t start(g.S.D);
    for (size_t d = 0; d < g.S.D; ++d)
        start[d] = marker[d][g.S.Qi[q][d] + 1];
    return g.CalcP2I(start, g.FG.Stride);
}

va_idx_t ROI::calcSize(const GridGho& g, size_t q)
{
    va_idx_t size(g.S.D);
    for (size_t d = 0; d < g.S.D; ++d)
        size[d] = g.S.Qi[q][d] == 0 ? g.Sz[d] : g.GhoSz;
    return size;
}

va_idx_t ROI::calcDivider(const Grid &g, const va_idx_t &sftSz, size_t q)
{
    va_idx_t div(static_cast<size_t>(0), g.S.D);
    for (size_t d = 0; d < g.S.D; ++d){
        if      (g.S.Qi[q][d] ==  1) div[d] = g.Sz[d] - sftSz[d];
        else if (g.S.Qi[q][d] == -1) div[d] = sftSz[d];
    }
    return div;
}

gslice ROI::calcSlice(const Grid &g, const va_idx_t &div, size_t qM, size_t qS)
{
    va_idx_t start(static_cast<size_t>(0), g.S.D);
    va_idx_t size(g.Sz);
    for (size_t d = 0; d < g.S.D; ++d){
        if (g.S.Qi[qM][d] ==  1){       // x+, y+, z+
            if (g.S.Qi[qS][d] == 0){        // from 0 to divider.
                start[d] = 0;
                size[d] = div[d];
            }
            else if (g.S.Qi[qS][d] == 1) {  // from divider to end.
                start[d] = div[d];
                size[d] = g.Sz[d] - div[d];
            }
            else throw invalid_argument
                    ("ROISft::calcSlice - qS is not subject to qM.");
        }
        else if (g.S.Qi[qM][d] == -1) { // x-, y-, z-
            if (g.S.Qi[qS][d] == 0){        //from divider to end.
                start[d] = div[d];
                size[d] = g.Sz[d] - div[d];
            }
            else if (g.S.Qi[qS][d] == -1) { // from 0 to divider.
                start[d] = 0;
                size[d] = div[d];
            }
            else throw invalid_argument
                    ("ROISft::calcSlice - qS is not subject to qM.");
        }
    }
    return gslice(g.CalcP2I(start, g.Stride), size, g.Stride);
}

} // namespace GRID

} // namespace DRL
