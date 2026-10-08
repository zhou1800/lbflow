#include "tgrid.h"
#include "roi.h"
#include <sstream>

using namespace std;

void tROI(ostream &os, const GridGho &g, const va_u8_t &raw)
{
    va_gslice ibd  = ROI::getIBd(g);
    va_gslice obd  = ROI::getOBd(g);
    va_gslice sft  = ROI::getSft(g);
    edge_gs_t  edge = ROI::getEdge(g, va_idx_t(g.GhoSz, g.S.D));
    map_idx_t  idx  = ROI::getRawIdx(raw);

    Debugger dbg(os);
    dbg.Add("ibd", ibd);
    dbg.Add("ibd", obd);
    dbg.Add("ibd", sft);
    for (size_t qm = 0; qm < edge.size(); ++qm)
        for (auto qs : edge[qm]){
            stringstream ss;
            ss << "edge(" << qm << ", " << qs.first << ")";
            dbg.Add(ss.str(), edge[qm][qs.first]);
        }
    for (auto ai : idx){
        stringstream ss;
        ss << "RawIdx(" << ai.first << ")";
        dbg.Add(ss.str(), idx[ai.first]);
    }

    va_idx_t ai1 = {1,3,5,6,7,9,14};
    va_idx_t ai2 = {2,4,5,7,10};
    va_idx_t pos;
    va_idx_t interIdxIdx = ROI::Intersect(ai1, ai2, &pos);
    dbg.Add("interIdxIdx",interIdxIdx);
    dbg.Add("interIdxIdxPos", pos);
    va_idx_t mergeIdxIdx = ROI::Merge(ai1, ai2);
    dbg.Add("mergeIdxIdx", mergeIdxIdx);

    va_idx_t interGsGs;
    if (g.S.Q == 19)
        interGsGs = ROI::Intersect(edge[10][4], g.N, edge[16][4], g.N);
    dbg.Add("interGsGs", interGsGs);

    va_idx_t interIdxGs;
    if (g.S.Q == 19)
        interIdxGs = ROI::Intersect(ai1, edge[16][4], g.N);
    dbg.Add("interIdxGs", interIdxGs);

    dbg.Show();
}
