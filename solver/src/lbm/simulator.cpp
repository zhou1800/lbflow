/*! \file simulator.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "simulator.h"
#include "fileraw.h"
#include "roi.h"

namespace DRL {

namespace LBM {

using namespace std;

Simulator::Simulator(const FilePrm &prm)
    : p(prm)
    , g(p.getD(), p.getQ())
    , raw(FileRaw<uint8_t>(p.getGFile()).Load(g))
    , gIdx(ROI::getRawIdx(raw))
    , ibI(g.LG.N), ibS(g.S.Q), ibA(g.S.Q), ibNorm(g.LG.N)
    , ibSinCA(g.LG.N), ibCosCA(g.LG.N)
    , ibF{veve_idx_t(g.S.Q), veve_idx_t(g.S.Q)}
    , output(p.getOutPath())
    , wi(getWi())
    , wiT(wi * 6.0 / p.getLSpeed())
{
    buildIb();
}

void Simulator::streaming(vava_val_t *f, vava_val_t *fOld, size_t tsp)
{
    static edge_gs_t sft = ROI::getEdge(g.LG, va_idx_t(g.GhoSz, g.S.D));

    for (size_t q = 1; q < g.S.Q; ++q) {
        size_t opq = g.S.Qop[q];
        // fill tmp with local processor
#if IS_MSVC
        va_val_t tmp1 = (*f)[q][sft[q].at(0)];
        (*fOld)[q][sft[opq].at(0)] = tmp1;
#else
        (*fOld)[q][sft[opq].at(0)] = (*f)[q][sft[q].at(0)];
#endif

        // fill tmp with adjacent processor
        for (size_t qi : g.S.Qadj[q]){    // streaming on direction i
            size_t opqi = g.S.Qop[qi];
            va_val_t sendData = (*f)[q][sft[q].at(qi)];
            double *pSend = &sendData[0];
            size_t size = sendData.size();
            va_val_t recvData(size);
            double *pRecv = &recvData[0];
            mpi().SendRecv(pSend, pRecv, size, MPI_DOUBLE, tsp,
                           g.AdjCi[opqi], g.AdjCi[qi]);
            (*fOld)[q][sft[opq].at(opqi)] = recvData;
        }
        (*f)[q].swap((*fOld)[q]);
    }
}

void Simulator::bounceBack(vava_val_t *f, const vava_val_t& fOld,
                           const veve_idx_t& ib)
{
    for (size_t q = 1; q < g.S.Q; ++q){
        size_t opq = g.S.Qop[q];
        for (size_t i : ib[q])
            (*f)[opq][i] = fOld[q][i];
    }
}

void Simulator::axisymmetric(vava_val_t *f)
{
    // Not Ready.
    //for (size_t q = 1; q < g.S.Q; ++q)
    //    for (size_t i : ibA[q])
    //        (*f)[g.S.Qop[q]][i] = (*f)[q][i];
}

void Simulator::initDist(vava_val_t *dist, double density,
                         const va_idx_t &ai) const
{
    va_val_t iDist = wi * density;
    for (size_t q = 0; q < g.S.Q; ++q)
        (*dist)[q][ai] = iDist[q];
}

va_idx_t Simulator::getElem(size_t d, size_t first, size_t last)
{
    static va_idx_t space = getSpace();
    va_idx_t start(size_t(0), g.S.D);
    va_idx_t size{g.Sz};
    start[d] = first;
    size[d] = last + 1 - first;
    gslice glbGs(Grid::CalcP2I(start, g.Stride), size, g.Stride);

    va_idx_t area;
    ROI::Intersect(g.LGs, g.N, glbGs, g.N, &area);
    return ROI::Intersect(space, area);
}

double Simulator::getElemCount(const va_idx_t &val)
{
    double cntT;
    double cnt = static_cast<double>(val.size());
    mpi().AllReduce(&cnt, &cntT, 1, MPI_DOUBLE, MPI_SUM);
    return cntT;
}

double Simulator::getKPar(size_t d)
{
    static double viscosity = (p.getRTime()[0] - 0.5) * p.getLSpeed() / 3;
    static double len = *(p.getReportSlice().end() - 1)
            - *(p.getReportSlice().begin());
    double area = 1;
    for (size_t di = 0; di < g.S.D; ++di)
        if (di != d) area *= g.Sz[di];
    return len * viscosity / area * pow(p.getLatticeLen(), 2) * 1.0e12;
}

va_val_t Simulator::getMRTM()
{
    if (g.S.Q == 9)
        return va_val_t(
        { 1, 1, 1, 1, 1, 1, 1, 1, 1,
         -4,-1,-1,-1,-1, 2, 2, 2, 2,
          4,-2,-2,-2,-2, 1, 1, 1, 1,
          0, 1, 0,-1, 0, 1,-1,-1, 1,
          0,-2, 0, 2, 0, 1,-1,-1, 1,
          0, 0, 1, 0,-1, 1, 1,-1,-1,
          0, 0,-2, 0, 2, 1, 1,-1,-1,
          0, 1,-1, 1,-1, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 1,-1, 1,-1});
    else if (g.S.Q == 19)
        return va_val_t(
        {  1,  1,  1,  1,  1,  1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         -30,-11,-11,-11,-11,-11,-11, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
          12, -4, -4, -4, -4, -4, -4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           0,  1, -1,  0,  0,  0,  0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0,
           0, -4,  4,  0,  0,  0,  0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0,
           0,  0,  0,  1, -1,  0,  0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1,
           0,  0,  0, -4,  4,  0,  0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1,
           0,  0,  0,  0,  0,  1, -1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1,
           0,  0,  0,  0,  0, -4,  4, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1,
           0,  2,  2, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1,-2,-2,-2,-2,
           0, -4, -4,  2,  2,  2,  2, 1, 1, 1, 1, 1, 1, 1, 1,-2,-2,-2,-2,
           0,  0,  0,  1,  1, -1, -1, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0,
           0,  0,  0, -2, -2,  2,  2, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0,
           0,  0,  0,  0,  0,  0,  0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
           0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1,
           0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0,
           0,  0,  0,  0,  0,  0,  0, 1,-1, 1,-1,-1, 1,-1, 1, 0, 0, 0, 0,
           0,  0,  0,  0,  0,  0,  0,-1,-1, 1, 1, 0, 0, 0, 0, 1,-1, 1,-1,
           0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 1, 1,-1,-1,-1,-1, 1, 1});
    return va_val_t();
}

va_val_t Simulator::getMRTS()
{
    static const double su = 1.0 / p.getRTime()[0]; // relaxtion frequency.
    static const double sq = 8.0 * (2.0 - su) / (8.0 - su);
    if (g.S.Q == 9)
        return va_val_t({1, su, su, 1, sq, 1, sq, su, su});
    else if (g.S.Q == 19)
        return va_val_t({ 0, su, su, 0, sq, 0, sq, 0, sq, su, su, su, su,
                           su, su, su, sq, sq, sq });
    return va_val_t();
}

void Simulator::buildIb()
{
    // 1.build extended raw;
    const va_idx_t& i2f = g.LG.getI2fMap();
    const vava_idx_t& adjI = g.LG.FG.getAdjI();
    va_u8_t exRaw(g.LG.FG.N);
    for (size_t i = 0; i < g.LG.N; ++i)
        exRaw[i2f[i]] = raw[i];
    syncBorder(&exRaw, 0);

    // 2.build ib, ibS, ibNorm
    size_t c = 1;
    static constexpr double pi = 3.1415926535897932384626433832795;
    for (size_t idx : p.getIdx({IT::VD, IT::PB, IT::VB, IT::GV, IT::HF})) {
        for (size_t i : gIdx[idx]) {
            bool isIB = false;
            ve_val_t norm(g.S.D, 0.0);
            double ca = 0;
            int nCA = 0;
            for (size_t q = 1; q < g.S.Q; ++q) {
                size_t ai = adjI[i2f[i]][q];
                FilePrm::IDXTYPE type = p.getIVal(exRaw[ai]).Type;
                if (type == IT::BB || type == IT::SD) { // solid boundary
                    ibS[q].push_back(i);
                    isIB = true;
                    for (size_t d = 0; d < g.S.D; ++d)
                        norm[d] -= g.S.Qi[q][d];
                    ca += p.getIVal(exRaw[ai]).Value;
                    nCA++;
                }
                else if (type == IT::AB)    // axisymmetrical boundary
                    ibA[q].push_back(i);
                else if (type == IT::HF){   // half permeable film
                    size_t phs = p.getIVal(exRaw[ai]).Phase;
                    ibF[phs][q].push_back(i);
                }
            }
            ca /= nCA;
            if (isIB) {
                ibI[i] = c;  // ib
                ibNorm[c] = move(norm);
                ibSinCA[c] = sin(ca / 180.0 * pi) ;
                ibCosCA[c] = cos(ca / 180.0 * pi) ;
                c++;
            }
        }
    }
    ibNorm.resize(c);
    ibSinCA.resize(c);
    ibCosCA.resize(c);
}

va_val_t Simulator::getWi()
{
    if (g.S.Q == 3)
        return va_val_t({4.0/6,
                        1.0/6, 1.0/6});
    else if (g.S.Q == 9)
        return va_val_t({4.0/9,
                        1.0/9,  1.0/9,  1.0/9,  1.0/9,
                        1.0/36, 1.0/36, 1.0/36, 1.0/36});

    else if (g.S.Q == 15)
        return va_val_t({1.0/8,
                        1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8, 1.0/8,
                        1.0/64, 1.0/64, 1.0/64, 1.0/64,
                        1.0/64, 1.0/64, 1.0/64, 1.0/64});

    else if (g.S.Q == 19)
        return va_val_t({1.0/3,
                        1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18,
                        1.0/36, 1.0/36, 1.0/36, 1.0/36,
                        1.0/36, 1.0/36, 1.0/36, 1.0/36,
                        1.0/36, 1.0/36, 1.0/36, 1.0/36});
    else
        throw invalid_argument("Simulator::geneWi - invalid number of"
                               " directions.");
}

va_idx_t Simulator::getSpace()
{
    va_idx_t space;
    for (auto& idx : p.getIdx({IT::VD, IT::GV, IT::VB, IT::PB, IT::HF}))
        space = ROI::Merge(space, gIdx[idx]);
    return space;
}

} // namespace LBM

} // namespace DRL
