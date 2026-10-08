/*! \file simulator.cpp
 * \author Yu Zhou
 */
#include "simulator.h"
#include "fileraw.h"
#include "roi.h"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace DRL {

namespace LBM {

using namespace std;

Simulator::Simulator(const FilePrm &prm)
    : p((prm.ValidateSimulation(), prm))
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
    for (size_t value : raw) p.getIVal(value);
    if (p.getReportDir() != 0) {
        const size_t axis = g.S.Qaxis[p.getReportDir()];
        for (size_t value : p.getReportSlice())
            if (value >= g.Sz[axis])
                throw invalid_argument("Report slice is outside the selected axis.");
    }
    buildIb();
}

void Simulator::writeTiming(ostream& log, const double* seconds, double wall,
                            size_t completedSteps) const
{
    double local[TIMING_CATEGORIES + 1], maximum[TIMING_CATEGORIES + 1];
    for (size_t i = 0; i < TIMING_CATEGORIES; ++i) local[i] = seconds[i];
    local[TIMING_CATEGORIES] = wall;
    mpi().AllReduce(local, maximum, TIMING_CATEGORIES + 1, MPI_DOUBLE, MPI_MAX);
    double fluidLocal = 0.0, fluidGlobal = 0.0;
    const auto& active = p.getFluidModel() == "color"
        ? p.getIdx({IT::VD, IT::GV, IT::VB, IT::PB, IT::HF})
        : p.getIdx({IT::VD, IT::GV, IT::VB, IT::PB});
    for (size_t idx : active) {
        const auto found = gIdx.find(idx);
        if (found != gIdx.end())
            fluidLocal += static_cast<double>(found->second.size());
    }
    mpi().AllReduce(&fluidLocal, &fluidGlobal, 1, MPI_DOUBLE, MPI_SUM);
    if (mpi().getId() != 0) return;
    // Build in a private stream so the caller's numeric format is preserved.
    ostringstream summary;
    summary << setprecision(17)
            << "LBFlowTiming {\"schema\":1,\"ranks\":" << mpi().getCnt()
            << ",\"completed_steps\":" << completedSteps
            << ",\"fluid_nodes\":" << fluidGlobal
            << ",\"total_voxels\":" << g.N
            << ",\"initialization_seconds_max\":" << maximum[INIT_TIME]
            << ",\"numerical_seconds_max\":" << maximum[NUMERICAL_TIME]
            << ",\"fields_diagnostics_seconds_max\":" << maximum[FIELD_TIME]
            << ",\"checkpoint_seconds_max\":" << maximum[CHECKPOINT_TIME]
            << ",\"run_seconds_max\":" << maximum[TIMING_CATEGORIES] << "}\n";
    log << summary.str();
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

va_val_t Simulator::getMRTS(double relaxationTime)
{
    const double rt = relaxationTime == 0.0 ? p.getRTime()[0] : relaxationTime;
    if (!std::isfinite(rt) || rt <= 0.5)
        throw invalid_argument("Relaxation time must be finite and greater than 0.5.");
    const double su = 1.0 / rt;
    const double sq = 8.0 * (2.0 - su) / (8.0 - su);
    if (g.S.Q == 9)
        return va_val_t({0, su, su, 0, sq, 0, sq, su, su});
    else if (g.S.Q == 19)
        return va_val_t({ 0, su, su, 0, sq, 0, sq, 0, sq, su, su, su, su,
                           su, su, su, sq, sq, sq });
    return va_val_t();
}

va_val_t Simulator::forceForIndex(size_t idx) const
{
    if (p.getIVal(idx).Type != IT::GV)
        return va_val_t(0.0, g.S.D);
    const auto& condition = p.getIVal(idx);
    if (condition.Dir >= g.S.Q || !std::isfinite(condition.Value))
        throw invalid_argument("Invalid body-force direction or magnitude.");
    if (p.getLSpeed() != 1.0)
        throw invalid_argument("The validated force scheme requires latticeSpeed=1.");
    return g.S.Qp[condition.Dir] * condition.Value;
}

va_val_t Simulator::forceSource(size_t i, const vava_val_t& velocity,
                                const va_val_t& force) const
{
    va_val_t source(0.0, g.S.Q);
    double uf = 0.0;
    for (size_t d = 0; d < g.S.D; ++d)
        uf += velocity[d][i] * force[d];
    for (size_t q = 0; q < g.S.Q; ++q) {
        double eu = 0.0, ef = 0.0;
        for (size_t d = 0; d < g.S.D; ++d) {
            eu += g.S.Qi[q][d] * velocity[d][i];
            ef += g.S.Qi[q][d] * force[d];
        }
        source[q] = wi[q] * (3.0 * (ef - uf) + 9.0 * eu * ef);
    }
    return source;
}

bool Simulator::updateConvergence(const vava_val_t& velocity, bool scalarsStable,
                                  double tolerance, ostream& os,
                                  const va_val_t* phaseIndicator)
{
    if (!hasReportVelocity) {
        previousReportVelocity = velocity;
        if (phaseIndicator) previousPhaseIndicator = *phaseIndicator;
        hasReportVelocity = true;
        stableReports = 0;
        return false;
    }
    double local[3] = {0.0, 0.0, 0.0};
    for (size_t d = 0; d < velocity.size(); ++d)
        for (size_t i = 0; i < velocity[d].size(); ++i) {
            const double value = velocity[d][i];
            if (!std::isfinite(value))
                throw runtime_error("Non-finite velocity in convergence report.");
            const double delta = value - previousReportVelocity[d][i];
            local[0] += delta * delta;
            local[1] += value * value;
            local[2] += 1.0;
        }
    double global[3];
    mpi().AllReduce(local, global, 3, MPI_DOUBLE, MPI_SUM);
    if (!std::isfinite(global[0]) || !std::isfinite(global[1]) || global[2] <= 0.0)
        throw runtime_error("Invalid global velocity residual.");
    const double rmsChange = sqrt(global[0] / global[2]);
    const double residual = global[1] > 0.0 ? sqrt(global[0] / global[1]) : 0.0;
    const bool velocityStable = rmsChange <= 1e-12 || residual < tolerance;
    double phaseChange = 0.0;
    if (phaseIndicator) {
        if (phaseIndicator->size() != previousPhaseIndicator.size())
            throw runtime_error("Phase-indicator convergence history has the wrong size.");
        double localChange = 0.0;
        for (size_t i = 0; i < phaseIndicator->size(); ++i) {
            if (!std::isfinite((*phaseIndicator)[i]))
                throw runtime_error("Non-finite phase indicator in convergence report.");
            localChange = max(localChange, abs((*phaseIndicator)[i] - previousPhaseIndicator[i]));
        }
        mpi().AllReduce(&localChange, &phaseChange, 1, MPI_DOUBLE, MPI_MAX);
        previousPhaseIndicator = *phaseIndicator;
    }
    previousReportVelocity = velocity;
    stableReports = scalarsStable && velocityStable && phaseChange < tolerance ? stableReports + 1 : 0;
    if (mpi().getId() == 0)
        os << " - velocityResidual = " << residual
           << ", velocityRmsChange = " << rmsChange
           << ", maxPhaseChange = " << phaseChange
           << ", stableReports = " << stableReports << "/3" << endl;
    return stableReports >= 3;
}

void Simulator::validatePopulations(const vava_val_t& populations) const
{
    for (const auto& direction : populations)
        for (double value : direction)
            if (!std::isfinite(value))
                throw runtime_error("Non-finite post-collision population; no restart will be written.");
}

void Simulator::prepareRestart(const string& path)
{
    if (mpi().getId() == 0) {
        ofstream data(path, ios::binary | ios::trunc);
        if (!data) throw runtime_error("Cannot create restart population file.");
        std::remove((path + ".meta").c_str());
    }
    mpi().Barrier();
}

void Simulator::saveRestartMetadata(const string& path, size_t time, size_t phases)
{
    mpi().Barrier();
    if (mpi().getId() == 0) {
        ofstream out(path + ".meta", ios::trunc);
        out << "DRLAB_OUTGOING_V1\n" << time << ' ' << phases << ' '
            << g.S.D << ' ' << g.S.Q << '\n';
        for (size_t n : g.Sz) out << n << ' ';
        out << '\n' << "NEXT_STREAM " << time + 1 << '\n'
            << "INITIAL_COLLISION_INCLUDED " << time + 1 << '\n'
            << "CONVERGENCE_HISTORY_RESET_ON_LOAD\n";
        if (!out) throw runtime_error("Cannot write restart metadata.");
    }
    mpi().Barrier();
}

void Simulator::validateRestart(const string& path, size_t time, size_t phases)
{
    if (mpi().getId() == 0) {
        ifstream in(path + ".meta");
        string format;
        size_t storedTime = 0, storedPhases = 0, storedD = 0, storedQ = 0;
        in >> format >> storedTime >> storedPhases >> storedD >> storedQ;
        if (!in || format != "DRLAB_OUTGOING_V1" || storedTime != time ||
            storedPhases != phases || storedD != g.S.D || storedQ != g.S.Q)
            throw runtime_error("Missing or incompatible restart metadata; historical unlabeled checkpoints are not accepted.");
        for (size_t n : g.Sz) {
            size_t storedN = 0;
            in >> storedN;
            if (!in || n != storedN) throw runtime_error("Restart grid dimensions differ.");
        }
        ifstream data(path, ios::binary | ios::ate);
        const size_t expected = g.N * g.S.Q * phases * sizeof(double);
        if (!data || data.tellg() != static_cast<streamoff>(expected))
            throw runtime_error("Restart population file has an invalid length.");
    }
    mpi().Barrier();
    hasReportVelocity = false;
    stableReports = 0;
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
            if (isIB) {
                ca /= nCA;
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
