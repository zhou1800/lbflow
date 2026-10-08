/*! \file scsp.cpp
 * \author Yu Zhou
 */
#include "scsp.h"
#include "logger.h"
#include "fileraw.h"
#include "fn.h"
#include <sstream>
#include <iomanip>
#include <cfloat>
#include <algorithm>

namespace DRL {

namespace LBM {

using namespace DRL::HLP;

using namespace DRL::ALG;
using namespace std;

SCSP::SCSP(const FilePrm &prm)
    : Simulator(prm)
    , f(g.S.Q), fOld(g.S.Q)
    , v(g.S.D)
    , dens(0.0, g.LG.N)
{
    if (p.getEquilModel() != "SRT" && p.getEquilModel() != "MRT")
        throw invalid_argument("Invalid equilibrium model.");
    if (g.S.Q != 9 && g.S.Q != 19)
        throw invalid_argument("Invalid DmQn model for SCSP.");

    for (size_t q = 0; q < g.S.Q; ++q){
        f[q] = va_val_t(0.0, g.LG.N);
        fOld[q] = va_val_t(0.0, g.LG.N);
    }
    for (size_t d = 0; d < g.S.D; ++d)
        v[d]  = va_val_t(0.0, g.LG.N);

    output.AddScalar("density", &dens);
    for (size_t d = g.S.D; d != 0; --d)
        output.AddVector("velocity", &v[d-1]);
    output.Init(g.S.D, g.N, g.Sz, g.Stride, g.LG.Sz, g.LStart);
}

void SCSP::Run(std::ostream &logProc)
{
    Logger procLog(logProc);
    procLog.Reset(true);

    procLog.Log("Initializing...", false);
    if (p.getTStart() != 0){
        loadRst(p.getTStart());
    }
    else {
        initialize();
        output.Save(p.getTStart());
    }
    procLog.Log("Done!\t");

    for (size_t t = p.getTStart() + 1; t <= p.getTEnd(); ++t) {
        streaming(&f, &fOld, t);
        bounceBack(&f, fOld, ibS);
        //axisymmetric(&f);
        zouHeBoundary();
        calcMacroVar();
        if (t % p.getTStep() == 0) {
            stringstream ss;
            ss << "Step " << t << " - ";
            procLog.Log(ss.str());
            output.Save(t);
            if (report(logProc, dens, v)){
                saveRst(t);
                break;
            }
        }
        else if (t % p.getTRpt() == 0){
            stringstream ss;
            ss << "Step " << t << " - ";
            procLog.Log(ss.str());
            report(logProc, dens, v);
        }
        collision();
        if (t % p.getTRst() == 0) saveRst(t);
    }
}

void SCSP::initialize()
{
    static const double rho = p.getDensity()[0];
    static const va_val_t iDist = wi * rho;
    for (auto& idx : p.getIdx({IT::VD, IT::GV, IT::VB, IT::PB}))
        for (size_t i : gIdx[idx]){
            for (size_t q = 0; q < g.S.Q; ++q)
                f[q][i] = iDist[q];
                dens[i] = rho;
            // vel = 0 by default.
        }
}

void SCSP::calcMacroVar()
{
    static const double dens0 = p.getDensity()[0];
    for (size_t idx : p.getIdx({IT::VD, IT::GV})) {
        for(size_t i : gIdx[idx]) {
            dens[i] = f[0][i];
            for (size_t q = 1; q < g.S.Q; ++q)
                dens[i] += f[q][i];

            for (size_t d = 0; d < g.S.D; ++d) {
                v[d][i] = 0;
                for (size_t q = 1; q < g.S.Q; ++q)   // Qi[0] = {0, 0, 0}
                    if (g.S.Qi[q][d] == 1) v[d][i] += f[q][i];
                    else if (g.S.Qi[q][d] == -1) v[d][i] -= f[q][i];
#if ISCOMPRESSABLE
                v[d][i] /= dens[i];
#else
                v[d][i] /= dens0;
#endif
            }
        }
    }
}

void SCSP::collision()
{
    // incorporate gravity
    for (size_t idx : p.getIdx({IT::GV})) {
        va_val_t fi = g.S.Qp[p.getIVal(idx).Dir]
                * p.getRTime()[0] * p.getIVal(idx).Value;
        for (size_t d = 0; d < g.S.D; ++d)
            for (size_t i : gIdx[idx])
                v[d][i] += fi[d] / dens[i];
    }
    static const EQUILFN equilfn = p.getEquilModel() == "SRT" ? SRT : MRT;

    switch (equilfn){
    case SRT:       equilSRT();  break;
    case MRT:       equilMRT();  break;
    }
}

void SCSP::zouHeBoundary()
{
    // Apply Flux Boundary
    for (size_t idx : p.getIdx({IT::VB})) {
        size_t dir = p.getIVal(idx).Dir;      // flow direction
        double velocity = p.getIVal(idx).Value;
        velocity = (g.S.Qsign[dir] == 1 ? velocity : - velocity);
        double vel_1 = 1.0 / (1.0 - velocity);

        size_t qAxis = g.S.Qaxis[dir];
        va_idx_t qOpRel = g.S.Qrel[g.S.Qop[dir]];
        va_idx_t qNorm = g.S.Qprep[dir];

        for (size_t i : gIdx[idx]){
            double sumF = 0;
            for (size_t q : qOpRel)
                sumF += f[q][i] + f[q][i];

            for (size_t q : qNorm)
                sumF += f[q][i];

            for (size_t d = 0; d < g.S.D; ++d)
                if (d != qAxis) v[d][i] = 0;
                else {
                    v[d][i] = velocity;
                    dens[i] = sumF * vel_1;
                }
        }
    }

    // Apply Pressure boundary
    for (size_t idx : p.getIdx({IT::PB})) {
        size_t dir = p.getIVal(idx).Dir;      // flow direction
        double density = p.getIVal(idx).Value;

        size_t qAxis = g.S.Qaxis[dir];
        va_idx_t qOpRel = g.S.Qrel[g.S.Qop[dir]];
        va_idx_t qNorm = g.S.Qprep[dir];
        int qSign = g.S.Qsign[dir];

        for (size_t i : gIdx[idx]){
            double sumF = 0;
            for (size_t q : qOpRel)
                sumF += f[q][i] + f[q][i];

            for (size_t q : qNorm)
                sumF += f[q][i];

            dens[i] = density;
            for (size_t d = 0; d < g.S.D; ++d)
                if (d != qAxis) v[d][i] = 0;
                else {
                    sumF /= density; sumF--;
                    v[d][i] = qSign == 1 ? - sumF : sumF;
                }
        }
    }

    // Update Unknown distribution functions
    for (size_t idx : p.getIdx({IT::PB, IT::VB})) {
        size_t dir = p.getIVal(idx).Dir;      // flow direction

        size_t qAxis = g.S.Qaxis[dir];
        va_idx_t qRel = g.S.Qrel[dir];
        va_idx_t qPrep = g.S.Qprep[dir];
        int qSign = g.S.Qsign[dir];

        for (size_t i : gIdx[idx])
            for (size_t q : qRel) {  // for each unknown dist[q]
                // calculate distribution function.
                // 1. dist = dist_opp
                f[q][i] = f[g.S.Qop[q]][i];
                // 2. dist += 6 * wi * rho * u_e / c
                f[q][i] += qSign * wiT[q] * dens[i] * v[qAxis][i];
                // 3. dist += 1/2 * sum(F_normal)
                double sumF = 0;
                for (size_t qn : qPrep)
                    for (size_t d = 0; d < g.S.D; ++d)
                        if (g.S.Qi[qn][d] == 0 || g.S.Qi[q][d] == 0)
                            continue;
                        else if (g.S.Qi[qn][d] == g.S.Qi[q][d])
                            sumF += f[qn][i];
                        else
                            sumF -= f[qn][i];
                f[q][i] -= sumF * 0.5;
            }
    }
}

void SCSP::loadRst(std::size_t t)
{
    // prepare filename
    ostringstream os;
    os << p.getOutPath() << "." << setw(10) << setfill('0') << t << ".rst";
    FileRaw<double> file(os.str(), true);

    // load distribution functions.
    for (size_t q = 0,s=0; q < g.S.Q; ++q,s += g.LG.N)
        f[q] = file.Load(g, t);
}

void SCSP::saveRst(std::size_t t)
{
    // prepare filename
    ostringstream os;
    os << p.getOutPath() << "." << setw(10) << setfill('0') << t << ".rst";
    FileRaw<double> file(os.str(), true);

    // save distribution function.
    for (size_t q = 0,s=0; q < g.S.Q; ++q,s += g.LG.N)
        file.Save(g, f[q], t);
}

bool SCSP::report(ostream &os, const va_val_t &dens, const vava_val_t &v)
{
    static const size_t rptDir = p.getReportDir();
    if (rptDir == 0)  // no report
        return false;
    try{
        static const size_t d = g.S.Qaxis[rptDir];  // report axis.
        static const ve_idx_t slices = p.getReportSlice();   // report slices.
        static const auto fl = minmax_element(slices.begin(), slices.end());
        static const va_idx_t lclRange = getElem(d, *fl.first, *fl.second);
        static const vava_idx_t lclSlice = {  // first and last slice.
            getElem(d, *fl.first, *fl.first),
            getElem(d, *fl.second, *fl.second)};
        static const double cntS[] = {      // number of elements in f&l slice.
            getElemCount(lclSlice[0]),
            getElemCount(lclSlice[1])};

        static const double MaxErr = p.getStopError();
        static const double kPar = 3.0 * getKPar(d);
        static double KOld = 0;                // Permeability of last step.

        // 1. Accumulate total velocity of range.
        double velR = 0;             // Total sum velocity of range.
        {
            double velLR = 0;             // local sum velocity of range.
            for (auto a : lclRange)
                velLR += v[d][a];
            mpi().AllReduce(&velLR, &velR, 1, MPI_DOUBLE, MPI_SUM);
        }

        // 2. Accumulate density of first and last slice.
        double rhoS[2] = {0, 0};    // Total rho[slice] of slice.
        for (size_t i = 0; i < 2; ++i){
            double rhoLS = 0;
            for (auto a : lclSlice[i])
                rhoLS += dens[a];
            mpi().AllReduce(&rhoLS, &rhoS[i], 1, MPI_DOUBLE, MPI_SUM);
        }

        // 3. Calculate K.
        double dP = rhoS[0] / cntS[0] - rhoS[1] / cntS[1];
        double Q = velR / double(*fl.second - *fl.first + 1);
        double K = kPar * Q / dP;
        if (K == 0) return false;
        double errK = abs(K - KOld) / K;
        KOld = K;

        if (std::isnan(dP) || std::isnan(Q) ||
            std::isnan(K)  || std::isnan(errK))
            errK = 0;       // end simulation!

        // 4. Report.
        if (mpi().getId() == 0){
            os << setprecision(6) << scientific
               << " - dP = " << dP << "(mu/lu/ts^2),"
               << " Q = " << Q << "(lu^3/ts)" << endl
               << " - K = " << K << "(um^2), "
               << fixed << "errK = " << errK * 100.0 << "%" << endl;
        }
        return abs(errK) < MaxErr;
    }
    catch (exception& e) { os << e.what(); }
    return false;
}

void SCSP::equilSRT()
{
    static const double c = p.getLSpeed();
    static const double c1 = 3.0 / (c * c );
    static const double c2 = 4.5 / (c * c * c * c);
    static const double c3 = 1.5 / (c * c);
    static const double omega = 1.0 / p.getRTime()[0]; // relaxtion frequency.
    static const double dens0 = p.getDensity()[0];

    for (size_t idx : p.getIdx({IT::VD, IT::GV, IT::VB, IT::PB})) {
        for (size_t i : gIdx[idx]) {     // array index
            double usq = 0;
            for (size_t d = 0; d < g.S.D; ++d)
                usq += v[d][i] * v[d][i];

            for (size_t q = 0; q < g.S.Q; ++q) {
                double ui = 0;
                for (size_t d = 0; d < g.S.D; ++d)
                    ui += v[d][i] * g.S.Qi[q][d];

#if ISCOMPRESSABLE
                double feq = wi[q] * dens[i] * (1 + c1*ui + c2*ui*ui- c3*usq);
#else
                double feq = wi[q] *
                     (dens[i] + dens0 * (c1*ui + c2*ui*ui- c3*usq));
#endif
                f[q][i] += (feq - f[q][i]) * omega;
            }
        }
    }
}

void SCSP::equilMRT()
{
    static const double c = p.getLSpeed();
    static const double c1 = 3.0 / (c * c );
    static const double c2 = 4.5 / (c * c * c * c);
    static const double c3 = 1.5 / (c * c);
    static const double dens0 = p.getDensity()[0];

    static const size_t Q = g.S.Q;  // for short
    static const va_val_t M {getMRTM()};
    static const va_val_t MInv = InvOrthoRow(M, Q, Q);
    static const va_val_t S {getMRTS()};
    static const va_val_t MInvS = MatMulDiag(MInv, S, Q);
    static const va_val_t MInvSM = MatMulMat(MInvS, M, Q, Q, Q);
    static const va_val_t MInvSMT = MatTranspose(MInvSM, Q, Q);

    for (size_t idx : p.getIdx({IT::VD, IT::GV, IT::VB, IT::PB})) {
        for (size_t i : gIdx[idx]) {     // array index

            // 1. pre calculated MInvSM.
            double usq = 0;
            for (size_t d = 0; d < g.S.D; ++d)
                usq += v[d][i] * v[d][i];

            va_val_t ff(0.0, Q);
            auto im = begin(MInvSMT);
            for (size_t q = 0; q < g.S.Q; ++q) {
                double ui = 0;
                for (size_t d = 0; d < g.S.D; ++d)
                    ui += v[d][i] * g.S.Qi[q][d];

#if ISCOMPRESSABLE
                double feq = wi[q] * dens[i] * (1 + c1*ui + c2*ui*ui- c3*usq);
#else
                double feq = wi[q] *
                     (dens[i] + dens0 * (c1*ui + c2*ui*ui- c3*usq));
#endif
                for (size_t qq = 0; qq < Q; ++qq)
                    ff[qq] += *im++ * (feq - f[q][i]);
            }
            for (size_t q = 0; q < g.S.Q; ++q)
                f[q][i] += ff[q];
        }
    }
}

} // namespace LBM

} // namespace DRL
