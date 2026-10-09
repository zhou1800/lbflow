/*! \file colorfluid.cpp
 * \author Yu Zhou
 */
#include "colorfluid.h"
#include "logger.h"
#include "fileraw.h"
#include "fn.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace DRL {

namespace LBM {

using namespace DRL::HLP;
using namespace DRL::ALG;
using namespace std;

ColorFluid::ColorFluid(const FilePrm &prm)
    : Simulator(prm)
    , renormalizedWall(p.getWallGradient()=="renormalized")
    , tanhWall(p.getWallGradient()=="tanh")
    , capillaryForce(p.getCapillaryModel()=="csf" ? g.S.D : 0)
    , capillaryGradient(capillaryForce.size())
    , interfaceNormal(capillaryForce.size())
    , SOLID_DENS_N(86)
    , f{vava_val_t(g.S.Q), vava_val_t(g.S.Q)}
    , fOld{vava_val_t(g.S.Q), vava_val_t(g.S.Q)}
    , v(g.S.D)
    , dens{va_val_t(0.0, g.LG.N), va_val_t(0.0, g.LG.N)}
    , densN(SOLID_DENS_N, g.LG.FG.N)
    , densT(0.0, g.LG.N)
    , tau(0.0, g.LG.N)
    , TibNorm(ibNorm.size()), TibNormInv(ibNorm.size())
    , wallGradientInv(ibNorm.size())
    , ai{va_val_t(0.0, g.S.Q), va_val_t(0.0, g.S.Q)}
{
    if (g.S.Q != 9 && g.S.Q != 19)
        throw invalid_argument("Invalid DmQn model for SCSP.");

    for (size_t ph = 0; ph < NP; ++ph)
        for (size_t q = 0; q < g.S.Q; ++q){
            f[ph][q] = va_val_t(0.0, g.LG.N);
            fOld[ph][q] = va_val_t(0.0, g.LG.N);
        }

    for (size_t d = 0; d < g.S.D; ++d)
        v[d]  = va_val_t(0.0, g.LG.N);

    for (size_t d=0; d<capillaryForce.size(); ++d) {
        capillaryForce[d] = capillaryGradient[d] = va_val_t(0.0,g.LG.N);
        interfaceNormal[d] = va_val_t(0.0,g.LG.FG.N);
    }
    initPhi();
    calcTibNorm();
    calcWallGradientInverse();

    output.AddScalar("r0", &dens[0]);
    output.AddScalar("r1", &dens[1]);
    for (size_t d = g.S.D; d != 0; --d)
        output.AddVector("v", &v[d-1]);
    output.Init(g.S.D, g.N, g.Sz, g.Stride, g.LG.Sz, g.LStart);
}

void ColorFluid::Run(std::ostream &logProc)
{
    const double runStart = MPI_Wtime();
    double seconds[TIMING_CATEGORIES] = {};
    size_t completedSteps = 0;
    Logger procLog(logProc);
    procLog.Reset(true);

    procLog.Log("Initializing...", false);
    if (p.getTStart() != 0) {
        const double begin = MPI_Wtime();
        loadRst(p.getTStart());
        seconds[CHECKPOINT_TIME] += MPI_Wtime() - begin;
    }
    else {
        double begin = MPI_Wtime();
        initialize();
        // Construct the t=0 incoming state, then use the normal collision
        // kernel to obtain outgoing populations whose next action is stream.
        calcMacroVar(0);
        if (capillaryForce.size()) initializeCapillaryForce();
        bodyForce(0);
        seconds[INIT_TIME] += MPI_Wtime() - begin;
        begin = MPI_Wtime();
        output.Save(p.getTStart());
        seconds[FIELD_TIME] += MPI_Wtime() - begin;
        begin = MPI_Wtime();
        collision();
        seconds[INIT_TIME] += MPI_Wtime() - begin;
    }
    procLog.Log("Done!\t");

    for (size_t t = p.getTStart() + 1; t <= p.getTEnd(); ++t) {
        bool converged = false;
        double begin = MPI_Wtime();
        // 1. Streaming.
        streaming(&(f[0]), &(fOld[0]),t);
        streaming(&(f[1]), &(fOld[1]),t);

        // 2. Apply doundary condition.
        bounceBack(&f[0], fOld[0], ibS);
        bounceBack(&f[1], fOld[1], ibS);

        // Halfpass film boundary condition.
        bounceBack(&f[0], fOld[0], ibF[1]);
        bounceBack(&f[1], fOld[1], ibF[0]);
        clearHalfFilm();
        //axisymmetric(&f[0]);
        //axisymmetric(&f[1]);
        zouHeBoundary();

        // 3. Calculate macroscopic variables.
        calcMacroVar(t);
        bodyForce(t);
        seconds[NUMERICAL_TIME] += MPI_Wtime() - begin;
        begin = MPI_Wtime();

        // 4. Saves macro variable.
        if (t % p.getTStep() == 0) {
            stringstream ss;
            ss << "Step " << t << " - ";
            procLog.Log(ss.str());
            output.Save(t);
            converged = report(logProc, dens[0], dens[1], v);
        }
        else if (t % p.getTRpt() == 0){
            stringstream ss;
            ss << "Step " << t << " - ";
            procLog.Log(ss.str());
            converged = report(logProc, dens[0], dens[1], v);
        }

        // 5. Collosion.
        seconds[FIELD_TIME] += MPI_Wtime() - begin;
        begin = MPI_Wtime();
        collision();
        seconds[NUMERICAL_TIME] += MPI_Wtime() - begin;

        // 6. Save f.
        if (converged || t % p.getTRst() == 0) {
            begin = MPI_Wtime();
            //procLog.Pause();
            saveRst(t);
            seconds[CHECKPOINT_TIME] += MPI_Wtime() - begin;
            //procLog.Continue();
        }
        ++completedSteps;
        if (converged) break;
    }
    writeTiming(logProc, seconds, MPI_Wtime() - runStart, completedSteps);
}

void ColorFluid::initialize()
{
    static const va_idx_t& i2f {g.LG.getI2fMap()};
    static const double rho[NP] {p.getDensity()[0], p.getDensity()[1]};
    static va_val_t dist[NP] {ai[0] * rho[0], ai[1] * rho[1]};

    for (const auto& idx : p.getIdx({IT::VD, IT::GV, IT::PB, IT::VB, IT::AB,
                                     IT::HF})) {
        size_t phs = p.getIVal(idx).Phase;
        const va_val_t force = forceForIndex(idx);
        for (size_t i : gIdx[idx]){
            for (size_t q = 0; q < g.S.Q; ++q) {
                f[phs][q][i] = dist[phs][q];
                double ef = 0.0;
                for (size_t d = 0; d < g.S.D; ++d)
                    ef += g.S.Qi[q][d] * force[d];
                f[phs][q][i] -= 1.5 * wi[q] * ef;
            }
            densT[i] = dens[phs][i] = rho[phs];
            densN[i2f[i]] = phs == 0 ? 1 : -1;  // (rho0 - rho1) / (rho0 + rho1)
        }
    }
}

void ColorFluid::calcMacroVar(size_t t)
{
    static const va_idx_t& i2f {g.LG.getI2fMap()};
    static const double relaxTime[NP] {p.getRTime()[0]/2.0,
                                       p.getRTime()[1]/2.0};
    static const bool harmonic = p.getViscosityInterpolation() == "harmonic";
    for (size_t idx : p.getIdx({IT::VD, IT::GV, IT::PB, IT::VB, IT::HF}))
        for (size_t i : gIdx[idx]){
            // 1. calculate density
            dens[0][i] = f[0][0][i];
            dens[1][i] = f[1][0][i];
            for (size_t q = 1; q < g.S.Q; ++q){
                dens[0][i] += f[0][q][i];
                dens[1][i] += f[1][q][i];
            }

            densT[i] = dens[0][i] + dens[1][i];
            if (!std::isfinite(dens[0][i]) || !std::isfinite(dens[1][i]) ||
                !std::isfinite(densT[i]) || densT[i] <= 0.0)
                throw runtime_error("Non-finite or non-positive total fluid density.");
            if (dens[0][i] < -1e-10 * densT[i] || dens[1][i] < -1e-10 * densT[i])
                throw runtime_error("Negative component density outside roundoff tolerance.");
            densN[i2f[i]] = (dens[0][i] - dens[1][i]) / densT[i];
            tau[i] = (1 + densN[i2f[i]]) * relaxTime[0] +
                     (1 - densN[i2f[i]]) * relaxTime[1];
            if (harmonic) {
                const double phi = dens[0][i] / densT[i];
                const double nu0 = (p.getRTime()[0] - 0.5) / 3.0;
                const double nu1 = (p.getRTime()[1] - 0.5) / 3.0;
                tau[i] = 0.5 + 3.0 / (phi / nu0 + (1.0 - phi) / nu1);
            }
            if (!std::isfinite(tau[i]) || tau[i] <= 0.5)
                throw runtime_error("Invalid local relaxation time.");

            // 2. calculate velocity.
            for (size_t d = 0; d < g.S.D; ++d) {
                v[d][i] = 0;
                for (size_t q = 1; q < g.S.Q; ++q)      // Qi[0] = {0, 0, 0}
                    if (g.S.Qi[q][d] == 1)
                        v[d][i] += f[0][q][i] + f[1][q][i];
                    else if (g.S.Qi[q][d] == -1)
                        v[d][i] -= f[0][q][i] + f[1][q][i];
                // compressible for all.
                v[d][i] /= densT[i];
                if (!std::isfinite(v[d][i]))
                    throw runtime_error("Non-finite fluid velocity.");
            }
        }

    syncBorder(&densN, t);
}

void ColorFluid::collision()
{

    // shorts
    static const size_t Q = g.S.Q;
    static const va_idx_t& i2f {g.LG.getI2fMap()};
    static const vava_idx_t& adjI {g.LG.FG.getAdjI()};

    // parameters
    static const double c = p.getLSpeed();
    static const double c1 = 3.0 / (c * c);
    static const double c2 = 4.5 / (c * c * c * c);
    static const double c3 = 1.5 / (c * c);
    static const double ak {2.25 * p.getSTension()};
    static const va_val_t bi = getBi();
    //static constexpr double pi = 3.1415926535897932384626433832795;
    //static const double cAngle = cos(p.getCAngle() / 180.0 * pi) ;
    //static const double sAngle = sin(p.getCAngle() / 180.0 * pi) ;
    static const double beta = p.getIThickness();
    const bool csf = capillaryForce.size() != 0;

    // MRT
    static const EQUILFN equilfn = p.getEquilModel() == "SRT" ? SRT : MRT;
    static const va_val_t M {getMRTM()};
    static const va_val_t MInv = InvOrthoRow(M, Q, Q);
    static const bool uniformTau = p.getRTime()[0] == p.getRTime()[1];
    static const va_val_t uniformRates = getMRTS(p.getRTime()[0]);
    static const va_val_t uniformCollision = MatMulMat(
        MatMulDiag(MInv, uniformRates, Q), M, Q, Q, Q);
    static const va_val_t uniformForce = MatMulMat(
        MatMulDiag(MInv, 1.0 - 0.5 * uniformRates, Q), M, Q, Q, Q);
    for (size_t idx : p.getIdx({IT::VD, IT::GV, IT::VB, IT::PB, IT::HF})) {
        const va_val_t external = forceForIndex(idx);
        const bool externallyForced = inner_product(begin(external),end(external),begin(external),0.0) != 0.0;
        va_val_t localForce(0.0,csf ? g.S.D : 0);
        for (size_t i : gIdx[idx]) {     // array index
            if (csf) {
                localForce = external;
                for (size_t d=0; d<g.S.D; ++d) localForce[d] += capillaryForce[d][i];
            }
            const va_val_t& force = csf ? localForce : external;
            const bool forced = csf ? inner_product(begin(force),end(force),begin(force),0.0) != 0.0 : externallyForced;
            // 2. B-G-K collision.
            double usq = 0;
            for (size_t d = 0; d < g.S.D; ++d)
                usq += v[d][i] * v[d][i];

            va_val_t feq_u0 (g.S.Q);       // save for later use;
            va_val_t fT (g.S.Q);           // f_total = f_r + f_b
            const va_val_t source = forced ? forceSource(i, v, force)
                                           : va_val_t(0.0, Q);

            if (equilfn == SRT){
                for (size_t q = 0; q < g.S.Q; ++q) {
                    fT[q] = f[0][q][i] + f[1][q][i];
                    double ui = 0;
                    for (size_t d = 0; d < g.S.D; ++d)
                        ui += v[d][i] * g.S.Qi[q][d];

                    feq_u0[q] = dens[0][i] * ai[0][q] + dens[1][i] * ai[1][q];
                    double feq = feq_u0[q]
                            + densT[i] * wi[q] * (c1*ui + c2*ui*ui - c3*usq);

                    fT[q] += (feq - fT[q]) / tau[i];
                    if (forced)
                        fT[q] += (1.0 - 0.5 / tau[i]) * source[q];
                }
            }
            else {
                va_val_t residual(0.0, Q);
                for (size_t q = 0; q < g.S.Q; ++q) {
                    fT[q] = f[0][q][i] + f[1][q][i];
                    double ui = 0;
                    for (size_t d = 0; d < g.S.D; ++d)
                        ui += v[d][i] * g.S.Qi[q][d];

                    feq_u0[q] = dens[0][i] * ai[0][q] + dens[1][i] * ai[1][q];
                    double feq = feq_u0[q]
                            + densT[i] * wi[q] * (c1*ui + c2*ui*ui - c3*usq);

                    residual[q] = feq - fT[q];
                }
                va_val_t correction(0.0, Q);
                if (uniformTau) {
                    correction = MatMulVec(uniformCollision, residual, Q, Q);
                    if (forced) correction += MatMulVec(uniformForce, source, Q, Q);
                } else {
                    const va_val_t rates = getMRTS(tau[i]);
                    correction = MatMulVec(M, residual, Q, Q) * rates;
                    if (forced)
                        correction += MatMulVec(M, source, Q, Q) * (1.0 - 0.5 * rates);
                    correction = MatMulVec(MInv, correction, Q, Q);
                }
                for (size_t q = 0; q < g.S.Q; ++q)
                    fT[q] += correction[q];
            }

            double norm = 0.0;
            va_val_t gDens;
            if (csf) {
                gDens = va_val_t(0.0,g.S.D);
                for (size_t d=0; d<g.S.D; ++d) gDens[d] = capillaryGradient[d][i];
                norm = sqrt(inner_product(begin(gDens),end(gDens),begin(gDens),0.0));
            } else {
                gDens = colorGradient(i,norm);
            }

            double dens0T = dens[0][i] / densT[i];
            if (norm > 1e-9){
                double betaDens01TT = beta * dens0T * dens[1][i] / densT[i];
                double AkNorm = ak * norm / tau[i];
                for (size_t q = 0; q < g.S.Q; ++q){
                    double angle = 0.0;
                    for (size_t d = 0; d < g.S.D; ++ d){
                        if (g.S.Qi[q][d] == 1)  angle += gDens[d];
                        else if (g.S.Qi[q][d] == -1) angle -= gDens[d];
                    }
                    angle /= norm;

                    if (!csf) fT[q] += AkNorm * (wi[q]*angle*angle - bi[q]);

            // 4. recolor.
                    if (q != 0) angle /= g.S.Qnorm[q];
                    f[0][q][i] = dens0T * fT[q]
                            + betaDens01TT * angle * feq_u0[q];
                    f[1][q][i] = fT[q] - f[0][q][i];
                }
            }
            else {
                for (size_t q = 0; q < g.S.Q; ++q){
                    f[0][q][i] = dens0T * fT[q];
                    f[1][q][i] = fT[q] - f[0][q][i];
                }
            }
        }
    }
    validatePopulations(f[0]);
    validatePopulations(f[1]);
}

va_val_t ColorFluid::colorGradient(size_t i, double& norm) const
{
    const size_t Q = g.S.Q;
    const auto& i2f = g.LG.getI2fMap();
    const auto& adjI = g.LG.FG.getAdjI();
    // Common color gradient and geometric wall projection.
    va_val_t gDens(0.0, g.S.D);
    norm = 0;
    if (!ibI[i]){
        for (size_t d = 0; d < g.S.D; ++ d){
            for (size_t q : g.S.QrelD[d]){
                if (g.S.Qi[q][d] == 1)
                    gDens[d] += wi[q] * densN[adjI[i2f[i]][q]];
                else// if (g.S.Qi[q][d] == -1)
                    gDens[d] -= wi[q] * densN[adjI[i2f[i]][q]];
            }
        }
        gDens *= 3.0;       // 3.0 / (c * t)
        norm = sqrt(inner_product(
                    begin(gDens), end(gDens), begin(gDens),
                               static_cast<double>(0.0)));
    }
    else{    // solid - liquid boundary.
        for (size_t d = 0; d < g.S.D; ++ d){
            for (size_t q : g.S.QrelD[d]){
                double densN0 = densN[i2f[i]];
                double densNQ = densN[adjI[i2f[i]][q]];
                if (densNQ == SOLID_DENS_N)
                    continue;   // skip solid boundary.
                if (g.S.Qi[q][d] == 1)
                    gDens[d] += wi[q] * (densNQ - densN0);
                else// if (g.S.Qi[q][d] == -1)
                    gDens[d] -= wi[q] * (densNQ - densN0);
            }
        }
        gDens *= 3.0;       // 3.0 / (c * t)
        norm = sqrt(inner_product(
                    begin(gDens), end(gDens), begin(gDens),
                               static_cast<double>(0.0)));

        double cAngle = ibCosCA[ibI[i]];
        double sAngle = ibSinCA[ibI[i]];
        if (renormalizedWall) {
            // gDens currently contains 3*b from the clipped stencil.
            gDens = MatMulVec(wallGradientInv[ibI[i]], gDens,
                              g.S.D, g.S.D) / 3.0;
            norm = sqrt(inner_product(begin(gDens), end(gDens),
                                       begin(gDens), 0.0));
        }
        if (tanhWall) {
            // An equilibrium tanh profile has locally affine atanh(chi).
            // Continue that fitted profile over a virtual full stencil.
            const auto transform = [](double color) {
                return atanh(max(-1.0 + 1e-8, min(1.0 - 1e-8, color)));
            };
            const double center = transform(densN[i2f[i]]);
            va_val_t moment(0.0, g.S.D);
            for (size_t q = 1; q < Q; ++q) {
                const double neighbor = densN[adjI[i2f[i]][q]];
                if (neighbor == SOLID_DENS_N) continue;
                const double difference = transform(neighbor) - center;
                for (size_t d = 0; d < g.S.D; ++d)
                    moment[d] += wi[q] * g.S.Qi[q][d] * difference;
            }
            const va_val_t slope = MatMulVec(wallGradientInv[ibI[i]], moment,
                                            g.S.D, g.S.D);
            gDens = va_val_t(0.0, g.S.D);
            for (size_t q = 1; q < Q; ++q) {
                double position = center;
                for (size_t d = 0; d < g.S.D; ++d)
                    position += g.S.Qi[q][d] * slope[d];
                const double color = tanh(position);
                for (size_t d = 0; d < g.S.D; ++d)
                    gDens[d] += 3.0 * wi[q] * g.S.Qi[q][d] * color;
            }
            norm = sqrt(inner_product(begin(gDens), end(gDens),
                                       begin(gDens), 0.0));
        }
        if (g.S.D == 3) {
            va_val_t gDensNorm = MatMulMat(TibNorm[ibI[i]], gDens,
                                            3, 3, 1);
            va_val_t txInv = RotateTx3D(-gDensNorm[0], gDensNorm[1]);
            va_val_t tInv = MatMulMat(TibNormInv[ibI[i]], txInv,
                    3, 3, 3);
            gDens = {0, norm * sAngle, norm * cAngle};
            gDens = MatMulMat(tInv, gDens, 3, 3, 1);
        }
        else if (g.S.D == 2) {
            va_val_t gDensNorm = MatMulMat(TibNorm[ibI[i]], gDens,
                                            2, 2, 1);
            if (gDensNorm[0] >= 0)
                gDens = {norm * sAngle, norm * cAngle};
            else
                gDens = {-norm * sAngle, norm * cAngle};
            gDens = MatMulMat(TibNormInv[ibI[i]], gDens, 2, 2, 1);
        }
    }
    return gDens;
}

void ColorFluid::calcCapillaryForce(size_t t)
{
    const auto& i2f = g.LG.getI2fMap();
    const auto& adj = g.LG.FG.getAdjI();
    const size_t D = g.S.D;
    for (size_t idx : p.getIdx({IT::VD,IT::GV})) {
        for (size_t i : gIdx[idx]) {
            double norm;
            const va_val_t gradient = colorGradient(i,norm);
            for (size_t d=0; d<D; ++d) {
                capillaryGradient[d][i] = gradient[d];
                interfaceNormal[d][i2f[i]] = norm > 1e-9 ? gradient[d]/norm : 0.0;
            }
        }
    }
    for (size_t d=0; d<D; ++d) syncBorder(&interfaceNormal[d],t);
    for (size_t idx : p.getIdx({IT::VD,IT::GV})) {
        for (size_t i : gIdx[idx]) {
            double curvature = 0.0, norm2 = 0.0;
            for (size_t d=0; d<D; ++d) norm2 += capillaryGradient[d][i]*capillaryGradient[d][i];
            if (norm2 > 1e-18) {
                for (size_t a=0; a<D; ++a) {
                    va_val_t moment(0.0,D);
                    const double center = interfaceNormal[a][i2f[i]];
                    for (size_t q=1; q<g.S.Q; ++q) {
                        if (densN[adj[i2f[i]][q]]==SOLID_DENS_N) continue;
                        const double difference = interfaceNormal[a][adj[i2f[i]][q]]-center;
                        for (size_t b=0; b<D; ++b)
                            moment[b] += wi[q]*g.S.Qi[q][b]*difference;
                    }
                    const va_val_t derivative = ibI[i] ?
                        MatMulVec(wallGradientInv[ibI[i]],moment,D,D) : moment*3.0;
                    for (size_t b=0; b<D; ++b)
                        curvature += (center*interfaceNormal[b][i2f[i]]-(a==b ? 1.0 : 0.0))*derivative[b];
                }
            }
            for (size_t d=0; d<D; ++d) {
                capillaryForce[d][i] = .5*p.getSTension()*curvature*capillaryGradient[d][i];
                if (!isfinite(capillaryForce[d][i])) throw runtime_error("Non-finite capillary force.");
            }
        }
    }
}

void ColorFluid::initializeCapillaryForce()
{
    calcCapillaryForce(0);
    for (size_t idx : p.getIdx({IT::VD,IT::GV})) {
        for (size_t i : gIdx[idx]) {
            for (size_t q=0; q<g.S.Q; ++q) {
                double ef = 0.0;
                for (size_t d=0; d<g.S.D; ++d) ef += g.S.Qi[q][d]*capillaryForce[d][i];
                for (size_t ph=0; ph<NP; ++ph)
                    f[ph][q][i] -= 1.5*wi[q]*ef*dens[ph][i]/densT[i];
            }
        }
    }
    calcMacroVar(0);
}

void ColorFluid::zouHeBoundary()
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
                sumF += f[0][q][i] + f[0][q][i] + f[1][q][i] + f[1][q][i];

            for (size_t q : qNorm)
                sumF += f[0][q][i] + f[1][q][i];

            for (size_t d = 0; d < g.S.D; ++d) {
                if (d != qAxis) v[d][i] = 0;
                else {
                    v[d][i] = velocity;
                    densT[i] = sumF * vel_1;
                }
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
                sumF += f[0][q][i] + f[0][q][i] + f[1][q][i] + f[1][q][i];

            for (size_t q : qNorm)
                sumF += f[0][q][i] + f[1][q][i];

            densT[i] = density;
            for (size_t d = 0; d < g.S.D; ++d) {
                if (d != qAxis) v[d][i] = 0;
                else {
                    sumF /= density; sumF--;
                    v[d][i] = qSign == 1 ? - sumF : sumF;
                }
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

        size_t phs = p.getIVal(idx).Phase;
        size_t phsOp = 1 - phs;

        for (size_t i : gIdx[idx])
            for (size_t q : qRel) {  // for each unknown dist[q]
                // calculate distribution function.
                // 1. dist = dist_opp
                double fT = f[0][g.S.Qop[q]][i] + f[1][g.S.Qop[q]][i];
                // 2. dist += 6 * wi * rho * u_e / c
                fT += qSign * wiT[q] * densT[i] * v[qAxis][i];
                // 3. dist += 1/2 * sum(F_normal)
                double sumF = 0;
                for (size_t qn : qPrep) {
                    for (size_t d = 0; d < g.S.D; ++d) {
                        if (g.S.Qi[qn][d] == 0 || g.S.Qi[q][d] == 0)
                            continue;
                        else if (g.S.Qi[qn][d] == g.S.Qi[q][d])
                            sumF += f[0][qn][i] + f[1][qn][i];
                        else
                            sumF -= f[0][qn][i] + f[1][qn][i];
                    }
                }
                fT -= sumF * 0.5;
                f[phs][q][i] = fT;
                f[phsOp][q][i] = 0;
            }
    }
}

void ColorFluid::bodyForce(size_t t)
{
    if (capillaryForce.size()) {
        calcCapillaryForce(t);
        for (size_t idx : p.getIdx({IT::VD,IT::GV}))
            for (size_t i : gIdx[idx])
                for (size_t d=0; d<g.S.D; ++d) v[d][i] += .5*capillaryForce[d][i]/densT[i];
    }
    for (size_t idx : p.getIdx({IT::GV}))
        for (size_t i : gIdx[idx]) {
            const va_val_t fi = forceForIndex(idx) * 0.5;
            for (size_t d = 0; d < g.S.D; ++d)
                v[d][i] += fi[d] / densT[i];
        }
}

void ColorFluid::clearHalfFilm()
{
    for (size_t idx : p.getIdx({IT::HF})){
        if (p.getIVal(idx).Phase == 0)
            for (size_t i : gIdx[idx])
                for (size_t q = 0; q < g.S.Q; ++q)
                    f[1][q][i] = 0;
        if (p.getIVal(idx).Phase == 1)
            for (size_t i : gIdx[idx])
                for (size_t q = 0; q < g.S.Q; ++q)
                    f[0][q][i] = 0;
    }
}

void ColorFluid::loadRst(std::size_t t)
{
    //static const va_idx_t& i2f {g.LG.getI2fMap()};

    // prepare filename and load data.
    ostringstream os;
    os << p.getOutPath() << "." << setw(10) << setfill('0') << t << ".rst";
    validateRestart(os.str(), t, NP);
    FileRaw<double> file(os.str(), true);

    // load distribution functions.
    for (size_t ph = 0; ph < NP; ++ph)
        for (size_t q = 0; q < g.S.Q; ++q)
            f[ph][q] = file.Load(g, t);

    //va_val_t data = FileRaw<double>(os.str()).Load(g.S.Q * NP, g, t);

    //// 2.organize data
    //for (size_t q = 0,s=0; q < g.S.Q; ++q)
    //    for (size_t ph = 0; ph < NP; ++ph,s += g.LG.N)
    //        f[ph][q] = data[slice(s, g.LG.N, 1)];
}

void ColorFluid::saveRst(std::size_t t)
{
    // prepare filename
    ostringstream os;
    os << p.getOutPath() << "." << setw(10) << setfill('0') << t << ".rst";
    prepareRestart(os.str());
    FileRaw<double> file(os.str(), true);

    // save distribution function.
    for (size_t ph = 0; ph < NP; ++ph)
        for (size_t q = 0; q < g.S.Q; ++q)
            file.Save(g, f[ph][q], t);
    saveRestartMetadata(os.str(), t, NP);

    //// 2.prepare data
    //va_val_t data(g.S.Q * g.LG.N * NP);
    //for (size_t q = 0,s=0; q < g.S.Q; ++q) {
    //    for (size_t ph = 0; ph < NP; ++ph,s += g.LG.N)
    //        data[slice(s, g.LG.N, 1)] = f[ph][q];
    //}

    //// 3.Save
    //FileRaw<double>(os.str()).Save(g.N * g.S.Q * NP, data, t);
}

bool ColorFluid::report(ostream &os, const va_val_t &r0, const va_val_t &r1,
                        const vava_val_t &v)
{
    static const size_t rptDir = p.getReportDir();
    if (rptDir == 0)  // no report
        return false;
    try{
        static const size_t d = g.S.Qaxis[rptDir];  // report axis.
        static const ve_idx_t slices = p.getReportSlice(); // report slices.
        static const auto fl = minmax_element(slices.begin(), slices.end());
        static const va_idx_t lclRange = getElem(d, *fl.first, *fl.second);
        static const vava_idx_t lclSlice = {  // first and last slice.
            getElem(d, *fl.first, *fl.first),
            getElem(d, *fl.second, *fl.second)};
        static const double cntS[] = {      // number of elements in f&l slice.
            getElemCount(lclSlice[0]),
            getElemCount(lclSlice[1])};

        static const auto dens0 = p.getDensity();
        static const double MaxErr = p.getStopError();
        static const double kPar = getKPar(d);
        static const double pressureFactor = g.S.Q == 9 ? 3.0 / 5.0 : 1.0 / 2.0;
        static const double r2p[] = {(1.0 - ai[0][0]) * pressureFactor,
                                      (1.0 - ai[1][0]) * pressureFactor};
        static va_val_t SOld = {0.0, 0.0};    // Saturation of last step.
        static double KOld = 0;                // Permeability of last step.

        // 1. Accumulate total density and velocity of range.
        double rhoR[] = {0, 0};     // Total sum density of range.
        double velR = 0;             // Total sum velocity of range.
        {
            double rhoLR[] = {0, 0};     // Local sum density of range.
            double velLR = 0;             // local sum velocity of range.
            for (auto a : lclRange){
                rhoLR[0] += r0[a];
                rhoLR[1] += r1[a];
                velLR += v[d][a];
            }
            mpi().AllReduce(&rhoLR[0], &rhoR[0], 1, MPI_DOUBLE, MPI_SUM);
            mpi().AllReduce(&rhoLR[1], &rhoR[1], 1, MPI_DOUBLE, MPI_SUM);
            mpi().AllReduce(&velLR,    &velR,    1, MPI_DOUBLE, MPI_SUM);
        }

        // 2. Accumulate density of first and last slice.
        double rhoS[2][2] = {0, 0, 0, 0};   // Total rho[slice][phase] of slice.
        for (size_t i = 0; i < 2; ++i){
            double rhoLS[2] = {0, 0};
            for (auto a : lclSlice[i]){
                rhoLS[0] += r0[a];
                rhoLS[1] += r1[a];
            }
            mpi().AllReduce(&rhoLS[0], &rhoS[i][0], 1, MPI_DOUBLE, MPI_SUM);
            mpi().AllReduce(&rhoLS[1], &rhoS[i][1], 1, MPI_DOUBLE, MPI_SUM);
        }

        // 3. Calculate Sw, So, K.
        double stcR[] = {rhoR[0] / dens0[0], rhoR[1] / dens0[1]};  // substance
        double stcRT = stcR[0] + stcR[1];       // total substance.
        va_val_t S = {stcR[0] / stcRT, stcR[1] / stcRT};   // saturation.
        va_val_t errS = {abs(S[0] - SOld[0]), abs(S[1] - SOld[1])};
        SOld = S;

        double dP = (rhoS[0][0] * r2p[0] + rhoS[0][1] * r2p[1]) / cntS[0] -
                    (rhoS[1][0] * r2p[0] + rhoS[1][1] * r2p[1]) / cntS[1];
        double Q = velR / double(*fl.second - *fl.first + 1);
        if (!std::isfinite(S[0]) || !std::isfinite(S[1]) ||
            !std::isfinite(errS[0]) || !std::isfinite(errS[1]) ||
            !std::isfinite(dP) || !std::isfinite(Q))
            throw runtime_error("Non-finite two-phase report state.");
        const bool hasPressureDrop = abs(dP) > 1e-14;
        const double K = hasPressureDrop ? kPar * Q / dP : 0.0;
        const double errK = hasPressureDrop ? abs(K - KOld) / std::max(abs(K), 1e-30) : 0.0;
        if (!std::isfinite(K) || !std::isfinite(errK))
            throw runtime_error("Non-finite apparent permeability.");
        KOld = K;

        // 4. Report.
        if (mpi().getId() == 0){
            os << setprecision(6) << fixed
               << " - S0 = " << S[0] * 100 << "%, "
               << "dS0 = " << errS[0] * 100 << "%" << endl
               << " - S1 = " << S[1] * 100 << "%, "
               << "dS1 = " << errS[1] * 100 << "%" << endl
               << scientific
               << " - dP = " << dP << "(mu/lu/ts^2), "
               << "Q = " << Q << "(lu^3/ts)" << endl
               << " - K = " << K << "(um^2), "
               << fixed << "errK = " << errK * 100.0 << "%" << endl;
        }
        if (!hasPressureDrop && mpi().getId() == 0)
            os << " - pressure-based permeability unavailable: zero pressure drop" << endl;
        const bool scalarsStable = std::max(errS[0], errS[1]) < MaxErr &&
                (!hasPressureDrop || errK < MaxErr);
        va_val_t indicator(0.0, r0.size());
        for (size_t i = 0; i < r0.size(); ++i)
            if (r0[i] + r1[i] > 0.0) indicator[i] = (r0[i] - r1[i]) / (r0[i] + r1[i]);
        return updateConvergence(v, scalarsStable, MaxErr, os, &indicator);
    }
    catch (exception& e) { os << e.what(); throw; }
    return false;
}

void ColorFluid::calcWallGradientInverse()
{
    if (p.getWallGradient() != "renormalized" && p.getWallGradient() != "tanh" &&
        !capillaryForce.size()) return;
    const auto& i2f = g.LG.getI2fMap();
    const auto& adj = g.LG.FG.getAdjI();
    va_u8_t extended(g.LG.FG.N);
    for (size_t i = 0; i < g.LG.N; ++i) extended[i2f[i]] = raw[i];
    syncBorder(&extended, 0);
    const size_t D = g.S.D;
    for (size_t i = 0; i < g.LG.N; ++i) {
        if (!ibI[i]) continue;
        double augmented[3][6] = {};
        for (size_t q = 1; q < g.S.Q; ++q) {
            const auto type = p.getIVal(extended[adj[i2f[i]][q]]).Type;
            if (type == IT::SD || type == IT::BB) continue;
            for (size_t a = 0; a < D; ++a)
                for (size_t b = 0; b < D; ++b)
                    augmented[a][b] += wi[q] * g.S.Qi[q][a] * g.S.Qi[q][b];
        }
        for (size_t a = 0; a < D; ++a) augmented[a][D+a] = 1.0;
        for (size_t col = 0; col < D; ++col) {
            size_t pivot = col;
            for (size_t row = col+1; row < D; ++row)
                if (abs(augmented[row][col]) > abs(augmented[pivot][col])) pivot = row;
            if (abs(augmented[pivot][col]) < 1e-12)
                throw runtime_error("Rank-deficient fluid wall stencil; reconstructed derivatives are unsupported here.");
            for (size_t k = 0; k < 2*D; ++k) swap(augmented[col][k],augmented[pivot][k]);
            const double diagonal = augmented[col][col];
            for (size_t k = 0; k < 2*D; ++k) augmented[col][k] /= diagonal;
            for (size_t row = 0; row < D; ++row) {
                if (row == col) continue;
                const double scale = augmented[row][col];
                for (size_t k = 0; k < 2*D; ++k) augmented[row][k] -= scale*augmented[col][k];
            }
        }
        va_val_t inverse(0.0,D*D);
        for (size_t a = 0; a < D; ++a)
            for (size_t b = 0; b < D; ++b) inverse[a*D+b] = augmented[a][D+b];
        wallGradientInv[ibI[i]] = inverse;
    }
}

void ColorFluid::calcTibNorm()
{
    for (size_t i = 1; i < ibNorm.size(); ++i){
        if (g.S.D == 3){
            va_val_t tx = RotateTx3D(ibNorm[i][0], ibNorm[i][1]);
            va_val_t txInv = RotateTx3D(-ibNorm[i][0], ibNorm[i][1]);
            double val = sqrt(ibNorm[i][0] * ibNorm[i][0]
                    + ibNorm[i][1] * ibNorm[i][1]);
            va_val_t tz = RotateTz3D(val, ibNorm[i][2]);
            va_val_t tzInv = RotateTz3D(-val, ibNorm[i][2]);
            TibNorm[i] = MatMulMat(tz, tx, 3, 3, 3);
            TibNormInv[i] = MatMulMat(txInv, tzInv, 3, 3, 3);
        }
        else if (g.S.D == 2){
            TibNorm[i] = RotateT2D(ibNorm[i][0], ibNorm[i][1]);
            TibNormInv[i] = RotateT2D(-ibNorm[i][0], ibNorm[i][1]);
        }
    }
}

void ColorFluid::initPhi()
{
    double alpha[NP];
    alpha[0] = p.getAlpha0();
    alpha[1] = 1 - p.getDensity()[0] / p.getDensity()[1] * (1 - alpha[0]);
    for (size_t ph = 0; ph < NP; ++ph) {
        if (g.S.Q == 9) {
            double a = (1 - alpha[ph]) / 5;
            double b = (1 - alpha[ph]) / 20;
            ai[ph] = va_val_t({alpha[ph], a, a, a, a, b, b, b, b});
        }

        else if (g.S.Q == 19) {
            double a = (1 - alpha[ph]) / 12;
            double b = (1 - alpha[ph]) / 24;
            ai[ph] = va_val_t({alpha[ph], a, a, a, a, a, a,
                              b, b, b, b, b, b, b, b, b, b, b, b});
        }
    }
}

va_val_t ColorFluid::getBi()
{
    constexpr double x = 2;
    static const double c2 = p.getLSpeed() * p.getLSpeed();
    static const double b0 = -(2.0 + 2.0*x) / (3.0*x + 12.0) * c2;
    static const double b1 = x / (6.0*x + 24.0) * c2;
    static const double b7 = 1.0 / (6.0*x + 24.0) * c2;

    if (g.S.Q == 9)
        return va_val_t({-4.0/27,
                        2.0/27,  2.0/27,  2.0/27,  2.0/27,
                        5.0/108, 5.0/108, 5.0/108, 5.0/108});
    else if (g.S.Q == 19)
        return va_val_t({b0, b1, b1, b1, b1, b1, b1,
                          b7, b7, b7, b7, b7, b7, b7, b7, b7, b7, b7, b7});
    throw invalid_argument("No capillary perturbation weights for this lattice.");
}

} // LBM

} // DRL
