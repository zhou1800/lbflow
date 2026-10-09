// Guarded flat, pure-collar CSF pressure experiment. Original Zou--He
// populations and closed numerical operators are retained in colorfluid.cpp.
#include "colorfluid.h"
#include <iomanip>
#include <sstream>

namespace DRL { namespace LBM {

void ColorFluid::initPressureReservoirs()
{
    const auto& i2f=g.LG.getI2fMap();
    bool first=true, defined[2]={false,false};
    size_t phase[2]={0,0}; double density[2]={0.,0.};
    for (size_t idx : p.getIdx({IT::PB})) {
        const auto& value=p.getIVal(idx);
        if (!value.Dir || g.S.Qnorm[value.Dir]!=1.0 || value.Value<=0.0)
            throw std::invalid_argument("CSF pressure reservoirs require positive density and an axial inward direction.");
        const size_t axis=g.S.Qaxis[value.Dir];
        if (first) {reservoirAxis=axis;first=false;}
        if (axis!=reservoirAxis)
            throw std::invalid_argument("CSF pressure reservoirs must share one axis.");
        const size_t side=g.S.Qsign[value.Dir]>0 ? 0 : 1;
        if (defined[side] && (phase[side]!=value.Phase || density[side]!=value.Value))
            throw std::invalid_argument("Each CSF pressure face must prescribe one phase and density.");
        defined[side]=true;phase[side]=value.Phase;density[side]=value.Value;
    }
    if (!defined[0] || !defined[1] || g.Sz[reservoirAxis]<8)
        throw std::invalid_argument("CSF pressure experiments need two opposite faces and distinct interior collars.");
    double localCount[2]={0.,0.},globalCount[2]={0.,0.};
    for (size_t i=0;i<g.LG.N;++i) {
        const auto& value=p.getIVal(raw[i]);
        reservoirTypes[i2f[i]]=value.Type;
        const bool fluid=value.Type==IT::VD || value.Type==IT::GV || value.Type==IT::PB;
        if (!fluid) continue;
        if (value.Type==IT::GV && value.Value!=0.0)
            throw std::invalid_argument("Guarded CSF pressure experiments do not combine external body forcing.");
        const size_t coordinate=(i/g.LG.Stride[reservoirAxis])%g.LG.Sz[reservoirAxis]+g.LStart[reservoirAxis];
        const bool end=coordinate==0 || coordinate+1==g.Sz[reservoirAxis];
        if ((value.Type==IT::PB)!=end)
            throw std::invalid_argument("CSF pressure nodes must cover both outer fluid faces and have no interior nodes.");
        if (value.Type==IT::PB) {
            const size_t side=coordinate==0 ? 0 : 1;
            const size_t directionSide=g.S.Qsign[value.Dir]>0 ? 0 : 1;
            if (side!=directionSide)
                throw std::invalid_argument("CSF pressure face direction must point inward.");
            ++localCount[side];
        }
        if (ibI[i] && ibNorm[ibI[i]][reservoirAxis]!=0.0)
            throw std::invalid_argument("CSF pressure experiments require walls parallel to the reservoir axis.");
        if (coordinate<=2) reservoirPhase[i]=phase[0];
        if (coordinate+3>=g.Sz[reservoirAxis]) reservoirPhase[i]=phase[1];
    }
    mpi().AllReduce(localCount,globalCount,2,MPI_DOUBLE,MPI_SUM);
    if (globalCount[0]==0.0 || globalCount[0]!=globalCount[1])
        throw std::invalid_argument("CSF pressure experiments require two complete equal-area fluid faces.");
    syncBorder(&reservoirTypes,0);
}

void ColorFluid::checkPressureCollar(bool outgoing) const
{
    for (size_t i=0;i<g.LG.N;++i) {
        if (reservoirPhase[i]>1) continue;
        double rho0=0.,rho1=0.;
        if (outgoing) {
            for (size_t q=0;q<g.S.Q;++q) {rho0+=f[0][q][i];rho1+=f[1][q][i];}
        } else {rho0=dens[0][i];rho1=dens[1][i];}
        if (!std::isfinite(rho0+rho1) || rho0+rho1<=0.0 ||
            std::abs((reservoirPhase[i]==0 ? rho1 : rho0)/(rho0+rho1))>1e-9)
            throw std::runtime_error("CSF pressure collar lost its prescribed pure phase; interface too close to reservoir.");
        if (!outgoing && reservoirTypes[g.LG.getI2fMap()[i]]!=IT::PB) {
            const size_t coordinate=(i/g.LG.Stride[reservoirAxis])%g.LG.Sz[reservoirAxis]+g.LStart[reservoirAxis];
            if (coordinate==1 || coordinate+2==g.Sz[reservoirAxis]) {
                double norm=0.;colorGradient(i,norm);
                if (norm>1e-9)
                    throw std::runtime_error("CSF pressure reconstruction requires a force-free first interior strip.");
            }
        }
    }
}

void ColorFluid::initializeBoundaryBudget()
{
    checkPressureCollar(true);
    for (size_t ph=0;ph<2;++ph) {
        double compensation=0.;
        for (size_t idx : p.getIdx({IT::VD,IT::GV}))
            for (size_t i : gIdx[idx])
                for (size_t q=0;q<g.S.Q;++q) {
                    const double increment=f[ph][q][i]-compensation;
                    const double sum=boundaryInitialMass[ph]+increment;
                    compensation=(sum-boundaryInitialMass[ph])-increment;
                    boundaryInitialMass[ph]=sum;
                }
    }
}

void ColorFluid::accumulateBoundaryBudget()
{
    const auto& i2f=g.LG.getI2fMap();const auto& adj=g.LG.FG.getAdjI();
    double net[2]={0.,0.};
    for (size_t idx : p.getIdx({IT::VD,IT::GV}))
        for (size_t i : gIdx[idx])
            for (size_t q=1;q<g.S.Q;++q) {
                const size_t source=adj[i2f[i]][g.S.Qop[q]],target=adj[i2f[i]][q];
                if (reservoirTypes[source]==IT::PB) {
                    const double zero=wallOutgoing0[q][source];
                    net[0]+=zero;net[1]+=wallOutgoingTotal[q][source]-zero;
                }
                if (reservoirTypes[target]==IT::PB)
                    for (size_t ph=0;ph<2;++ph) net[ph]-=f[ph][q][i];
            }
    for (size_t ph=0;ph<2;++ph) {
        const double increment=net[ph]-boundaryCompensation[ph];
        const double sum=boundaryNetMass[ph]+increment;
        boundaryCompensation[ph]=(sum-boundaryNetMass[ph])-increment;
        boundaryNetMass[ph]=sum;
    }
}

void ColorFluid::reportBoundaryBudget(size_t t,std::ostream& os) const
{
    double local[6]={boundaryInitialMass[0],boundaryInitialMass[1],boundaryNetMass[0],boundaryNetMass[1],0.,0.};
    for (size_t idx : p.getIdx({IT::VD,IT::GV})) {
        const auto found=gIdx.find(idx);
        if (found==gIdx.end()) continue;
        for (size_t i : found->second) {local[4]+=dens[0][i];local[5]+=dens[1][i];}
    }
    double global[6]={};mpi().AllReduce(local,global,6,MPI_DOUBLE,MPI_SUM);
    const double residual[2]={global[4]-global[0]-global[2],global[5]-global[1]-global[3]};
    if (std::abs(residual[0])+std::abs(residual[1])>1e-9*std::max(1.0,global[0]+global[1]))
        throw std::runtime_error("CSF pressure interior mass does not balance the recorded reservoir links.");
    if (mpi().getId()==0) {
        std::ostringstream line;line<<std::setprecision(17)
            <<"LBFlowBoundaryBudget {\"step\":"<<t<<",\"reference_step\":"<<boundaryBudgetReference
            <<",\"initial_mass\":["<<global[0]<<','<<global[1]
            <<"],\"integrated_net_mass\":["<<global[2]<<','<<global[3]
            <<"],\"current_mass\":["<<global[4]<<','<<global[5]
            <<"],\"residual\":["<<residual[0]<<','<<residual[1]<<"]}";
        os<<line.str()<<'\n';
    }
}

}}
