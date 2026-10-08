#include "colorfluid.h"
#include "fn.h"
#include <iomanip>
#include <sstream>

namespace DRL { namespace LBM {
class NumericalAudit {
    static void array(std::ostream& os, const va_val_t& values) {
        os << '[';
        for (size_t i=0;i<values.size();++i) {if(i) os<<',';os<<values[i];}
        os << ']';
    }
public:
    static void dump(ColorFluid& s,std::ostream& os) {
        const size_t d=s.g.S.D,q=s.g.S.Q;
        va_val_t matrix=s.getMRTM();
        va_val_t inverse=ALG::InvOrthoRow(matrix,q,q);
        vava_val_t velocity(d);
        va_val_t force(0.0,d),sampleVelocity(0.0,d),directions(0.0,q*d);
        for(size_t k=0;k<d;++k) {
            sampleVelocity[k]=k%2 ? -.013 : .021;
            velocity[k]=va_val_t(sampleVelocity[k],s.g.LG.N);
            force[k]=(k+1)*1e-6;
            for(size_t i=0;i<q;++i) directions[i*d+k]=s.g.S.Qi[i][k];
        }
        os<<std::setprecision(17)<<"MOMENT_AUDIT {\"dimensions\":"<<d<<",\"q\":"<<q;
        os<<",\"directions\":";array(os,directions);
        os<<",\"weights\":";array(os,s.wi);
        os<<",\"matrix\":";array(os,matrix);
        os<<",\"inverse\":";array(os,inverse);
        os<<",\"rates\":";array(os,s.getMRTS(1.2));
        os<<",\"phase0\":";array(os,s.ai[0]);
        os<<",\"phase1\":";array(os,s.ai[1]);
        os<<",\"perturbation_weights\":";array(os,s.getBi());
        os<<",\"velocity\":";array(os,sampleVelocity);
        os<<",\"force\":";array(os,force);
        os<<",\"force_source\":";array(os,s.forceSource(0,velocity,force));
        // A translating phase field can have constant saturation and velocity.
        va_val_t indicator(0.0,s.g.LG.N);
        bool stoppedWhileMoving=false;
        std::ostringstream log;
        for(size_t frame=0;frame<8;++frame) {
            for(size_t i=0;i<indicator.size();++i)
                indicator[i]=std::sin(2.0*3.14159265358979323846*(double(i%8)-.1*frame)/8.0);
            stoppedWhileMoving = s.updateConvergence(velocity,true,1e-6,log,&indicator) || stoppedWhileMoving;
        }
        bool stoppedAfterSettling=false;
        for(size_t frame=0;frame<3;++frame)
            stoppedAfterSettling=s.updateConvergence(velocity,true,1e-6,log,&indicator);
        os<<",\"stopped_while_phase_moving\":"<<(stoppedWhileMoving ? "true":"false")
          <<",\"stopped_after_three_stable_reports\":"<<(stoppedAfterSettling ? "true":"false");
        os<<"}\n";
    }
};
}}

int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::invalid_argument("Provide a ColorFluid parameter file.");
        DRL::IO::FilePrm p(argv[1]);
        DRL::LBM::ColorFluid s(p);
        if(DRL::mpi().getId()==0) DRL::LBM::NumericalAudit::dump(s,std::cout);
        else {std::ostringstream sink;DRL::LBM::NumericalAudit::dump(s,sink);}
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        MPI_Abort(MPI_COMM_WORLD,1);
        return 1;
    }
}
