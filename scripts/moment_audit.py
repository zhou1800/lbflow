#!/usr/bin/env python3
"""Check analytic moments against coefficients emitted by the real C++ solver."""
import argparse,json,subprocess
from pathlib import Path
import numpy as np
from benchmarks import Study,force_cases


def audit(data):
    q,d=data['q'],data['dimensions']
    e=np.array(data['directions']).reshape(q,d)
    w=np.array(data['weights']);M=np.array(data['matrix']).reshape(q,q)
    inv=np.array(data['inverse']).reshape(q,q);rates=np.array(data['rates'])
    B=np.array(data['perturbation_weights']);u=np.array(data['velocity']);F=np.array(data['force'])
    G=np.array(data['force_source']);I=np.eye(d)
    errors={}
    if data.get('stopped_while_phase_moving',False):raise AssertionError('Moving phase falsely converged')
    if not data.get('stopped_after_three_stable_reports',True):raise AssertionError('Settled phase did not converge')
    def check(name,actual,expected,tolerance=3e-14):
        error=float(np.max(np.abs(np.asarray(actual)-np.asarray(expected))))
        errors[name]=error
        if error>tolerance:raise AssertionError(f'D{d}Q{q}: {name}: {error} > {tolerance}')
    check('inverse',M@inv,np.eye(q))
    check('weight_mass',w.sum(),1)
    check('weight_momentum',w@e,np.zeros(d))
    check('weight_stress',np.einsum('q,qa,qb->ab',w,e,e),I/3)
    check('source_mass',G.sum(),0)
    check('source_force',G@e,F)
    check('source_stress',np.einsum('q,qa,qb->ab',G,e,e),np.outer(u,F)+np.outer(F,u))
    conserved=[0]
    for axis in range(d):
        rows=[r for r in range(q) if np.array_equal(M[r],e[:,axis])]
        if len(rows)!=1:raise AssertionError('Momentum row cannot be identified')
        conserved+=rows
    check('conserved_rates',rates[conserved],np.zeros(len(conserved)))
    check('B_mass',B.sum(),1/3)
    check('B_momentum',B@e,np.zeros(d))
    check('B_stress',np.einsum('q,qa,qb->ab',B,e,e),I/3)
    phi0=np.array(data['phase0']);phi1=np.array(data['phase1'])
    rho0,rho1=.8,1.1;rho=rho0+rho1
    eq0=rho0*phi0+rho1*phi1
    eu=e@u;eq=eq0+rho*w*(3*eu+4.5*eu**2-1.5*(u@u))
    cp0=(.6 if d==2 else .5)*(1-phi0[0]);cp1=(.6 if d==2 else .5)*(1-phi1[0])
    check('equilibrium_mass',eq.sum(),rho)
    check('equilibrium_momentum',eq@e,rho*u)
    check('equilibrium_stress',np.einsum('q,qa,qb->ab',eq,e,e),
          (rho0*cp0+rho1*cp1)*I+rho*np.outer(u,u))
    coupled=set()
    rng=np.random.default_rng(1800)
    for _ in range(32):
        n=rng.normal(size=d);n/=np.linalg.norm(n)
        sigma,gradient,tau=.005,.17,1.2
        delta=2.25*sigma*gradient/tau*(w*(e@n)**2-B)
        check('perturbation_mass',delta.sum(),0)
        check('perturbation_momentum',delta@e,np.zeros(d))
        check('perturbation_stress',np.einsum('q,qa,qb->ab',delta,e,e),
              sigma*gradient/(2*tau)*(np.outer(n,n)-I))
        transformed=M@delta
        rows=np.flatnonzero(abs(transformed)>1e-14)
        coupled.update(map(int,rows))
        check('capillary_relaxation_rates',rates[rows],np.full(len(rows),1/tau))
        lengths=np.linalg.norm(e,axis=1)
        cosine=np.divide(e@n,lengths,out=np.zeros(q),where=lengths!=0)
        R=.7*rho0*rho1/rho**2*cosine*eq0
        f0=rho0/rho*(eq+delta)+R;f1=eq+delta-f0
        check('recoloring_phase0_mass',f0.sum(),rho0)
        check('recoloring_phase1_mass',f1.sum(),rho1)
        check('recoloring_total_momentum',(f0+f1)@e,rho*u)
        # Independent linear-field reconstruction at a flat wall.
        fluid=e[:,-2]>=0 # vertical coordinate is y in the internal reversed order.
        A=np.einsum('q,qa,qb->ab',w[fluid],e[fluid],e[fluid])
        b=np.einsum('q,qa,q->a',w[fluid],e[fluid],e[fluid]@n)
        check('wall_linear_reconstruction',np.linalg.solve(A,b),n)
    return dict(dimensions=d,q=q,errors=errors,capillary_coupled_rows=sorted(coupled),
                conserved_rows=conserved,maximum_error=max(errors.values()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--driver',required=True);p.add_argument('--root',required=True)
    args=p.parse_args();study=Study(args.driver,args.root)
    records=[]
    for ndim in (2,3):
        s=next(s for s in force_cases() if s['model']=='color' and s['dimensions']==ndim)
        s['name']=f'moments-d{ndim}'
        # Study's run contract is an LBSim simulation; the audit driver emits JSON instead.
        folder=study.root/'cases'/s['name'];folder.mkdir(parents=True,exist_ok=False)
        nx,ny,nz=s['size'];raw=folder/'geometry.raw';raw.write_bytes(bytes([75])*(nx*ny*nz))
        dims=(ny,nx) if ndim==2 else (nz,ny,nx)
        text=f'''gridFile "{raw}";
dimension {' '.join(map(str,dims))};
directions {9 if ndim==2 else 19};
fileType "uint8";
latticeLen 1;
startTime 0;
endTime 100;
writeInterval 50;
reportInterval 50;
restartInterval 100;
stopError 0.000001;
outPath "{folder/'field'}";
logPath "{folder}";
reportDir 0;
reportSlice 1 6;
fluidModel "color";
equilModel "MRT";
rexlationTime 1 1;
latticeSpeed 1;
alpha0 {4/9 if ndim==2 else 1/3:.17g};
density 1 1;
contactAngle 90;
surfaceTension 0.005;
interfaceThickness 0.7;
75 "GV" "0.000001" "1" "1"
'''
        inp=folder/'input.prm';inp.write_text(text)
        r=subprocess.run([str(study.solver),str(inp)],env=study.env,text=True,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
        (folder/'operators.log').write_text(r.stdout)
        if r.returncode:raise RuntimeError(r.stdout)
        line=next(line for line in r.stdout.splitlines() if line.startswith('MOMENT_AUDIT '))
        data=json.loads(line.removeprefix('MOMENT_AUDIT '));record=audit(data)
        (folder/'operators.json').write_text(json.dumps(data,indent=2));records.append(record)
    (study.root/'moments.json').write_text(json.dumps(records,indent=2))
    print(json.dumps(records,indent=2))


if __name__=='__main__':main()
