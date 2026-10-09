#!/usr/bin/env python3
"""Generate, run and independently measure the LBFlow numerical benchmarks.

Requires Python 3 and NumPy. The executable is supplied explicitly; no build or
download occurs here. Output directories are never silently overwritten.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

import numpy as np
from vti import read_vti


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def spec(name, nx=8, ny=8, nz=4, dimensions=3, model='color', eq='MRT',
         tau=(1., 1.), steps=1000, write=250, label=None, indices=None, **kw):
    result = dict(name=name, size=(nx, ny, nz), dimensions=dimensions,
                  model=model, eq=eq, tau=tau, steps=steps, write=write,
                  restart=steps, density=(1., 1.) if model == 'color' else (1.,),
                  alpha=4/9 if dimensions == 2 else 1/3, sigma=.005, beta=.7,
                  label=label or (lambda x, y, z: 0),
                  indices=indices or ['0 "VD" "N/A" "N/A" "0"'])
    result.update(kw)
    return result


class Study:
    def __init__(self, solver, root, mpiexec="mpiexec", numproc_flag="-n", launch_prefix=()):
        self.solver = Path(solver).resolve()
        self.root = Path(root).resolve()
        # Keep every run tied to an immutable executable even if another build occurs.
        data=self.solver.read_bytes()
        self.solver_sha=hashlib.sha256(data).hexdigest()
        executable=self.root/'.bin'/self.solver_sha/self.solver.name
        executable.parent.mkdir(parents=True,exist_ok=True)
        if not executable.exists():
            executable.write_bytes(data)
            executable.chmod(self.solver.stat().st_mode)
        self.solver=executable
        self.env = os.environ.copy()
        self.mpiexec = mpiexec
        self.numproc_flag = numproc_flag
        self.launch_prefix = list(launch_prefix)
        for d in ('cases', 'results'): (self.root/d).mkdir(parents=True, exist_ok=True)

    def run(self, s, ranks=1, start=0, checkpoint=None):
        name = s['name'] + f'-{ranks}rank'
        folder = self.root/'cases'/name
        folder.mkdir(parents=True, exist_ok=False)
        nx, ny, nz = s['size']
        raw = folder/'geometry.raw'
        if s.get('geometry_file'):
            source = Path(s['geometry_file'])
            if source.stat().st_size != nx*ny*nz:
                raise ValueError('Raw geometry size differs from the lattice dimensions')
            shutil.copyfile(source, raw)
        else:
            raw.write_bytes(bytes(s['label'](x,y,z) for z in range(nz)
                                  for y in range(ny) for x in range(nx)))
        dims = (ny,nx) if s['dimensions'] == 2 else (nz,ny,nx)
        params = dict(gridFile=f'"{raw}"', dimension=' '.join(map(str,dims)),
                      directions=9 if s['dimensions']==2 else 19, fileType='"uint8"',
                      latticeLen=1, startTime=start, endTime=s['steps'],
                      writeInterval=s['write'], reportInterval=s.get('report_interval',s['write']),
                      restartInterval=s['restart'], stopError=s.get('stop_error',1e-8),
                      outPath=f'"{folder/"field"}"', logPath=f'"{folder}"',
                      reportDir=s.get('report_dir',0), reportSlice=f'1 {nx-2}',
                      fluidModel=f'"{s["model"]}"', equilModel=f'"{s["eq"]}"',
                      rexlationTime=' '.join(format(t,'.17g') for t in s['tau']),
                      latticeSpeed=1, alpha0=format(s['alpha'],'.17g'),
                      density=' '.join(map(str,s['density'])), contactAngle=90,
                      surfaceTension=s['sigma'], interfaceThickness=s['beta'],
                      wallGradient=f'"{s.get("wall_gradient","legacy")}"')
        params['viscosityInterpolation']=f'"{s.get("viscosity_interpolation","arithmetic")}"'
        if 'capillary_model' in s:
            params['capillaryModel']=f'"{s["capillary_model"]}"'
        if 'wall_recoloring' in s:
            params['wallRecoloring']=f'"{s["wall_recoloring"]}"'
        if 'curvature_stencil' in s:
            params['curvatureStencil']=f'"{s["curvature_stencil"]}"'
        params.update(s.get('overrides',{}))
        text = '\n'.join(f'{key} {val};' for key,val in params.items())
        text += '\n'+'\n'.join(s['indices'])+'\n'
        inp = folder/'input.prm'
        inp.write_text(text)
        if checkpoint:
            for suffix in ('.rst','.rst.meta'):
                shutil.copyfile(str(checkpoint)+suffix, folder/f'field.{start:010d}{suffix}')
        command = self.launch_prefix+[str(self.solver), str(inp)]
        if ranks>1: command = [self.mpiexec,self.numproc_flag,str(ranks)]+command
        begun = time.monotonic()
        try:
            p = subprocess.Popen(command,cwd=folder,env=self.env,stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,text=True,start_new_session=True)
            deadline=begun+s.get('timeout',1200)
            while True:
                remaining=deadline-time.monotonic()
                if remaining<=0:raise subprocess.TimeoutExpired(command,s.get('timeout',1200))
                try:
                    transcript,_=p.communicate(timeout=min(30,remaining))
                    break
                except subprocess.TimeoutExpired:
                    frames=sorted(folder.glob('field.*.vti'))
                    print(json.dumps(dict(running=name,seconds=round(time.monotonic()-begun),
                                          last_frame=frames[-1].name if frames else None)),flush=True)
            code=p.returncode
        except subprocess.TimeoutExpired as e:
            os.killpg(p.pid, signal.SIGKILL)
            transcript,_=p.communicate()
            code='timeout'
            if isinstance(transcript,bytes): transcript=transcript.decode(errors='replace')
        (folder/'run.log').write_text(transcript)
        timing_lines = [line[len('LBFlowTiming '):] for line in transcript.splitlines()
                        if line.startswith('LBFlowTiming ')]
        result=dict(name=name,exit=code,seconds=time.monotonic()-begun,
                    solver_sha256=self.solver_sha,input_sha256=sha(inp),geometry_sha256=sha(raw),
                    parameters={k:v for k,v in s.items() if k!='label'},ranks=ranks,start=start,
                    finished='Simulation finished!' in transcript,command=command,
                    output_files=[str(f) for f in sorted(folder.glob('field.*.vti'))],
                    timing=json.loads(timing_lines[0]) if len(timing_lines)==1 else None)
        if checkpoint:
            result['checkpoint_sha256']=sha(folder/f'field.{start:010d}.rst')
            result['checkpoint_metadata_sha256']=sha(folder/f'field.{start:010d}.rst.meta')
        (self.root/'results'/f'{name}.json').write_text(json.dumps(result,indent=2))
        print(json.dumps({k:result[k] for k in ('name','exit','finished','seconds')}),flush=True)
        if code=='timeout': raise RuntimeError(f'{name}: timed out; inspect {folder/"run.log"}')
        if s.get('expect_failure',False):
            if code==0 or result['finished']: raise AssertionError(f'{name}: invalid input succeeded')
            return result
        if code!=0 or not result['finished']: raise RuntimeError(f'{name}: failed; inspect {folder/"run.log"}')
        return result


def fields(result, at=None):
    files=result['output_files']
    f=Path(files[-1]) if at is None else next((Path(f) for f in files if f'.{at:010d}.' in f),None)
    if f is None and at==0 and result['parameters'].get('initial_output'):
        f=Path(result['parameters']['initial_output'])
    if f is None:raise ValueError(f'Missing requested snapshot at {at}')
    (nx,ny,nz), data=read_vti(f)
    arrays={key:np.array(v['values']).reshape(nz,ny,nx,v['components'])
            for key,v in data.items()}
    if not all(np.isfinite(a).all() for a in arrays.values()):
        raise AssertionError(f'{f}: non-finite output')
    return arrays


def geometry_fields(result, at=None):
    a=fields(result,at)
    if 'v' in a:
        rho=a['r0'][...,0]+a['r1'][...,0]
        chi=np.divide(a['r0'][...,0]-a['r1'][...,0],rho,
                      out=np.zeros_like(rho),where=rho!=0)
        return a,rho,chi,a['v']
    return a,a['density'][...,0],None,a['velocity']


def force_cases():
    cases=[]
    for model in ('scsp','color'):
        for ndim in (2,3):
            for eq in ('SRT','MRT'):
                for tau in (.6,1.4):
                    name=f'force-{model}-d{ndim}-{eq}-tau{tau}'
                    cases.append(spec(name,ny=6,nz=1 if ndim==2 else 4,
                        dimensions=ndim,model=model,eq=eq,
                        tau=(.8,tau) if model=='color' else (tau,),steps=100,write=50,
                        label=lambda x,y,z:75,
                        indices=[f'75 "GV" "0.000001" "1" "{1 if model=="color" else 0}"']))
    return cases


def channel_cases():
    cases=[]
    for model,ndim in (('scsp',3),('color',2)):
        for eq in ('SRT','MRT'):
            for tau in (.6,1.4):
                for h in (8,16,32):
                    nu=(tau-.5)/3
                    n=math.ceil(18*h*h/(nu*math.pi**2))
                    steps=((n+499)//500)*500
                    force=1e-6*(8/h)**3
                    cases.append(spec(f'channel-{model}-d{ndim}-{eq}-tau{tau}-H{h}',
                        nx=4,ny=h+2,nz=1 if ndim==2 else 4,dimensions=ndim,
                        model=model,eq=eq,tau=(.8,tau) if model=='color' else (tau,),
                        steps=steps,write=steps//2,force=force,height=h,
                        label=lambda x,y,z,h=h:255 if y in (0,h+1) else 75,
                        indices=[f'75 "GV" "{force:.17g}" "1" "{1 if model=="color" else 0}"',
                                 '255 "SD" "90" "N/A" "N/A"']))
    return cases


def compare(left,right,checkpoint=False):
    a,b=fields(left),fields(right)
    if a.keys()!=b.keys(): raise AssertionError('Field names differ')
    maxima={key:float(np.max(np.abs(a[key]-b[key]))) for key in a}
    if any(value!=0 for value in maxima.values()):
        raise AssertionError(f'Field comparison differs: {maxima}')
    return dict(left=left['name'],right=right['name'],max_absolute_differences=maxima,
                bitwise_fields=all(np.array_equal(a[key].view(np.uint64),b[key].view(np.uint64)) for key in a))


def verification(study,jobs):
    cases=force_cases()+channel_cases()
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        results=list(pool.map(study.run,cases))
    measured=[]
    for r in results:
        s=r['parameters']; a,rho,chi,v=geometry_fields(r)
        if s['name'].startswith('force-'):
            error=float(np.max(np.abs(v[...,0]-1e-4)))
            transverse=float(np.max(np.abs(v[...,1:])))
            if error>1e-12 or transverse>1e-12: raise AssertionError(f'{r["name"]}: force error {error}')
            measured.append(dict(name=r['name'],kind='force',tau=s['tau'][-1],
                                 absolute_error=error,transverse_max=transverse))
        else:
            ux=v[...,0].mean(axis=(0,2))[1:-1]
            y=np.arange(s['height'])+.5
            mu=(s['tau'][-1]-.5)/3
            ideal=s['force']/(2*mu)*y*(s['height']-y)
            l2=float(np.linalg.norm(ux-ideal)/np.linalg.norm(ideal))
            measured.append(dict(name=r['name'],kind='channel',height=s['height'],tau=s['tau'][-1],
                                 relative_l2_error=l2,relative_flow_error=float(abs(ux.sum()-ideal.sum())/ideal.sum()),
                                 measured_profile=ux.tolist(),analytic_profile=ideal.tolist()))
    restart=[]
    force=force_cases()[-1].copy();force.update(name='restart-force-continuous',restart=50)
    drop=spec('restart-drop-continuous',nx=16,ny=16,nz=16,steps=400,write=200,restart=200,
              label=lambda x,y,z:100 if (x-7.5)**2+(y-7.5)**2+(z-7.5)**2<=16 else 0,
              indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"'])
    for s,cut in ((force,50),(drop,200)):
        continuous=study.run(s,1)
        parallel=study.run(dict(s,name=s['name']+'-parallel'),2)
        restart.append(compare(continuous,parallel))
        cp=Path(continuous['output_files'][-1]).parent/f'field.{cut:010d}'
        for ranks in (1,2):
            resumed=study.run(dict(s,name=s['name'].replace('continuous','resumed')+f'-to{ranks}'),ranks,cut,cp)
            restart.append(compare(continuous,resumed))
    invalid=[]
    base=force_cases()[-1]
    for name,overrides in [('zero-write',{'writeInterval':0}),('invalid-tau',{'rexlationTime':'0.8 0.5'}),
                           ('invalid-alpha',{'alpha0':1.1}),('invalid-speed',{'latticeSpeed':2})]:
        invalid.append(study.run(dict(base,name='guard-'+name,overrides=overrides,expect_failure=True)))
    summary=dict(measurements=measured,restart_and_mpi=restart,
                 invalid_input_failures=[dict(name=r['name'],exit=r['exit']) for r in invalid])
    (study.root/'verification.json').write_text(json.dumps(summary,indent=2))
    return summary


def capillary_cases():
    cases=[]
    for ndim,radii in ((2,(12,20,28)),(3,(8,12,16))):
        for radius in radii:
            n=3*radius+8;center=(n-1)/2
            label=lambda x,y,z,c=center,r=radius,d=ndim:100 if ((x-c)**2+(y-c)**2+
                        ((z-c)**2 if d==3 else 0))<=r*r else 0
            cases.append(spec(f'laplace-d{ndim}-R{radius}',nx=n,ny=n,nz=n if ndim==3 else 1,
                dimensions=ndim,steps=4000,write=1000,radius=radius,
                label=label,indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"']))
    for tau in ((.6,1.4),(1.4,.6)):
        base=cases[1].copy();base.update(name=f'laplace-d2-R20-tau{tau[0]}-{tau[1]}',tau=tau)
        cases.append(base)
    for beta in (.35,.7,.9):
        cases.append(spec(f'planar-d2-beta{beta}',nx=128,ny=8,nz=1,dimensions=2,
            steps=4000,write=1000,beta=beta,label=lambda x,y,z:0 if x<64 else 100,
            indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"']))
    return cases


def linear_crossings(profile,level):
    return [float(i+(level-profile[i])/(profile[i+1]-profile[i])) for i in range(len(profile)-1)
            if (profile[i]-level)*(profile[i+1]-level)<0]


def interface_points(chi,level=0.,fluid_mask=None):
    """Subgrid intersections on every Cartesian edge, in x,y[,z] order."""
    ndim=2 if chi.shape[0]==1 else 3
    points=[]
    for axis in range(3-ndim,3):
        left=[slice(None)]*3;right=left.copy()
        left[axis]=slice(None,-1);right[axis]=slice(1,None)
        a,b=chi[tuple(left)],chi[tuple(right)]
        crossed=(a-level)*(b-level)<0
        if fluid_mask is not None:crossed &= fluid_mask[tuple(left)] & fluid_mask[tuple(right)]
        coordinates=np.column_stack(np.nonzero(crossed)).astype(float)
        if not len(coordinates):continue
        coordinates[:,axis]+=(level-a[crossed])/(b[crossed]-a[crossed])
        points.append(coordinates[:,::-1][:,:ndim])
    return np.concatenate(points) if points else np.empty((0,ndim))


def fit_sphere(points):
    """Orthogonal-distance fit, with algebraic initialization and Gauss-Newton."""
    if len(points)<points.shape[1]+2:raise ValueError('Insufficient contour points')
    design=np.column_stack((2*points,np.ones(len(points))))
    solution=np.linalg.lstsq(design,(points*points).sum(axis=1),rcond=None)[0]
    center=solution[:-1]
    radius=math.sqrt(max(0.,solution[-1]+center@center))
    for _ in range(30):
        delta=points-center
        distances=np.linalg.norm(delta,axis=1)
        residual=distances-radius
        jac=np.column_stack((-delta/distances[:,None],-np.ones(len(points))))
        update=np.linalg.lstsq(jac,-residual,rcond=None)[0]
        center+=update[:-1];radius+=update[-1]
        if np.linalg.norm(update)<1e-11:break
    residual=np.linalg.norm(points-center,axis=1)-radius
    return center,radius,float(np.sqrt(np.mean(residual*residual)))


def measure_capillary(r):
    s=r['parameters'];a,rho,chi,v=geometry_fields(r)
    initial=fields(r,0)
    masses={key:float(a[key].sum()) for key in ('r0','r1')}
    mass_drift={key:float(abs(a[key].sum()-initial[key].sum())/initial[key].sum()) for key in masses}
    summary=dict(name=r['name'],mass_drift=mass_drift,
                 max_spurious_speed=float(np.linalg.norm(v,axis=-1).max()),
                 beta=s['beta'],tau=s['tau'],dimensions=s['dimensions'])
    if s['name'].startswith('planar'):
        profile=chi[0,chi.shape[1]//2]
        low=linear_crossings(profile,-.9);high=linear_crossings(profile,.9)
        middle=linear_crossings(profile,0)[0]
        a=min(low,key=lambda x:abs(x-middle));b=min(high,key=lambda x:abs(x-middle))
        summary.update(width_09=abs(a-b),profile=profile.tolist())
    else:
        nz,ny,nx=chi.shape
        z,y,x=np.indices(chi.shape)
        d2=(x-(nx-1)/2)**2+(y-(ny-1)/2)**2
        if s['dimensions']==3:d2=d2+(z-(nz-1)/2)**2
        center,radius,fit_rms=fit_sphere(interface_points(chi))
        inside=(chi<-.99)&(np.sqrt(d2)<radius-3)
        outside=(chi>.99)&(np.sqrt(d2)>radius+3)
        if not inside.any() or not outside.any():raise AssertionError('No bulk sampling region')
        pressure=rho/3 # All Laplace cases have equal rho and standard pressure weights.
        dp=float(pressure[inside].mean()-pressure[outside].mean())
        analytic=(s['dimensions']-1)*s['sigma']/radius
        previous=s['steps']-s['write']
        prev=geometry_fields(r,previous)[1]/3
        old_dp=float(prev[inside].mean()-prev[outside].mean())
        summary.update(measured_radius=radius,pressure_jump=dp,analytic_jump=analytic,
                       fitted_center=center.tolist(),interface_fit_rms=fit_rms,
                       relative_laplace_error=abs(dp-analytic)/analytic,
                       effective_sigma=dp*radius/(s['dimensions']-1),
                       late_pressure_relative_change=abs(dp-old_dp)/max(abs(dp),1e-30),
                       bulk_sample_counts=[int(inside.sum()),int(outside.sum())])
    return summary


def capillary(study,jobs):
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        results=list(pool.map(study.run,capillary_cases()))
    measured=[measure_capillary(r) for r in results]
    summary=dict(measurements=measured)
    (study.root/'capillary.json').write_text(json.dumps(summary,indent=2))
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('suite',choices=('verification','capillary'))
    parser.add_argument('--solver',required=True)
    parser.add_argument('--root',required=True)
    parser.add_argument('--jobs',type=int,default=4)
    args=parser.parse_args()
    study=Study(args.solver,args.root)
    result=globals()[args.suite](study,args.jobs)
    print(json.dumps({k:len(v) for k,v in result.items()}))


if __name__=='__main__':main()
