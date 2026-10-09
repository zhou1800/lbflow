#!/usr/bin/env python3
"""Production neutral bridges and transverse phase controls in common units."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np

from benchmarks import Study, geometry_fields, sha, spec
from front_mobility import match_beta, WIDTH_REFERENCE
from phase_recurrence import bulk_profile, coefficients, VELOCITIES, WEIGHTS


def gradient(phi):
    color=2*phi-1;extended=np.vstack((color[0],color,color[-1]));g=np.zeros((2,)+phi.shape)
    for (dx,dy),w in zip(VELOCITIES,WEIGHTS):
        value=np.roll(extended,-dx,axis=1)[1+dy:1+dy+phi.shape[0]]
        g[0]+=3*w*dx*value;g[1]+=3*w*dy*value
    norm=np.linalg.norm(g,axis=0)
    return np.divide(g,norm[None],out=np.zeros_like(g),where=norm[None]>1e-9)


def seed(path,h,length=4,amplitude=0.,mode=1,width=WIDTH_REFERENCE):
    if not np.isfinite(width) or width<=0:raise ValueError('Use a finite positive physical width')
    nx=length*h;beta=match_beta(width*h/16)
    ratio=coefficients()['bulk_segregation']/coefficients()['bulk_diffusion']
    profile=bulk_profile(nx,beta/ratio);dx=16/h
    derivative=(np.roll(profile,-1)-np.roll(profile,1))/(2*dx)
    y=(np.arange(h)+.5)*dx
    phi=np.broadcast_to(profile,(h,nx)).copy()-amplitude*np.cos(mode*np.pi*y[:,None]/16)*derivative[None]
    if phi.min()<0 or phi.max()>1:raise ValueError('Perturbation violates phase fractions')
    normal=gradient(phi);pop=np.zeros((2,9,h+2,nx))
    for q,((ex,ey),w) in enumerate(zip(VELOCITIES,WEIGHTS)):
        phase=w*(phi+beta*phi*(1-phi)*(ex*normal[0]+ey*normal[1]))
        pop[0,q,1:-1]=phase;pop[1,q,1:-1]=w-phase
    Path(str(path)+'.rst').write_bytes(pop.astype('=f8').tobytes())
    Path(str(path)+'.rst.meta').write_text(f'DRLAB_OUTGOING_V1\n1 2 2 9\n{h+2} {nx} \nNEXT_STREAM 2\nINITIAL_COLLISION_INCLUDED 2\nCONVERGENCE_HISTORY_RESET_ON_LOAD\n')
    np.savez_compressed(str(path)+'.phase.npz',phi=phi,profile=profile,amplitude=amplitude,dx=dx)
    return beta,phi


def config(h,force=1e-7,sigma=.005,length=4,end=8000.,amplitude=0.,mode=1):
    dx=16/h;nx=length*h;force_lattice=force*dx**3
    return spec(f'periodic-H{h}-L{length}-F{force:g}-sigma{sigma:g}-amp{amplitude:g}-mode{mode}',
        nx=nx,ny=h+2,nz=1,dimensions=2,tau=(1.,1.),sigma=sigma*dx,
        capillary_model='csf',wall_gradient='ghost',phase_boundary='flux',
        recoloring_model='projection',mrt_spectrum='legacy',height=h,length_factor=length,
        physical_force=force,physical_sigma=sigma,physical_end=end,physical_amplitude=amplitude,
        transverse_mode=mode,dx=dx,steps=1+round(end/dx**2),write=max(1,round(end/dx**2)//8),
        restart=1+round(end/dx**2),report_interval=max(1,round(end/dx**2)//8),
        label=lambda x,y,z:255 if y in (0,h+1) else 0 if nx//4<=x<3*nx//4 else 100,
        indices=[f'0 "GV" "{force_lattice:.17g}" "1" "0"',f'100 "GV" "{force_lattice:.17g}" "1" "1"',
            '255 "SD" "90" "N/A" "N/A"'])


def measure(r):
    cfg=r['parameters'];h=cfg['height'];dx=cfg['dx'];nx=cfg['size'][0];history=[]
    for filename in r['output_files']:
        t=int(Path(filename).name.split('.')[1]);_,rho,chi,v=geometry_fields(r,t);phi=(1+chi[0,1:-1])/2
        z=np.sum(phi*np.exp(2j*np.pi*(np.arange(nx)+.5)/nx)[None])
        row=dict(step=t,time=(t-r['start'])*dx**2,mean_speed=float(v[0,1:-1,:,0].mean()/dx),
            max_speed=float(np.linalg.norm(v[0,1:-1],axis=-1).max()/dx),
            centroid_modulo=float(np.angle(z)*nx*dx/(2*np.pi)),
            phase_mass=float((phi*rho[0,1:-1]).sum()*dx**2))
        if cfg['physical_amplitude']:
            profile=np.load(cfg['seed_phase_file'])['profile'];delta=phi[:,:nx//2]-profile[None,:nx//2]
            c=np.cos(cfg['transverse_mode']*np.pi*(np.arange(h)+.5)/h)
            row['projected_displacement']=float(-2/h*dx*(delta.sum(axis=1)*c).sum())
        history.append(row)
    pos=np.unwrap([row['centroid_modulo']*2*np.pi/(nx*dx) for row in history])*nx*dx/(2*np.pi)
    for row,x in zip(history,pos):row['centroid_unwrapped']=float(x)
    a,b=history[len(history)//2],history[-1]
    return dict(name=r['name'],parameters=cfg,history=history,
        late_centroid_speed=(b['centroid_unwrapped']-a['centroid_unwrapped'])/(b['time']-a['time']),
        late_mean_speed=b['mean_speed'],solver_sha256=r['solver_sha256'],
        scope='Periodic neutral bridge; circular centroid measures translation. Not the theta60 open-front geometry or a general wetting law.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['fronts','pure','waves','hold'])
    parser.add_argument('--solver',required=True);parser.add_argument('--root',required=True)
    parser.add_argument('--heights',type=int,nargs='+',default=[16,32,48]);parser.add_argument('--jobs',type=int,default=3)
    parser.add_argument('--end',type=float,default=8000.);parser.add_argument('--forces',type=float,nargs='+',default=[-2e-7,-1e-7,1e-7,2e-7])
    parser.add_argument('--length',type=int,default=4);parser.add_argument('--sigma',type=float,default=.005)
    parser.add_argument('--slip',type=float,default=0.,help='Prescribed slip in H16 physical length units; default zero.')
    parser.add_argument('--modes',type=int,nargs='+',default=[1])
    parser.add_argument('--seeds',nargs='+')
    parser.add_argument('--width',type=float,default=WIDTH_REFERENCE,help='Independently prescribed physical W09.')
    args=parser.parse_args()
    if args.width!=WIDTH_REFERENCE and args.stage not in ('fronts','waves'):
        parser.error('--width applies to newly seeded fronts/waves; holds inherit their recorded donor width')
    study=Study(args.solver,args.root);items=[]
    if args.stage=='hold':
        if not args.seeds:parser.error('Holds require recorded periodic-front receipts')
        def resume(filename):
            donor=json.loads(Path(filename).read_text());cfg=dict(donor['parameters']);cut=cfg['steps'];folder=Path(donor['output_files'][-1]).parent
            cfg.update(name=cfg['name']+'-held',steps=2*cut-1,write=max(1,(cut-1)//4),restart=2*cut-1,geometry_file=str(folder/'geometry.raw'))
            return study.run(cfg,start=cut,checkpoint=folder/f'field.{cut:010d}')
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(resume,args.seeds))
        rows=[measure(r) for r in receipts]
        (study.root/'measurements.json').write_text(json.dumps(dict(stage=args.stage,measurements=rows),indent=2)+'\n')
        (study.root/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
        print(json.dumps(dict(stage=args.stage,cases=len(receipts),root=str(study.root))));return
    for h in args.heights:
        if args.stage=='waves':
            choices=[(0.,0.,a,m) for m in args.modes for a in ([-.02,-.01,.01,.02] if h==16 else [-.01,.01])]
        else:choices=[(f,args.sigma,0.,1) for f in args.forces]
        for force,sigma,amplitude,mode in choices:
            cfg=config(h,force,sigma,args.length,args.end,amplitude,mode)
            if args.width!=WIDTH_REFERENCE:
                cfg['physical_width']=args.width;cfg['name']+=f'-W{args.width:g}'
            cfg['physical_slip']=args.slip
            if args.slip:
                cfg['name']+=f'-slip{args.slip:g}'
                cfg['overrides']={'wallSlipLength':format(args.slip/cfg['dx'],'.17g')}
            if args.stage=='pure':
                cfg.update(name=cfg['name']+'-pure',label=lambda x,y,z,hh=h:255 if y in (0,hh+1) else 100)
                items.append((cfg,None));continue
            path=study.root/(cfg['name']+'-seed');beta,phi=seed(path,h,args.length,amplitude,mode,width=args.width)
            cfg.update(beta=beta,seed_phase_file=str(path)+'.phase.npz',seed_phase_sha256=sha(str(path)+'.phase.npz'))
            items.append((cfg,path))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        receipts=list(pool.map(lambda item:study.run(item[0],start=1 if item[1] else 0,checkpoint=item[1]),items))
    rows=[measure(r) for r in receipts]
    (study.root/'measurements.json').write_text(json.dumps(dict(stage=args.stage,measurements=rows),indent=2)+'\n')
    (study.root/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print(json.dumps(dict(stage=args.stage,cases=len(receipts),root=str(study.root))))


if __name__=='__main__':main()
