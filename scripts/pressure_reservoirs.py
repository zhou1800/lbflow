#!/usr/bin/env python3
"""Guarded pressure-driven slit controls, exact mass budgets and front drift."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path

import numpy as np

from benchmarks import Study, geometry_fields, linear_crossings, sha, spec
from capillary_response import measure as closed_measure


def channel(height=16, spectrum='legacy', phase=1, sign=1, tau=1.55, length_factor=4):
    if length_factor<4:raise ValueError('Pressure collars require a channel length of at least 4H')
    length=length_factor*height;dp=sign*1e-5;steps=500*math.ceil(max(5000,20*length,24*height**2/(((tau-.5)/3)*math.pi**2))/500)
    suffix='' if length_factor==4 else f'-L{length_factor}H'
    return spec(f'pressure-H{height}-{spectrum}{suffix}-phase{phase}-sign{sign}',nx=length,ny=height+2,nz=1,
        dimensions=2,tau=(tau,tau),steps=steps,write=steps//5,restart=steps//2,
        capillary_model='csf',wall_gradient='ghost',phase_boundary='flux',recoloring_model='projection',
        mrt_spectrum=spectrum,height=height,pressure_difference=dp,length_factor=length_factor,
        label=lambda x,y,z:255 if y in (0,height+1) else 50 if x==0 else 51 if x==length-1 else 75,
        indices=[f'75 "VD" "N/A" "N/A" "{phase}"',f'50 "PB" "{1+1.5*dp:.17g}" "1" "{phase}"',
                 f'51 "PB" "{1-1.5*dp:.17g}" "3" "{phase}"','255 "SD" "60" "N/A" "N/A"'])


def budget(receipt):
    folder=Path(receipt['output_files'][-1]).parent
    rows=[json.loads(line[len('LBFlowBoundaryBudget '):]) for line in (folder/'run.log').read_text().splitlines()
          if line.startswith('LBFlowBoundaryBudget ')]
    if not rows:raise AssertionError('Missing reservoir kinetic mass budget')
    nx,ny,nz=receipt['parameters']['size']
    maximum=0.
    for row in rows:
        arrays,rho,chi,v=geometry_fields(receipt,row['step'])
        actual=[float(arrays[k][0,1:-1,1:-1,0].sum()) for k in ('r0','r1')]
        error=float(abs(np.array(actual)-row['current_mass']).max())
        scale=max(1.,sum(row['initial_mass']))
        residual=float(abs(np.array(row['current_mass'])-row['initial_mass']-np.array(row['integrated_net_mass'])).max())
        if error>2e-11*scale or residual>2e-11*scale:raise AssertionError('Independent interior mass/kinetic-flux budget failed')
        maximum=max(maximum,residual/scale)
    return dict(max_relative_budget_residual=maximum,history=rows,
                convention='VD/GV interior excludes PB and solid nodes; restart establishes a new outgoing baseline.')


def channel_measure(receipt):
    cfg=receipt['parameters'];h=cfg['height'];length=cfg['size'][0];tau=cfg['tau'][1]
    _,rho,chi,v=geometry_fields(receipt);x=length//2
    gradient=float((rho[0,1:-1,x-1].mean()-rho[0,1:-1,x+1].mean())/6)
    nominal=cfg['pressure_difference']/(length-1);nu=(tau-.5)/3
    sq=1.92 if cfg['mrt_spectrum']=='liu2015' else 8*(2-1/tau)/(8-1/tau)
    offset=4*(tau-.5)*(1/sq-.5)/3-.25;y=np.arange(h)+.5
    reference=gradient*(y*(h-y)+offset)/(2*nu*rho[0,1:-1,x])
    profile=v[0,1:-1,x,0]
    error=float(abs(profile-reference).max()/abs(reference).max())
    return dict(name=receipt['name'],height=h,tau=tau,mrt_spectrum=cfg['mrt_spectrum'],
        measured_pressure_gradient=gradient,nominal_pressure_gradient=nominal,
        relative_gradient_difference=gradient/nominal-1,constant_offset=offset,
        profile=profile.tolist(),reference=reference.tolist(),relative_profile_error=error,budget=budget(receipt),
        max_transverse_speed=float(abs(v[0,1:-1,x,1]).max()))


def front_seed(root, donor, length_factor=4):
    cfg=donor['parameters'];nx,ny,_=cfg['size'];h=ny-2
    if nx!=4*h:raise ValueError('Recorded donor must be a closed 4H slit')
    folder=Path(donor['output_files'][-1]).parent;source=folder/f'field.{cfg["steps"]:010d}.rst'
    original=np.fromfile(source,dtype='=f8').reshape(2,9,ny,nx)
    if length_factor<4:raise ValueError('Single-front collars require a channel length of at least 4H')
    target_length=length_factor*h;center=target_length//2
    target=np.zeros((2,9,ny,target_length))
    for x in range(target_length):
        donor_x=2*h if x<center-h else x-center+3*h if x<center+h else 0
        target[:,:,:,x]=original[:,:,:,donor_x]
        if x<center-h or x>=center+h:
            phase=0 if x<center-h else 1
            total=target[:,:,:,x].sum(axis=0)
            target[phase,:,:,x]=total;target[1-phase,:,:,x]=0.
    spectrum=cfg.get('mrt_spectrum','legacy');model=cfg.get('recoloring_model','cosine')
    suffix='' if model=='projection' else '-'+model
    if length_factor!=4:suffix+=f'-L{length_factor}H'
    path=root/f'single-front-H{h}-{spectrum}{suffix}'
    Path(str(path)+'.rst').write_bytes(target.astype('=f8').tobytes())
    Path(str(path)+'.rst.meta').write_text(f'DRLAB_OUTGOING_V1\n1 2 2 9\n{ny} {target_length} \n'
        'NEXT_STREAM 2\nINITIAL_COLLISION_INCLUDED 2\nCONVERGENCE_HISTORY_RESET_ON_LOAD\n')
    mapping=dict(donor_name=donor['name'],donor_checkpoint_sha256=sha(source),target_sha256=sha(str(path)+'.rst'),
        mapping='x<C-H: donor 2H and pure phase0; C-H<=x<C+H: donor x-C+3H; x>=C+H: donor 0 and pure phase1. C=L//2; each donor total population is preserved.',
        length_factor=length_factor,target_front_center=center,
        changed_clock='New outgoing step 1; no velocity damping/filter or surface-tension fit.')
    Path(str(path)+'.mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
    return path,mapping


def front_config(donor, ratio, steps, length_factor=4):
    if length_factor<4:raise ValueError('Single-front collars require a channel length of at least 4H')
    source=donor['parameters'];h=source['height'];length=length_factor*h;theta=source['angle'];sigma=source['sigma']
    pc=2*sigma*math.cos(math.radians(theta))/h;dp=ratio*pc
    cfg=dict(source);spectrum=source.get('mrt_spectrum','legacy');model=source.get('recoloring_model','cosine')
    suffix='' if model=='projection' else '-'+model
    if length_factor!=4:suffix+=f'-L{length_factor}H'
    cfg.update(name=f'front-H{h}-{spectrum}{suffix}-ratio{ratio:g}',steps=steps,write=max(1,steps//8),restart=steps,
        mrt_spectrum=spectrum,recoloring_model=model,size=(length,h+2,1),length_factor=length_factor,
        height=h,pressure_ratio=ratio,young_pressure=pc,pressure_difference=dp,
        label=lambda x,y,z:255 if y in (0,h+1) else 50 if x==0 else 51 if x==length-1 else 0 if x<length//2 else 100,
        indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"',
                 f'50 "PB" "{1+1.5*dp:.17g}" "1" "0"',f'51 "PB" "{1-1.5*dp:.17g}" "3" "1"',
                 f'255 "SD" "{theta:.17g}" "N/A" "N/A"'])
    cfg.pop('geometry_file',None)
    return cfg


def front_measure(receipt):
    cfg=receipt['parameters'];h=cfg['height'];center=cfg['size'][0]//2;history=[]
    for filename in receipt['output_files']:
        t=int(Path(filename).name.split('.')[1]);_,rho,chi,v=geometry_fields(receipt,t)
        phi=(chi[0,1:-1,1:-1]+1)/2
        x=.5+float(phi.sum())/h
        center_crossings=linear_crossings(chi[0,h//2+1],0)
        wall_crossings=linear_crossings(chi[0,1],0)
        center_cross=min(center_crossings,key=lambda a:abs(a-center))
        wall=min(wall_crossings,key=lambda a:abs(a-center))
        history.append(dict(step=t,volume_front=x,centerline_front=center_cross,first_fluid_row_front=wall,
            max_speed=float(np.linalg.norm(v,axis=-1).max()),mean_speed=float(v[0,1:-1,1:-1,0].mean())))
    left=history[len(history)//2];right=history[-1]
    speed=(right['volume_front']-left['volume_front'])/(right['step']-left['step'])
    return dict(name=receipt['name'],height=h,length_factor=cfg.get('length_factor',4),mrt_spectrum=cfg['mrt_spectrum'],pressure_ratio=cfg['pressure_ratio'],
        recoloring_model=cfg['recoloring_model'],
        pressure_difference=cfg['pressure_difference'],young_pressure=cfg['young_pressure'],history=history,
        late_volume_front_speed=speed,budget=budget(receipt),
        scope='Finite-time drift in a flat single-front slit. Volume-equivalent position is not the contact-line position; no throat threshold or hysteresis claim.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['channels','fronts','hold-fronts'])
    parser.add_argument('--solver',required=True);parser.add_argument('--root',required=True)
    parser.add_argument('--seeds',nargs='+');parser.add_argument('--ratios',nargs='+',type=float,default=[.9,.95,1.,1.05,1.1])
    parser.add_argument('--steps',type=int,default=8000);parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--tau',type=float,default=1.55)
    parser.add_argument('--length-factor',type=int,default=4)
    parser.add_argument('--closed-centered-halfwidth',type=float,
        help='Prespecified pressure ratios centered on the donor closed pressure, with this halfwidth in Young-pressure units.')
    parser.add_argument('--force-removed',action='store_true',help='Changed-physics sigma=0 diagnostic at unchanged reservoir pressures.')
    args=parser.parse_args();study=Study(args.solver,args.root)
    if args.length_factor<4:parser.error('Pressure channel length must be at least 4H')
    if args.stage=='channels':
        configs=[channel(h,spectrum,phase,sign,args.tau,args.length_factor) for h in (8,16,32)
                 for spectrum in ('legacy','liu2015') for phase,sign in ((0,-1),(1,1))]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(study.run,configs))
        measurements=[channel_measure(r) for r in receipts]
    elif args.stage=='fronts':
        if not args.seeds:parser.error('Fronts require held closed-slit donor receipts')
        configs=[]
        for filename in args.seeds:
            donor=json.loads(Path(filename).read_text());seed,mapping=front_seed(study.root,donor,args.length_factor)
            ratios=args.ratios
            if args.closed_centered_halfwidth is not None:
                if args.closed_centered_halfwidth<=0:parser.error('Pressure halfwidth must be positive')
                observed=closed_measure(donor)
                center=1+observed['relative_pressure_error']
                ratios=[center-args.closed_centered_halfwidth,center+args.closed_centered_halfwidth]
                mapping['pressure_protocol']=dict(closed_pressure_ratio=center,halfwidth=args.closed_centered_halfwidth,
                    rule='Both signs tested around the independently measured closed donor pressure; input sigma and beta unchanged.')
            for ratio in ratios:
                cfg=front_config(donor,ratio,args.steps*(donor['parameters']['height']//16)**2,args.length_factor)
                if args.force_removed:
                    cfg['sigma']=0.;mapping['force_removal']='Diagnostic continuation: sigma=0, unchanged donor populations and reservoir pressure difference.'
                cfg['donor_mapping']=mapping;configs.append((cfg,seed))
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(lambda item:study.run(item[0],start=1,checkpoint=item[1]),configs))
        measurements=[front_measure(r) for r in receipts]
    else:
        if not args.seeds:parser.error('Front holds require open-front receipts')
        configs=[]
        for filename in args.seeds:
            donor=json.loads(Path(filename).read_text());cfg=dict(donor['parameters']);cut=cfg['steps']
            folder=Path(donor['output_files'][-1]).parent
            cfg.update(name=cfg['name']+'-held',steps=4*cut,write=cut//2,restart=4*cut,geometry_file=str(folder/'geometry.raw'))
            configs.append((cfg,cut,folder/f'field.{cut:010d}',sha(filename)))
        def run(item):
            cfg,cut,path,receipt_hash=item;r=study.run(cfg,start=cut,checkpoint=path);r['seed_receipt_sha256']=receipt_hash
            (study.root/'results'/f'{r["name"]}.json').write_text(json.dumps(r,indent=2)+'\n')
            return r
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(run,configs))
        measurements=[front_measure(r) for r in receipts]
    (study.root/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
    (study.root/'measurements.json').write_text(json.dumps(dict(stage=args.stage,measurements=measurements,solver_sha256=study.solver_sha),indent=2)+'\n')
    print(json.dumps(dict(stage=args.stage,cases=len(measurements),root=str(study.root),
        maximum_budget_residual=max(r['budget']['max_relative_budget_residual'] for r in measurements)),indent=2))


if __name__=='__main__':main()
