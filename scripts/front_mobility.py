#!/usr/bin/env python3
"""Fixed-width, length and rate controls for capillary-front resistance."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path

import numpy as np

from benchmarks import Study, geometry_fields, linear_crossings, sha, spec
from capillary_response import cases, measure as closed_measure
from phase_recurrence import bulk_profile, coefficients
from pressure_reservoirs import channel, channel_measure, front_config, front_measure, front_seed

WIDTH_REFERENCE=4.765235021846117


def zero_flux_profile(beta):
    """Symmetric descending projection profile, independent of time marching."""
    if not 0<beta<1:raise ValueError('The monotone quadratic branch requires 0<beta<1')
    high=(1+(math.sqrt(1+beta*beta)-1)/beta)/2
    phi=[1-high]
    while phi[-1]>.01:
        current=phi[-1];constant=(1-beta)*current+beta*current*current
        discriminant=(1+beta)**2-4*beta*constant
        phi.append(2*constant/(1+beta+math.sqrt(discriminant)))
        if len(phi)>10000:raise ValueError('Requested interface exceeds calibration support')
    x=np.arange(len(phi))+.5
    width=2*float(np.interp(.05,np.array(phi)[::-1],x[::-1]))
    return width,np.array(phi)


def match_beta(width):
    low,high=.01,.9
    if not zero_flux_profile(high)[0]<width<zero_flux_profile(low)[0]:
        raise ValueError('Target width outside monotone calibration range')
    for _ in range(56):
        middle=(low+high)/2
        if zero_flux_profile(middle)[0]>width:low=middle
        else:high=middle
    return (low+high)/2


def save(study,filename,data):
    (study.root/filename).write_text(json.dumps(data,indent=2)+'\n')


def calibration(study,args):
    configs=[]
    for height in args.heights:
        length=max(128,8*height);target=WIDTH_REFERENCE*height/16;beta=match_beta(target)
        configs.append(spec(f'mobility-planar-H{height}',nx=length,ny=4,nz=1,dimensions=2,
            recoloring_model='projection',capillary_model='csf',mrt_spectrum='legacy',beta=beta,
            steps=8000,write=4000,height=height,width_target=target,
            label=lambda x,y,z,n=length:0 if x<n//2 else 100,
            indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"']))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(study.run,configs))
    rows=[]
    ratio=coefficients()['bulk_segregation']/coefficients()['bulk_diffusion']
    for r in receipts:
        cfg=r['parameters'];_,rho,chi,v=geometry_fields(r);profile=chi[0,0]
        middle=linear_crossings(profile,0)[0]
        a=min(linear_crossings(profile,.9),key=lambda x:abs(x-middle))
        b=min(linear_crossings(profile,-.9),key=lambda x:abs(x-middle));width=abs(a-b)
        marched=bulk_profile(cfg['size'][0],cfg['beta']/ratio)
        marched=np.roll(marched,-cfg['size'][0]//4)
        scalar_error=float(abs((profile+1)/2-marched).max())
        if abs(width/cfg['width_target']-1)>2e-8 or scalar_error>2e-11:
            raise AssertionError('Full-solver width does not verify independent zero-flux matching')
        rows.append(dict(name=r['name'],height=cfg['height'],beta=cfg['beta'],width_09=width,
            normalized_width=width/cfg['height'],independent_width=zero_flux_profile(cfg['beta'])[0],
            max_scalar_profile_error=scalar_error,max_speed=float(np.linalg.norm(v,axis=-1).max()),
            late_color_change=float(abs(chi-geometry_fields(r,4000)[2]).max()),profile=profile.tolist()))
    save(study,'calibration.json',dict(measurements=rows,solver_sha256=study.solver_sha,
        protocol='Independent zero-flux recurrence chooses beta for constant measured W09/H; full solver and time-marched scalar verify the choice. No pressure or speed fitting.'))
    return receipts


def closed(study,args):
    if not args.calibration:raise ValueError('Closed study requires independent width calibration')
    cal=json.loads(Path(args.calibration).read_text())['measurements'];configs=[]
    for row in cal:
        h=row['height']
        if h not in args.heights:continue
        cfg=next(cases([h],[60],['ghost'],row['beta'],['csf']))
        cfg.update(name=f'mobility-closed-H{h}',steps=int(32000*(h/16)**2),
            write=int(8000*(h/16)**2),restart=int(32000*(h/16)**2),
            phase_boundary='flux',recoloring_model='projection',mrt_spectrum='legacy',
            planar_width_09=row['width_09'],calibration_sha256=sha(args.calibration))
        configs.append(cfg)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(study.run,configs))
    rows=[closed_measure(r) for r in receipts];save(study,'closed.json',dict(measurements=rows,solver_sha256=study.solver_sha))
    return receipts


def channels(study,args):
    configs=[]
    for h in args.heights:
        lengths=(4,8) if h==16 else (4,)
        for length in lengths:
            for phase,sign in ((0,-1),(1,1)):
                cfg=channel(h,'legacy',phase,sign,tau=1.,length_factor=length)
                cfg['name']='mobility-'+cfg['name'];configs.append(cfg)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(study.run,configs))
    rows=[]
    for r in receipts:
        row=channel_measure(r);_,rho,chi,v=geometry_fields(r);cfg=r['parameters']
        row.update(length_factor=cfg['length_factor'],pressure_difference=cfg['pressure_difference'],
            volume_conductance=float(v[0,1:-1,1:-1,0].mean()/cfg['pressure_difference']),
            center_conductance=float(np.mean(row['profile'])/cfg['pressure_difference']))
        if row['relative_profile_error']>1e-4:raise AssertionError('Pure hydraulic control did not reach its independent reference')
        rows.append(row)
    save(study,'channels.json',dict(measurements=rows,solver_sha256=study.solver_sha));return receipts


def donors(root):
    for path in sorted((Path(root)/'results').glob('*.json')):
        data=json.loads(path.read_text())
        if isinstance(data,dict) and 'parameters' in data:yield path,data


def front_runs(study,args):
    if not args.closed_root:raise ValueError('Front study requires independent closed donors')
    configs=[]
    for path,donor in donors(args.closed_root):
        h=donor['parameters']['height']
        if h not in args.heights:continue
        observation=closed_measure(donor);center=1+observation['relative_pressure_error']
        folder=study.root/'seeds'/donor['parameters']['name'];folder.mkdir(parents=True)
        lengths=(4,8) if h==16 else (4,)
        for length in lengths:
            seed,mapping=front_seed(folder,donor,length)
            offsets=args.offsets if h==16 else [-.02,0.,.02]
            for offset in offsets:
                cfg=front_config(donor,center+offset,int(args.steps*(h/16)**2),length)
                cfg.update(name=f'mobility-front-H{h}-L{length}H-offset{offset:g}',
                    pressure_offset=offset,closed_pressure_ratio=center,donor_receipt_sha256=sha(path),
                    donor_mapping=mapping,write=max(1,int(args.steps*(h/16)**2)//8))
                configs.append((cfg,seed))
    def run(item):
        cfg,seed=item
        try:return study.run(cfg,start=1,checkpoint=seed)
        except RuntimeError as error:
            # Failed regimes are evidence, not measurements or automatic guards.
            r=json.loads((study.root/'results'/f'{cfg["name"]}-1rank.json').read_text())
            r['measurement_failure']=str(error);return r
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(run,configs))
    rows=[front_measure(r) for r in receipts if r['exit']==0 and r['finished']]
    save(study,'fronts.json',dict(measurements=rows,failures=[r['name'] for r in receipts if r['exit']!=0],solver_sha256=study.solver_sha))
    return receipts


def holds(study,args):
    if not args.front_root:raise ValueError('Holds require recorded front receipts')
    seeds=[(path,r) for path,r in donors(args.front_root) if r['exit']==0 and r['finished'] and r['parameters']['height'] in args.heights]
    def run(item):
        path,r=item;cfg=dict(r['parameters']);cut=cfg['steps'];folder=Path(r['output_files'][-1]).parent
        cfg.update(name=cfg['name']+'-held',steps=2*cut,write=max(1,cut//4),restart=2*cut,geometry_file=str(folder/'geometry.raw'))
        try:
            fresh=study.run(cfg,start=cut,checkpoint=folder/f'field.{cut:010d}')
            fresh['seed_receipt_sha256']=sha(path)
            save(study,f'results/{fresh["name"]}.json',fresh);return fresh
        except RuntimeError as error:
            fresh=json.loads((study.root/'results'/f'{cfg["name"]}-1rank.json').read_text());fresh['measurement_failure']=str(error);return fresh
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(run,seeds))
    save(study,'fronts.json',dict(measurements=[front_measure(r) for r in receipts if r['exit']==0 and r['finished']],
        failures=[r['name'] for r in receipts if r['exit']!=0],solver_sha256=study.solver_sha));return receipts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['calibrate','closed','channels','fronts','hold'])
    parser.add_argument('--solver',required=True);parser.add_argument('--root',required=True)
    parser.add_argument('--jobs',type=int,default=3);parser.add_argument('--heights',type=int,nargs='+',default=[16,32,48])
    parser.add_argument('--calibration');parser.add_argument('--closed-root');parser.add_argument('--front-root')
    parser.add_argument('--steps',type=int,default=8000)
    parser.add_argument('--offsets',type=float,nargs='+',default=[-.04,-.02,-.01,0.,.01,.02,.04])
    args=parser.parse_args();study=Study(args.solver,args.root)
    stages={'calibrate':calibration,'closed':closed,'channels':channels,'fronts':front_runs,'hold':holds}
    receipts=stages[args.stage](study,args);save(study,'receipts.json',receipts)
    print(json.dumps(dict(stage=args.stage,cases=len(receipts),failed=sum(r['exit']!=0 for r in receipts),root=str(study.root))))


if __name__=='__main__':main()
