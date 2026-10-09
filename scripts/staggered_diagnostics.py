#!/usr/bin/env python3
"""Measure adjacent physical fields without filtering or changing the solver."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np

from benchmarks import Study,geometry_fields,sha


def measure(original,resumed):
    cut=original['parameters']['steps']
    _,rho0,chi0,v0=geometry_fields(original)
    _,rho1,chi1,v1=geometry_fields(resumed,cut+1)
    _,rho2,chi2,v2=geometry_fields(resumed,cut+2)
    fluid=rho0>0
    z,y,x=np.indices(rho0.shape)
    fractions=[];amplitudes=[]
    for v,rho in [(v0,rho0),(v1,rho1),(v2,rho2)]:
        momentum=v*rho[...,None]
        a=[];f=[]
        for d,coord in enumerate((x,y)):
            component=momentum[...,d][fluid];pattern=(-1.)**coord[fluid]
            amplitude=float(np.mean(component*pattern));energy=float(np.sum(component**2))
            a.append(amplitude);f.append(float(len(component)*amplitude**2/energy) if energy else 0.)
        amplitudes.append(a);fractions.append(f)
    pair=(v0+v1)/2
    return dict(name=original['name'],parameters=original['parameters'],
                raw_speed=float(np.linalg.norm(v0[fluid],axis=-1).max()),
                adjacent_pair_average_speed=float(np.linalg.norm(pair[fluid],axis=-1).max()),
                adjacent_velocity_rms=float(np.sqrt(np.mean((v1[fluid]-v0[fluid])**2))),
                two_step_velocity_rms=float(np.sqrt(np.mean((v2[fluid]-v0[fluid])**2))),
                adjacent_max_color_change=float(abs(chi1[fluid]-chi0[fluid]).max()),
                two_step_max_color_change=float(abs(chi2[fluid]-chi0[fluid]).max()),
                staggered_momentum_amplitudes=amplitudes,staggered_component_energy_fractions=fractions,
                tangent_velocity_rows=v0[0,1:5,:,0].tolist(),
                solver_sha256=resumed['solver_sha256'],seed_solver_sha256=original['solver_sha256'],
                checkpoint_sha256=resumed['checkpoint_sha256'],
                scope='Adjacent-step averages and staggered projections are diagnostics; the raw solver fields are retained and not filtered.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',required=True);p.add_argument('--root',required=True)
    p.add_argument('--receipts',nargs='+',required=True);p.add_argument('--jobs',type=int,default=2)
    args=p.parse_args();study=Study(args.solver,args.root)
    def run(path):
        original=json.loads(Path(path).read_text());cfg=original['parameters'];cut=cfg['steps']
        folder=Path(original['output_files'][-1]).parent
        new=dict(cfg,name=cfg['name']+'-adjacent',steps=cut+2,write=1,restart=cut+2,
                 geometry_file=str(folder/'geometry.raw'),report_dir=0)
        resumed=study.run(new,start=cut,checkpoint=folder/f'field.{cut:010d}')
        resumed['seed_receipt_sha256']=sha(path)
        (study.root/'results'/f'{resumed["name"]}.json').write_text(json.dumps(resumed,indent=2)+'\n')
        m=measure(original,resumed);m['seed_receipt']=str(Path(path).resolve())
        return m
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:measurements=list(pool.map(run,args.receipts))
    summary=dict(protocol='Replay two adjacent steps from identified outgoing checkpoints; compare the physical half-force velocities.',measurements=measurements)
    (study.root/'staggered.json').write_text(json.dumps(summary,indent=2)+'\n')
    for m in measurements:print(json.dumps({k:m[k] for k in ('name','raw_speed','adjacent_pair_average_speed','staggered_component_energy_fractions')}))


if __name__=='__main__':main()
