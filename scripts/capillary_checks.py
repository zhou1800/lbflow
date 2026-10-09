#!/usr/bin/env python3
"""Physical force, normal-exchange and replay contracts for experimental CSF."""
import argparse
import json
from pathlib import Path
import uuid

import numpy as np

from benchmarks import Study, compare, geometry_fields, spec
from porous_media import closed_history


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True)
    parser.add_argument('--root',required=True)
    parser.add_argument('--mpiexec',default='mpiexec')
    parser.add_argument('--numproc-flag',default='-n')
    args=parser.parse_args()
    study=Study(args.solver,Path(args.root)/('capillary-checks-'+uuid.uuid4().hex[:10]),
                args.mpiexec,args.numproc_flag)
    controls=[]
    for dimension in (2,3):
        for collision in ('SRT','MRT'):
            common=dict(nx=32,ny=4,nz=1 if dimension==2 else 4,dimensions=dimension,
                        eq=collision,capillary_model='csf',steps=600,write=300)
            plane=spec(f'csf-plane-d{dimension}-{collision}',**common,
                       label=lambda x,y,z:0 if x<16 else 100,
                       indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"'])
            result=study.run(plane)
            _,rho,chi,velocity=geometry_fields(result)
            speed=float(np.linalg.norm(velocity,axis=-1).max())
            initial_speed=float(np.linalg.norm(geometry_fields(result,0)[3],axis=-1).max())
            mass=closed_history(result)['max_relative_mass_drift']
            if speed>2e-14 or initial_speed>2e-14:
                raise AssertionError('Bulk planar interface failed mechanical equilibrium')
            controls.append(dict(name=result['name'],max_speed=speed,initial_speed=initial_speed,mass=mass))
            forced=spec(f'csf-force-d{dimension}-{collision}',**common,
                        label=lambda x,y,z:75,indices=['75 "GV" "0.000001" "1" "0"'])
            result=study.run(forced)
            velocity=geometry_fields(result)[3]
            error=float(np.max(abs(velocity[...,0]-6e-4)))
            transverse=float(np.max(abs(velocity[...,1:])))
            if max(error,transverse)>1e-12:raise AssertionError('CSF changed uniform force integration')
            controls.append(dict(name=result['name'],velocity_error=error,transverse_speed=transverse))
    drop=spec('csf-drop-continuous',nx=12,ny=12,nz=12,steps=40,write=20,restart=20,
              capillary_model='csf',wall_gradient='tanh',
              label=lambda x,y,z:100 if (x-5.5)**2+(y-5.5)**2+(z-5.5)**2<=9 else 0,
              indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"'])
    continuous=study.run(drop)
    parallel=study.run(dict(drop,name='csf-drop-parallel'),ranks=2)
    comparisons=[compare(continuous,parallel)]
    checkpoint=Path(continuous['output_files'][-1]).parent/'field.0000000020'
    for ranks in (1,2):
        resumed=study.run(dict(drop,name=f'csf-drop-resumed-{ranks}'),ranks,20,checkpoint)
        comparisons.append(compare(continuous,resumed))
    # A contact line on a rank boundary exercises the reconstructed normal halo.
    wall=spec('csf-wall-continuous',nx=16,ny=10,nz=1,dimensions=2,steps=100,write=50,
              capillary_model='csf',wall_gradient='tanh',
              label=lambda x,y,z:255 if y in (0,9) else (0 if 4<=x<12 else 100),
              indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"',
                       '255 "SD" "60" "N/A" "N/A"'])
    continuous=study.run(wall)
    parallel=study.run(dict(wall,name='csf-wall-parallel'),ranks=2)
    comparisons.append(compare(continuous,parallel))
    if not all(c['bitwise_fields'] for c in comparisons):
        raise AssertionError('CSF MPI/restart fields differ')
    failures=[]
    base=spec('guard',nx=4,ny=4,nz=1,dimensions=2,steps=2,write=1,capillary_model='csf')
    for name,overrides in [('invalid-capillary',{'capillaryModel':'"missing"'}),
                           ('invalid-wall',{'wallGradient':'"missing"'}),
                           ('unsupported-density',{'density':'1 2'}),
                           ('unsupported-pressure-weights',{'alpha0':.5})]:
        result=study.run(dict(base,name=name,overrides=overrides,expect_failure=True))
        failures.append(dict(name=name,exit=result['exit']))
    for name,indices in [('pressure-boundary',['0 "PB" "1" "1" "0"']),
                          ('velocity-boundary',['0 "VB" "0" "1" "0"']),
                          ('selective-film',['0 "HF" "N/A" "N/A" "0"'])]:
        result=study.run(dict(base,name=name,indices=indices,expect_failure=True))
        failures.append(dict(name=name,exit=result['exit']))
    summary=dict(controls=controls,comparisons=comparisons,failures=failures,
                 solver_sha256=study.solver_sha,
                 scope='Bulk plane/forcing and short-time replay tests; these do not validate curved 3D wetting, contact-line accuracy or high density ratios.')
    (study.root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
