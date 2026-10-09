#!/usr/bin/env python3
"""Bounded CSF Laplace and zero-surface-tension wall-profile controls."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np
from benchmarks import Study,geometry_fields,measure_capillary,spec
from porous_media import closed_history


def cases():
    for dimension,radius in [(2,12),(2,24),(3,8)]:
        size=3*radius+8
        center=(size-1)/2
        yield spec(f'laplace-csf-d{dimension}-R{radius}',nx=size,ny=size,
                   nz=size if dimension==3 else 1,dimensions=dimension,radius=radius,
                   steps=4000,write=1000,capillary_model='csf',timeout=1800,
                   label=lambda x,y,z,c=center,r=radius,d=dimension:100 if (
                       (x-c)**2+(y-c)**2+((z-c)**2 if d==3 else 0)<=r*r) else 0,
                   indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"'])
    for wall in ['legacy','renormalized','tanh']:
        yield spec(f'neutral-zero-sigma-{wall}',nx=64,ny=18,nz=1,dimensions=2,
                   steps=8000,write=2000,sigma=0,wall_gradient=wall,capillary_model='csf',
                   label=lambda x,y,z:255 if y in (0,17) else (0 if 16<=x<48 else 100),
                   indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"',
                            '255 "SD" "90" "N/A" "N/A"'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True)
    parser.add_argument('--root',required=True)
    parser.add_argument('--initial-root',help='Reuse an existing matching initial campaign')
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    study=Study(args.solver,args.root)
    inputs=list(cases())
    if args.initial_root:
        receipts=[json.loads((Path(args.initial_root)/'results'/f'{s["name"]}-1rank.json').read_text()) for s in inputs]
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(study.run,inputs))
    initial=receipts.copy()

    def extend(pair):
        s,original=pair
        folder=Path(original['output_files'][-1]).parent
        new=dict(s,name=s['name']+'-relaxed',steps=12000,write=2000,restart=12000,
                 geometry_file=str(folder/'geometry.raw'))
        r=study.run(new,1,4000,folder/'field.0000004000')
        r['output_files']=original['output_files']+r['output_files']
        r['seed_solver_sha256']=original['solver_sha256']
        (study.root/'results'/f'{r["name"]}.json').write_text(json.dumps(r,indent=2)+'\n')
        return r

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        relaxed=list(pool.map(extend,[(s,r) for s,r in zip(inputs,receipts)
                                     if s['name'].startswith('laplace') and s['dimensions']==2]))
    replacements={r['parameters']['name'].removesuffix('-relaxed'):r for r in relaxed}
    receipts=[replacements.get(s['name'],r) for s,r in zip(inputs,receipts)]
    measurements=[]
    for r in receipts:
        _,rho,chi,v=geometry_fields(r)
        if r['name'].startswith('laplace'):
            m=measure_capillary(r)
        else:
            m=dict(name=r['name'],max_spurious_speed=float(np.linalg.norm(v,axis=-1).max()),
                   max_wall_color_difference_from_center=float(np.max(abs(chi[0,1]-chi[0,8]))))
        m.update(**closed_history(r),solver_sha256=r['solver_sha256'],
                 input_sha256=r['input_sha256'],geometry_sha256=r['geometry_sha256'])
        measurements.append(m)
    summary=dict(scope='Finite-time bulk Laplace and zero-sigma wall-color diagnostics; D2 drops continued from 4000 to 12000, D3 drop measured at 4000. No three-resolution 3D convergence claim.',
                 measurements=measurements,
                 initial_drop_measurements=[measure_capillary(r) for r in initial if r['name'].startswith('laplace')])
    (study.root/'bulk-and-recolor.json').write_text(json.dumps(summary,indent=2)+'\n')
    for m in measurements:
        print(json.dumps({k:m[k] for k in ('name','relative_laplace_error','late_pressure_relative_change',
                                          'max_spurious_speed','max_wall_color_difference_from_center') if k in m}))


if __name__=='__main__':main()
