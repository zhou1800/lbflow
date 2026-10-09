#!/usr/bin/env python3
"""Measure wall-phase consistency and mechanical equilibrium on flat slits."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import product
import json
from pathlib import Path

import numpy as np
from benchmarks import Study,geometry_fields,spec
from porous_media import closed_history
from capillary_response import bulk_pressure


def cases(heights,betas,methods,sigmas,models,collisions=('MRT',),axes=('x',),curvatures=('normal',)):
    for h,beta,method,sigma0,model,collision,axis,curvature in product(
            heights,betas,methods,sigmas,models,collisions,axes,curvatures):
        if model!='csf' and curvature!='normal':continue
        length=4*h
        size=(length,h+2) if axis=='x' else (h+2,length)
        steps=int(8000*(h/16)**2)
        def label(x,y,z,h=h,axis=axis):
            tangent,normal=(x,y) if axis=='x' else (y,x)
            return 255 if normal in (0,h+1) else (0 if h<=tangent<3*h else 100)
        name=f'wall-balance-H{h}-beta{beta:g}-{method}-sigma{sigma0:g}-{model}-{collision}-{axis}'
        if curvature!='normal':name+='-'+curvature
        yield spec(name,nx=size[0],ny=size[1],nz=1,dimensions=2,
                   height=h,beta=beta,sigma=sigma0/(h/16),steps=steps,write=steps//4,
                   wall_gradient='renormalized',wall_recoloring=method,
                   capillary_model=model,curvature_stencil=curvature,eq=collision,axis=axis,label=label,timeout=1800,
                   indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"',
                            '255 "SD" "90" "N/A" "N/A"'])


def measure(receipt):
    p=receipt['parameters'];h=p['height']
    _,rho,chi,v=geometry_fields(receipt)
    previous=geometry_fields(receipt,p['steps']-p['write'])
    fluid=rho>0
    drift=float(np.sqrt(np.mean((chi[fluid]-previous[2][fluid])**2)))
    if p['axis']=='y':rho,chi,v=rho.transpose(0,2,1),chi.transpose(0,2,1),v.transpose(0,2,1,3)[...,[1,0]]
    phi=(chi[0,1:-1]+1)/2
    center=(phi[h//2-1]+phi[h//2])/2
    wall=(phi[0]+phi[-1])/2
    speed=float(np.linalg.norm(v,axis=-1).max())
    mass=closed_history(receipt)
    return dict(name=receipt['name'],height=h,beta=p['beta'],axis=p['axis'],collision=p['eq'],
                wall_recoloring=p['wall_recoloring'],capillary_model=p['capillary_model'],sigma=p['sigma'],steps=p['steps'],
                curvature_stencil=p.get('curvature_stencil','normal'),
                max_speed=speed,residual_capillary_number=(speed/6/p['sigma']) if p['sigma'] else None,
                max_wall_color_difference=float(np.max(abs(2*(wall-center)))),
                max_transverse_color_difference=float(np.max(abs(2*(phi-center)))),
                wall_phase_profile=wall.tolist(),center_phase_profile=center.tolist(),
                late_color_rms_change=drift,bulk_pressure=bulk_pressure(rho,chi,h,transverse_margin=h//4),
                mass=mass['max_relative_mass_drift'],history=mass['history'],
                solver_sha256=receipt['solver_sha256'],input_sha256=receipt['input_sha256'],geometry_sha256=receipt['geometry_sha256'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True);parser.add_argument('--root',required=True)
    parser.add_argument('--heights',type=int,nargs='+',default=[16])
    parser.add_argument('--betas',type=float,nargs='+',default=[.35,.7,.9])
    parser.add_argument('--methods',nargs='+',choices=['legacy','tangential'],default=['legacy','tangential'])
    parser.add_argument('--sigmas',type=float,nargs='+',default=[0,.005])
    parser.add_argument('--models',nargs='+',choices=['csf','perturbation'],default=['csf'])
    parser.add_argument('--collisions',nargs='+',choices=['MRT','SRT'],default=['MRT'])
    parser.add_argument('--axes',nargs='+',choices=['x','y'],default=['x'])
    parser.add_argument('--curvature-stencils',nargs='+',choices=['normal','gradient'],default=['normal'])
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    study=Study(args.solver,args.root)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        receipts=list(pool.map(study.run,cases(args.heights,args.betas,args.methods,args.sigmas,args.models,args.collisions,args.axes,args.curvature_stencils)))
    measurements=[measure(r) for r in receipts]
    result=dict(protocol='Neutral flat Cartesian walls; closed periodic two-phase slit, equal unit densities and viscosities, tau=1; step initialization; diffusive sigma/time scaling.',
                scope='Discrete phase-profile and mechanical-equilibrium controls. Joint grid/interface sequence; no nonneutral, curved-wall, dynamic contact-line or trapping claim.',measurements=measurements)
    (study.root/'wall-consistency.json').write_text(json.dumps(result,indent=2)+'\n')
    for m in measurements:print(json.dumps({k:m[k] for k in ['name','max_speed','max_wall_color_difference','late_color_rms_change']}))


if __name__=='__main__':main()
