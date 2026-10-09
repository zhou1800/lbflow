#!/usr/bin/env python3
"""Direct pressure and equilibrium diagnostics for closed two-phase slit menisci."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path

import numpy as np

from benchmarks import Study, geometry_fields, interface_points, fit_sphere, spec
from porous_media import closed_history


def cases(heights, angles, walls, beta=.7, capillary_models=('perturbation',)):
    for height in heights:
        scale = height/16
        length = 4*height
        steps = int(8000*scale**2)
        for angle in angles:
            for wall in walls:
                for capillary in capillary_models:
                    suffix='' if capillary=='perturbation' else f'-{capillary}'
                    yield spec(f'menisci-H{height}-theta{angle:g}-{wall}-beta{beta:g}{suffix}',
                               nx=length, ny=height+2, nz=1, dimensions=2,
                               tau=(1., 1.), sigma=.005/scale, beta=beta,
                               steps=steps, write=steps//4, restart=steps,
                               height=height, angle=angle, wall_gradient=wall,capillary_model=capillary,
                               label=lambda x,y,z,h=height: 255 if y in (0,h+1) else (
                                   0 if h <= x < 3*h else 100),
                               indices=['0 "VD" "N/A" "N/A" "0"',
                                        '100 "VD" "N/A" "N/A" "1"',
                                        f'255 "SD" "{angle:.17g}" "N/A" "N/A"'], timeout=3600)


def bulk_pressure(rho, chi, height, level=.99, transverse_margin=0):
    z,y,x = np.indices(rho.shape)
    fluid = rho > 0
    transverse = (y >= 1+transverse_margin) & (y <= height-transverse_margin)
    phase0 = fluid & transverse & (x >= 1.5*height) & (x < 2.5*height) & (chi > level)
    phase1 = fluid & transverse & ((x < .5*height) | (x >= 3.5*height)) & (chi < -level)
    if phase0.sum() < 8 or phase1.sum() < 8:
        raise ValueError('Insufficient separated bulk samples')
    p0,p1 = [float(rho[mask].mean()/3) for mask in (phase0,phase1)]
    return dict(pressure_jump=p0-p1, pressure0=p0, pressure1=p1,
                sample_counts=[int(phase0.sum()),int(phase1.sum())],
                color_level=level, transverse_margin=transverse_margin)


def measure(result):
    p = result['parameters']
    height, sigma, angle = p['height'], p['sigma'], p['angle']
    audit = closed_history(result)
    history=[]
    for frame_index,snapshot in enumerate(audit['history'][1:],1):
        at=snapshot['step']
        _,rho,chi,velocity=geometry_fields(result,at)
        prior=geometry_fields(result,audit['history'][frame_index-1]['step'])
        fluid=rho>0
        pressure=bulk_pressure(rho,chi,height,transverse_margin=height//4)
        entry=dict(audit['history'][len(history)+1])
        entry.update(pressure,
                     pressure_coefficient=pressure['pressure_jump']*height/(2*sigma),
                     color_rms_change=float(np.sqrt(np.mean((chi[fluid]-prior[2][fluid])**2))),
                     rms_speed=float(np.sqrt(np.mean(np.sum(velocity[fluid]**2,axis=-1)))))
        history.append(entry)
    _,rho,chi,velocity=geometry_fields(result)
    sensitivity=[]
    for level in (.98,.99,.999):
        for margin in (0,height//4):
            try:
                sensitivity.append(bulk_pressure(rho,chi,height,level,margin))
            except ValueError as error:
                sensitivity.append(dict(color_level=level,transverse_margin=margin,error=str(error)))
    # Fit only the left interface and retain its extrapolation sensitivity.
    circles=[]
    if abs(angle-90)>1e-8:
        width={.7:4.765235021842287,.35:9.376436760977107,.9:3.929187708575192}.get(p['beta'])
        if width is not None:
            points=interface_points(chi,0.,rho>0)
            for exclusion in (.5,.75,1.):
                selected=points[(points[:,0]<2*height) & (points[:,1]>.5+exclusion*width)
                                & (points[:,1]<height+.5-exclusion*width)]
                try:
                    if len(selected)<8: raise ValueError('Insufficient interface points')
                    center,radius,rms=fit_sphere(selected)
                    sign=1 if center[0]>selected[:,0].mean() else -1
                    cosines=[sign*(center[1]-.5)/radius,sign*(height+.5-center[1])/radius]
                    if max(map(abs,cosines))>1: raise ValueError('Fit misses physical walls')
                    circles.append(dict(exclusion_widths=exclusion,radius=float(radius),
                                        fit_rms=rms,center=center.tolist(),points=len(selected),
                                        angles=[math.degrees(math.acos(x)) for x in cosines],
                                        laplace_from_fit=sign*sigma/radius))
                except (ValueError,np.linalg.LinAlgError) as error:
                    circles.append(dict(exclusion_widths=exclusion,error=str(error)))
    reference=2*sigma*math.cos(math.radians(angle))/height
    final=history[-1]
    return dict(name=result['name'],height=height,angle=angle,wall_gradient=p['wall_gradient'],
                capillary_model=p.get('capillary_model','perturbation'),
                wall_recoloring=p.get('wall_recoloring','legacy'),
                curvature_stencil=p.get('curvature_stencil','normal'),
                beta=p['beta'],sigma=sigma,steps=p['steps'],ranks=result['ranks'],
                reference_pressure_jump=reference,reference_pressure_coefficient=math.cos(math.radians(angle)),
                signed_pressure_coefficient_error=final['pressure_coefficient']-math.cos(math.radians(angle)),
                relative_pressure_error=(final['pressure_jump']/reference-1) if abs(reference)>1e-12 else None,
                late_pressure_coefficient_change=final['pressure_coefficient']-history[-2]['pressure_coefficient'],
                history=history,pressure_sampling_sensitivity=sensitivity,circle_fits=circles,
                mass=audit['max_relative_mass_drift'],solver_sha256=result['solver_sha256'],
                input_sha256=result['input_sha256'],geometry_sha256=result['geometry_sha256'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True)
    parser.add_argument('--root',required=True)
    parser.add_argument('--heights',nargs='+',type=int,default=[16,32,48])
    parser.add_argument('--angles',nargs='+',type=float,default=[60,90,120])
    parser.add_argument('--walls',nargs='+',choices=['legacy','renormalized','tanh'],default=['legacy','renormalized'])
    parser.add_argument('--beta',type=float,default=.7)
    parser.add_argument('--capillary-models',nargs='+',choices=['perturbation','csf'],default=['perturbation'])
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--ranks',type=int,default=1)
    args=parser.parse_args()
    if any(h<16 or h%16 for h in args.heights): parser.error('Heights must be positive multiples of 16')
    study=Study(args.solver,args.root)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        receipts=list(pool.map(lambda p:study.run(p,args.ranks),cases(args.heights,args.angles,args.walls,args.beta,args.capillary_models)))
    measurements=[measure(r) for r in receipts]
    summary=dict(protocol='Closed periodic slit with two menisci, no body force, halfway walls and contact angle through phase 1.',
                 scope='Direct finite-time pressure and field measurements; static Young-Laplace controls. No entry-pressure, hysteresis or curved-wall accuracy claim.',
                 measurements=measurements)
    (study.root/'capillary-response.json').write_text(json.dumps(summary,indent=2)+'\n')
    for m in measurements:
        print(json.dumps(dict(name=m['name'],relative_pressure_error=m['relative_pressure_error'],
                              pressure_coefficient=m['history'][-1]['pressure_coefficient'],
                              speed=m['history'][-1]['max_speed'],late_color=m['history'][-1]['color_rms_change'])))


if __name__=='__main__':main()
