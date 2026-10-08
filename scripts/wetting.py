#!/usr/bin/env python3
"""Flat-wall wetting size/interface-width sensitivity; no new-method claim."""
import argparse,json,math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from benchmarks import Study,spec,fields,geometry_fields,interface_points,fit_sphere


WIDTHS={.7:4.765235021842287,.35:9.376436760977107}


def cases():
    runs=[]
    for radius,beta in [(12,.7),(24,.7),(36,.7),(24,.35),(48,.35)]:
        nx,ny=6*radius+8,3*radius+8
        center=(nx-1)/2
        steps=24000*(radius//12)
        for angle in (60,90,120):
            for method in ('legacy','renormalized'):
                label=lambda x,y,z,n=ny,c=center,r=radius:255 if y in (0,n-1) else (
                    100 if (x-c)**2+(y-.5)**2<=r*r else 0)
                runs.append(spec(f'wetting-{method}-theta{angle}-R{radius}-beta{beta}',
                    nx=nx,ny=ny,nz=1,dimensions=2,steps=steps,write=steps//4,
                    radius=radius,angle=angle,beta=beta,wall_gradient=method,label=label,
                    indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"',
                             f'255 "SD" "{angle}" "N/A" "N/A"'],timeout=1800))
    return runs


def circle_measure(chi,rho,width,level=0.,exclusion=.75,side=None):
    points=interface_points(chi,level,rho>0)
    points=points[points[:,1]>.5+exclusion*width]
    center_x=(chi.shape[-1]-1)/2
    if side=='left':points=points[points[:,0]<=center_x]
    if side=='right':points=points[points[:,0]>=center_x]
    if len(points)<12:raise ValueError('Insufficient contour above wall exclusion')
    center,radius,rms=fit_sphere(points)
    cosine=(.5-center[1])/radius
    if abs(cosine)>1:raise ValueError('Fitted circle does not intersect physical wall')
    return dict(angle_degrees=math.degrees(math.acos(cosine)),radius=float(radius),
                center=center.tolist(),fit_rms=rms,points=len(points),
                contour_level=level,exclusion_widths=exclusion)


def measure(result):
    s=result['parameters'];a,rho,chi,v=geometry_fields(result)
    width=WIDTHS[s['beta']]
    principal=circle_measure(chi,rho,width)
    sensitivity=[]
    for level in (-.1,0.,.1):
        for exclusion in (.5,.75,1.,1.5):
            try:entry=circle_measure(chi,rho,width,level,exclusion)
            except ValueError as e:entry=dict(contour_level=level,exclusion_widths=exclusion,error=str(e))
            sensitivity.append(entry)
    angles=[e['angle_degrees'] for e in sensitivity if 'angle_degrees' in e]
    initial=fields(result,0)
    prev=geometry_fields(result,s['steps']-s['write'])
    prior=circle_measure(prev[2],prev[1],width)
    return dict(name=result['name'],method=s['wall_gradient'],target_angle=s['angle'],
                initial_radius=s['radius'],beta=s['beta'],width_09=width,
                width_over_initial_radius=width/s['radius'],principal=principal,
                signed_angle_error=principal['angle_degrees']-s['angle'],
                sensitivity_min=min(angles),sensitivity_max=max(angles),
                left=circle_measure(chi,rho,width,side='left'),
                right=circle_measure(chi,rho,width,side='right'),
                late_angle_change=principal['angle_degrees']-prior['angle_degrees'],
                max_spurious_speed=float(np.linalg.norm(v,axis=-1).max()),
                mass_relative_drift={key:float(abs(a[key].sum()-initial[key].sum())/initial[key].sum())
                                     for key in ('r0','r1')},sensitivity=sensitivity)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',required=True);p.add_argument('--root',required=True)
    p.add_argument('--jobs',type=int,default=4)
    a=p.parse_args();study=Study(a.solver,a.root)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:results=list(pool.map(study.run,cases()))
    measurements=[measure(r) for r in results]
    summary=dict(measurements=measurements,
                 scope='Size and diffuse-interface sensitivity at fixed lattice sigma, not a formal matched-physics grid convergence study.',
                 width_calibration='Planar D2Q9 W(-0.9,+0.9); both gradients share the same bulk model.')
    (study.root/'wetting.json').write_text(json.dumps(summary,indent=2))
    for m in measurements:print(m['name'],m['principal']['angle_degrees'],m['late_angle_change'])


if __name__=='__main__':main()
