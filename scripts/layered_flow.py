#!/usr/bin/env python3
"""Compare two-phase variable-viscosity channels with a sharp-interface solution."""
import argparse,json,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from benchmarks import Study,spec,geometry_fields


def cases():
    runs=[]
    for tau0,tau1 in ((.6,1.4),(1.4,.6)):
        for h in (16,32,64):
            force=1e-6*(16/h)**3
            steps=math.ceil(20*h*h/((min(tau0,tau1)-.5)/3*math.pi**2))
            steps=((steps+999)//1000)*1000
            runs.append(spec(f'layered-tau{tau0}-{tau1}-H{h}',nx=4,ny=h+2,nz=1,
                dimensions=2,tau=(tau0,tau1),steps=steps,write=steps//2,height=h,force=force,
                label=lambda x,y,z,h=h:255 if y in (0,h+1) else 75 if y<=h//2 else 76,
                indices=[f'75 "GV" "{force:.17g}" "1" "0"',f'76 "GV" "{force:.17g}" "1" "1"',
                         '255 "SD" "90" "N/A" "N/A"']))
    return runs


def measure(r):
    s=r['parameters'];a,rho,chi,v=geometry_fields(r)
    h=s['height'];mid=h/2;force=s['force'];mu0,mu1=[(t-.5)/3 for t in s['tau']]
    matrix=np.array([[0,1,0,0],[0,0,h,1],[mid,1,-mid,-1],[mu0,0,-mu1,0]],float)
    right=np.array([0,force*h*h/(2*mu1),force*mid*mid/2*(1/mu0-1/mu1),0])
    a0,b0,a1,b1=np.linalg.solve(matrix,right)
    y=np.arange(h)+.5
    ideal=np.where(y<mid,-force*y*y/(2*mu0)+a0*y+b0,-force*y*y/(2*mu1)+a1*y+b1)
    observed=v[...,0].mean(axis=(0,2))[1:-1]
    # Independently solve d/dy(mu(y) du/dy)=-F using the measured phase profile.
    # This separates diffuse-interface constitutive bias from the LB discretization.
    phi=a['r0'][...,0].mean(axis=(0,2))[1:-1]/rho.mean(axis=(0,2))[1:-1]
    harmonic=s.get('viscosity_interpolation','arithmetic')=='harmonic'
    mu=1/(phi/mu0+(1-phi)/mu1) if harmonic else phi*mu0+(1-phi)*mu1
    grid=np.linspace(0,h,h*64+1)
    mu_grid=np.interp(grid,y,mu)
    compliance=1/mu_grid
    stress_constant=force*np.trapezoid(grid*compliance,grid)/np.trapezoid(compliance,grid)
    integrand=(stress_constant-force*grid)*compliance
    continuum=np.concatenate(([0],np.cumsum((integrand[1:]+integrand[:-1])/2*np.diff(grid))))
    diffuse=np.interp(y,grid,continuum)
    return dict(name=r['name'],tau=s['tau'],height=h,
        interpolation=s.get('viscosity_interpolation','arithmetic'),
        relative_l2_error=float(np.linalg.norm(observed-ideal)/np.linalg.norm(ideal)),
        relative_diffuse_l2_error=float(np.linalg.norm(observed-diffuse)/np.linalg.norm(diffuse)),
        diffuse_continuum_profile=diffuse.tolist(),
        relative_flow_error=float(abs(observed.sum()-ideal.sum())/ideal.sum()),
        profile=observed.tolist(),analytic_profile=ideal.tolist(),
        maximum_speed=float(np.linalg.norm(v,axis=-1).max()),
        late_velocity_change=float(np.max(np.abs(v-geometry_fields(r,s['steps']-s['write'])[3]))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',required=True);p.add_argument('--root',required=True)
    p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--interpolation',choices=('arithmetic','harmonic'),default='arithmetic')
    a=p.parse_args();study=Study(a.solver,a.root)
    runs=[dict(s,name=s['name']+'-'+a.interpolation,viscosity_interpolation=a.interpolation) for s in cases()]
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:results=list(pool.map(study.run,runs))
    data=[measure(r) for r in results]
    (study.root/'layered-flow.json').write_text(json.dumps({'measurements':data},indent=2))
    for r in data:print(r['name'],r['relative_l2_error'],r['late_velocity_change'])


if __name__=='__main__':main()
