#!/usr/bin/env python3
"""Continuous-source manufactured checks of coupled Stokes/phase evolution.

The two smooth interfaces have neutral wall slopes. Analytic phase sources
and external forces prescribe their motion and a nonuniform Navier flow.
They are independent of the implemented curvature and face-flux formulas.
This verifies consistency; it is not an unforced contact-line benchmark.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time

import numpy as np

from stokes_phase_reference import StokesSlit, capillary_faces, phase_rhs


def target(x,y,t,b=.62,height=16.,length=128.):
    """Analytic logistic bridge and its continuous differential quantities."""
    ky=np.pi/height;omega=.05;amplitude=.25;translation=.0002
    eta=amplitude*np.cos(ky*y)*np.cos(omega*t)+translation*t
    etat=-amplitude*omega*np.cos(ky*y)*np.sin(omega*t)+translation
    etay=-amplitude*ky*np.sin(ky*y)*np.cos(omega*t)
    etayy=-amplitude*ky*ky*np.cos(ky*y)*np.cos(omega*t)
    offset=(x-length/2-eta+length/2)%length-length/2
    sign=-np.sign(offset);distance=length/4-abs(offset)
    phi=1/(1+np.exp(-2*b*distance));first=2*b*phi*(1-phi)
    second=2*b*(1-2*phi)*first
    # C2 joins to exact pure plateaus remove the absolute-distance cusps.
    # The joins affect only logistic tails below 5e-13, not the interface.
    lo,hi=-27.,-23.;span=hi-lo;value=1/(1+np.exp(-2*b*hi))
    slope=2*b*value*(1-value);curvature=2*b*(1-2*value)*slope
    a3=10*value-4*span*slope+.5*span**2*curvature
    a4=-15*value+7*span*slope-span**2*curvature
    a5=6*value-3*span*slope+.5*span**2*curvature
    tail=-abs(distance);r=np.maximum(0,np.minimum(1,(tail-lo)/span))
    p=a3*r**3+a4*r**4+a5*r**5
    p1=(3*a3*r*r+4*a4*r**3+5*a5*r**4)/span
    p2=(6*a3*r+12*a4*r*r+20*a5*r**3)/span**2
    active=tail<hi;upper=distance>0
    phi=np.where(active,np.where(upper,1-p,p),phi)
    first=np.where(active,p1,first);second=np.where(active,np.where(upper,-p2,p2),second)
    norm=np.sqrt(1+etay*etay);product=2*b*phi*(1-phi)
    divergence=second*norm**2-sign*first*etayy-2*b*(1-2*phi)*first*norm+sign*product*etayy/norm**3
    return dict(phi=phi,dx=sign*first,dy=-sign*first*etay,dt=-sign*first*etat,
        phase_divergence=divergence,capillary_x=etayy/norm**3*first,
        capillary_y=-etayy*etay/norm**3*first,interface_derivative=first)


def flow(x,y,slip,height=16.,length=128.):
    k=2*np.pi/length;nu=1/6;scale=4e-8;pressure=1e-5;ky=np.pi/height
    c=2*slip*height**2/(height+2*slip)
    a=y*y*(height-y)**2+c*y*(height-y)
    a1=2*height*height*y-6*height*y*y+4*y**3+c*(height-2*y)
    a2=2*height*height-12*height*y+12*y*y-2*c
    a3=-12*height+24*y
    return dict(u=scale*a1*np.sin(k*x),v=-scale*k*a*np.cos(k*x),
        fx=-nu*scale*(a3-k*k*a1)*np.sin(k*x)-pressure*k*np.sin(k*x)*np.cos(ky*y),
        fy=nu*scale*k*(a2-k*k*a)*np.cos(k*x)-pressure*ky*np.cos(k*x)*np.sin(ky*y))


def run(ny,slip,folder,end=16.,dt_factor=.125,sigma=.005):
    folder.mkdir(parents=True,exist_ok=False);solver=StokesSlit(ny,nx=8*ny,slip=slip)
    h=solver.h;b=.62;nu=1/6;diffusion=1/6
    xc=(np.arange(solver.nx)+.5)[None]*h;xf=np.arange(solver.nx)[None]*h
    yc=(np.arange(ny)+.5)[:,None]*h;yf=np.arange(ny+1)[:,None]*h
    cell_flow=flow(xc,yc,slip);x_flow=flow(xf,yc,slip);y_flow=flow(xc,yf,slip)
    phi=target(xc,yc,0.,b)['phi'];dt=dt_factor*h*h;steps=int(np.ceil(end/dt));dt=end/steps
    mass0=phi.sum()*h*h;started=time.monotonic();history=[]
    def rhs(field,clock):
        exact=target(xc,yc,clock,b);x_exact=target(xf,yc,clock,b);y_exact=target(xc,yf,clock,b)
        fx,fy,_=capillary_faces(field,h,sigma,b)
        fx+=x_flow['fx']-sigma*x_exact['capillary_x']
        fy+=y_flow['fy']-sigma*y_exact['capillary_y']
        u,v=solver.solve(fx,fy)
        source=exact['dt']+cell_flow['u']*exact['dx']+cell_flow['v']*exact['dy']-diffusion*exact['phase_divergence']
        rate=phase_rhs(field,u,v,h,b,transformed=True)+source
        return rate,u,v,float(source.sum()*h*h)
    budget=0.
    for step in range(steps+1):
        clock=step*dt;rate,u,v,source_integral=rhs(phi,clock)
        if step in (0,steps//2,steps):
            exact=target(xc,yc,clock,b);expected=exact['phi']
            pe=float(np.linalg.norm(phi-expected)/np.linalg.norm(expected))
            displacement=float(np.linalg.norm(phi-expected)/np.linalg.norm(exact['interface_derivative']))
            ve=float(np.sqrt(np.sum((u-x_flow['u'])**2)+np.sum((v-y_flow['v'])**2))/
                np.sqrt(np.sum(x_flow['u']**2)+np.sum(y_flow['v']**2)))
            history.append(dict(time=clock,relative_phase_error=pe,interface_displacement_error=displacement,
                relative_velocity_error=ve,max_divergence=float(abs(solver.divergence(u,v)).max()),
                mass_budget_residual=float(phi.sum()*h*h-mass0-budget),phase_range=[float(phi.min()),float(phi.max())]))
            np.savez_compressed(folder/f'field.{step:010d}.npz',phi=phi,u=u,v=v,exact_phi=expected,
                exact_u=x_flow['u'],exact_v=y_flow['v'],time=clock,h=h)
        if step==steps:break
        trial=phi+dt*rate;second,_,_,source_next=rhs(trial,clock+dt)
        phi+=dt*(rate+second)/2;budget+=dt*(source_integral+source_next)/2
        if not np.isfinite(phi).all() or phi.min()<-2e-11 or phi.max()>1+2e-11:
            raise ArithmeticError(f'Manufactured state outside phase bounds at step {step+1}')
    result=dict(parameters=dict(ny=ny,nx=8*ny,h=h,b=b,slip=slip,nu=nu,diffusion=diffusion,sigma=sigma,end=end,dt=dt,dt_factor=dt_factor),
        history=history,seconds=time.monotonic()-started,
        scope='Coupled continuous-source manufactured consistency; not unforced contact-line accuracy or rock validation.')
    (folder/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(ny=ny,slip=slip,dt_factor=dt_factor,final=history[-1])),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True);parser.add_argument('--heights',type=int,nargs='+',default=[16,32,64])
    parser.add_argument('--slips',type=float,nargs='+',default=[1.,4.]);parser.add_argument('--end',type=float,default=16.)
    args=parser.parse_args();root=Path(args.root);root.mkdir(parents=True,exist_ok=False)
    source=root/'source';source.mkdir();hashes={}
    for p in Path(__file__).parent.glob('*.py'):
        shutil.copy2(p,source/p.name);hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    (root/'source-manifest.json').write_text(json.dumps(hashes,indent=2)+'\n');rows=[];orders=[]
    for slip in args.slips:
        group=[run(h,slip,root/f'H{h}-slip{slip:g}',end=args.end) for h in args.heights];rows.extend(group)
        for key in ('relative_phase_error','interface_displacement_error','relative_velocity_error'):
            errors=[r['history'][-1][key] for r in group]
            orders.append(dict(slip=slip,observable=key,errors=errors,
                orders=[float(np.log(a/b)/np.log(hb/ha)) for a,b,ha,hb in zip(errors,errors[1:],args.heights,args.heights[1:])]))
    half=run(32,1.,root/'H32-slip1-halfstep',end=args.end,dt_factor=.0625);rows.append(half)
    result=dict(runs=rows,spatial_orders=orders,source_sha256=hashes,
        scope='Independent continuous-source coupled verification; observed orders are reported, not assumed.')
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(orders),flush=True)


if __name__=='__main__':main()
