#!/usr/bin/env python3
"""Analytic and manufactured checks of the independent reference, not LBM tests."""
import argparse
import json
from pathlib import Path

import numpy as np

from front_mobility import match_beta, WIDTH_REFERENCE
from stokes_phase_reference import (StokesSlit, capillary_faces, phase_rhs, planar_bridge,
    logistic_bridge, match_logistic_width, observation)


def check():
    rows=[];manufactured=[]
    for ny in (16,32,64):
        s=StokesSlit(ny);h=s.h;nx=s.nx;H=16.;L=64.;nu=1/6
        fx=np.ones((ny,nx))*1e-7;fy=np.zeros((ny+1,nx));u,v=s.solve(fx,fy)
        y=(np.arange(ny)+.5)*h
        exact=1e-7*(y*(H-y)+h*h/4)/(2*nu)
        pure_error=float(abs(u-exact[:,None]).max()/exact.max())
        p=np.cos(2*np.pi*(np.arange(nx)+.5)[None]/nx)*np.cos(np.pi*(np.arange(ny)+.5)[:,None]/ny)
        fx=(p-np.roll(p,1,axis=1))/h;fy[1:-1]=np.diff(p,axis=0)/h
        ug,vg=s.solve(fx,fy);gradient_speed=float(max(abs(ug).max(),abs(vg).max()))
        k=2*np.pi/L;xface=np.arange(nx)*h;xcell=(np.arange(nx)+.5)*h
        Y=y[:,None];A=Y**2*(H-Y)**2;A1=2*H*H*Y-6*H*Y*Y+4*Y**3;A3=-12*H+24*Y
        fx=-nu*(A3-k*k*A1)*np.sin(k*xface)[None]
        yf=np.arange(ny+1)[:,None]*h;Af=yf*yf*(H-yf)**2;A2=2*H*H-12*H*yf+12*yf*yf
        fy=nu*k*(A2-k*k*Af)*np.cos(k*xcell)[None]
        um,vm=s.solve(fx,fy)
        targetu=A1*np.sin(k*xface)[None];targetv=-k*Af*np.cos(k*xcell)[None]
        error=float(np.sqrt(np.sum((um-targetu)**2)+np.sum((vm-targetv)**2))/np.sqrt(np.sum(targetu**2)+np.sum(targetv**2)))
        manufactured.append(error)
        beta=match_beta(WIDTH_REFERENCE/h);phi=np.broadcast_to(planar_bridge(nx,beta),(ny,nx)).copy()
        f0x,f0y,kappa=capillary_faces(phi,h,.005)
        rate=phase_rhs(phi,np.zeros_like(um),np.zeros_like(vm),h,beta/h)
        planar_error=float(abs(rate).max());mass_rate=float(abs(rate.sum()*h*h))
        row=dict(ny=ny,h=h,pure_relative_profile_error=pure_error,
            pure_relative_continuum_mean_error=float(u.mean()/(1e-7*H*H/(12*nu))-1),
            gradient_force_max_speed=gradient_speed,manufactured_relative_velocity_error=error,
            max_manufactured_divergence=float(abs(s.divergence(um,vm)).max()),
            planar_phase_residual=planar_error,planar_capillary_force=float(max(abs(f0x).max(),abs(f0y).max())),
            planar_mass_rate=mass_rate)
        if pure_error>2e-12 or gradient_speed>2e-12 or row['max_manufactured_divergence']>2e-10 or planar_error>2e-9 or row['planar_capillary_force']>2e-14 or mass_rate>2e-13:
            raise AssertionError(row)
        rows.append(row)
    orders=[float(np.log(a/b)/np.log(2)) for a,b in zip(manufactured,manufactured[1:])]
    if min(orders)<1.8:raise AssertionError('Manufactured Stokes flow does not approach its continuum solution quadratically')
    slip_rows=[]
    for ny in (16,32,64):
        for slip in (.5,1.,2.):
            s=StokesSlit(ny,slip=slip);h=s.h;fx=np.full((ny,s.nx),1e-7);fy=np.zeros((ny+1,s.nx));u,v=s.solve(fx,fy)
            y=(np.arange(ny)+.5)*h;target=3e-7*(y*(16-y)+16*slip+h*h/4)
            error=float(abs(u-target[:,None]).max()/target.max())
            if error>2e-12:raise AssertionError('Independent Navier hydraulic formula failed')
            slip_rows.append(dict(ny=ny,slip=slip,relative_profile_error=error))
    scalar=[]
    for ny in (16,32):
        h=16/ny;nx=4*ny;b=match_logistic_width(WIDTH_REFERENCE,h)
        phi=np.broadcast_to(logistic_bridge(nx,h,b),(ny,nx)).copy();u=np.full_like(phi,.001);v=np.zeros((ny+1,nx))
        rate=phase_rhs(phi,u*0,v,h,b,transformed=True)
        if abs(rate).max()>1e-10:raise AssertionError('Transformed planar equilibrium is not stationary')
        first=observation(phi,u,v,h,0.)['centroid_modulo'];mass=phi.sum();dt=h*h
        for _ in range(round(100/dt)):
            a=phase_rhs(phi,u,v,h,b,transformed=True)
            z=phi+dt*a;phi+=dt*(a+phase_rhs(z,u,v,h,b,transformed=True))/2
        last=observation(phi,u,v,h,100.)['centroid_modulo']
        speed=(((last-first)+32)%64-32)/100
        if abs(speed/.001-1)>.003 or abs(phi.sum()-mass)/mass>1e-12:raise AssertionError('Prescribed scalar translation fails')
        scalar.append(dict(ny=ny,relative_translation_error=speed/.001-1,relative_mass_change=(phi.sum()-mass)/mass))
    return dict(controls=rows,manufactured_orders=orders,navier_controls=slip_rows,scalar_translation=scalar,
        scope='Independent hydraulic solve verified by analytic/manufactured fields and Navier profiles; scalar planar equilibrium and prescribed translation verified. Coupled dynamic accuracy needs separate resolution/time controls.')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    args=parser.parse_args();result=check();p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
