#!/usr/bin/env python3
"""Independent scalar wall transport and conservative population transfers.

This is a phase-only diagnostic. A zero-velocity recurrence is not a
hydrodynamic wetting validation.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from phase_recurrence import VELOCITIES,WEIGHTS,OPPOSITE,coefficients


def outgoing(phi,beta,angle=90,gradient='renormalized',neutral_mask=True):
    color=2*phi-1
    g=np.zeros((2,)+phi.shape)
    wall=np.zeros_like(phi,dtype=bool);wall[[0,-1]]=True
    for (dx,dy),weight in zip(VELOCITIES[1:],WEIGHTS[1:]):
        value=np.roll(color,(-dy,-dx),axis=(0,1))
        valid=np.ones_like(phi,dtype=bool)
        if dy>0:valid[-1]=False
        if dy<0:valid[0]=False
        difference=np.where(wall,value-color,value)
        g[0]+=3*weight*dx*difference*valid
        g[1]+=3*weight*dy*difference*valid
    g[0,[0,-1]]*=6/5;g[1,[0,-1]]*=2
    sine,cosine=np.sin(np.deg2rad(angle)),np.cos(np.deg2rad(angle))
    if gradient=='tangent':
        if abs(sine)<1e-8:raise ValueError('Singular geometric wall angle')
        g[0,[0,-1]]=(np.roll(color[[0,-1]],-1,axis=1)-np.roll(color[[0,-1]],1,axis=1))/2
        g[1,0]=abs(g[0,0])*cosine/sine;g[1,-1]=-abs(g[0,-1])*cosine/sine
    elif gradient!='renormalized':raise ValueError('Unknown wall gradient')
    if neutral_mask and abs(cosine)<1e-14:g[1,[0,-1]]=0
    norm=np.linalg.norm(g,axis=0)
    n=np.divide(g,norm[None],out=np.zeros_like(g),where=norm[None]>1e-9)
    for row,sign in [(0,1),(-1,-1)]:
        active=norm[row]>1e-9
        n[0,row]=np.where(active,np.where(g[0,row]>=0,1.,-1.)*sine,0)
        n[1,row]=np.where(active,sign*cosine,0)
    result=np.zeros((9,)+phi.shape)
    for q,((dx,dy),weight) in enumerate(zip(VELOCITIES,WEIGHTS)):
        projection=(dx*n[0]+dy*n[1])/np.hypot(dx,dy) if q else np.zeros_like(phi)
        result[q]=weight*(phi+beta*phi*(1-phi)*projection)
    return result


def wall_transfer(phase,total,phi):
    """Return antisymmetric relative-phase transfers on two flat walls."""
    delta=np.zeros_like(phase)
    for row,inward in [(0,1),(-1,-1)]:
        for q,(dx,dy) in enumerate(VELOCITIES):
            if dy!=inward or not dx:continue
            mirror=int(np.where(np.all(VELOCITIES==(dx,-dy),axis=1))[0][0])
            neighbor=np.roll(phi[row],dx)
            difference=np.roll(phase[mirror,row],dx)-phase[OPPOSITE[q],row]
            total_difference=np.roll(total[mirror,row],dx)-total[OPPOSITE[q],row]
            delta[q,row]=difference-.5*(phi[row]+neighbor)*total_difference
    return delta


def stream(phase):
    incoming=[]
    for q,(dx,dy) in enumerate(VELOCITIES):
        moved=np.roll(phase[q],(dy,dx),axis=(0,1))
        if dy>0:moved[0]=phase[OPPOSITE[q],0]
        if dy<0:moved[-1]=phase[OPPOSITE[q],-1]
        incoming.append(moved)
    return np.array(incoming)


def step(phi,beta,boundary='flux',angle=90,gradient='renormalized'):
    phase=outgoing(phi,beta,angle,gradient,neutral_mask=boundary!='bounce')
    if boundary=='balance':
        for q,(dx,dy) in enumerate(VELOCITIES):
            if dx and not dy:
                phase[q,[0,-1]]=WEIGHTS[q]*(phi[[0,-1]]+
                    coefficients()['tangential_factor']*(phase[q,[0,-1]]/WEIGHTS[q]-phi[[0,-1]]))
    elif boundary not in ('bounce','flux'):raise ValueError('Unknown phase boundary')
    incoming=stream(phase)
    if boundary=='flux':incoming+=wall_transfer(phase,np.broadcast_to(WEIGHTS[:,None,None],phase.shape),phi)
    return incoming.sum(axis=0)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True);p.add_argument('--steps',type=int,default=2000)
    args=p.parse_args();transient=[]
    a,b=coefficients()['bulk_diffusion'],coefficients()['bulk_segregation']
    for beta in (.35,.7,.9):
        for length in (64,65):
            phi=np.zeros((8,length));phi[:,length//4:3*length//4]=1
            bulk=phi[0].copy();maximum=0.;mass0=phi.sum()
            for _ in range(args.steps):
                color=2*bulk-1;g=(np.roll(color,-1)-np.roll(color,1))/2
                flux=np.where(abs(g)>1e-9,np.sign(g),0)*bulk*(1-bulk)
                bulk=bulk+a*(np.roll(bulk,1)+np.roll(bulk,-1)-2*bulk)+beta*b*(np.roll(flux,1)-np.roll(flux,-1))
                phi=step(phi,beta);maximum=max(maximum,float(abs(phi-bulk).max()))
            transient.append(dict(beta=beta,length=length,max_transient_error=maximum,mass_drift=float(abs(phi.sum()-mass0)/mass0)))
            if maximum>2e-12:raise AssertionError('Wall flux does not match the bulk transient')
    rng=np.random.default_rng(1800);height,length=8,65;y,x=np.indices((height,length))
    rho=1+.02*np.cos(2*np.pi*x/length)+.01*y/height
    u=.01*np.sin(2*np.pi*x/length)
    total=np.array([w*rho*(1+3*dx*u+4.5*(dx*u)**2-1.5*u*u) for (dx,dy),w in zip(VELOCITIES,WEIGHTS)])
    phi=rng.uniform(.05,.95,(height,length));phase=phi[None]*total
    phase+=outgoing(phi,.7)-WEIGHTS[:,None,None]*phi[None]
    delta=wall_transfer(phase,total,phi)
    conservative_error=float(abs(delta.sum()))
    pure_errors=[float(abs(wall_transfer(value*total,total,np.full_like(phi,value))).max()) for value in (0.,1.)]
    complement_error=float(abs((stream(phase)+delta)+(stream(total-phase)-delta)-stream(total)).max())
    if conservative_error>2e-14 or max(pure_errors)>2e-14 or complement_error>2e-14:raise AssertionError('Wall transfer conservation/pure-phase contract failed')
    result=dict(protocol='Independent D2Q9 scalar specular-forward/relative-flux control; random mixed fraction and nonuniform hydrodynamic equilibria; no flow solver.',
                neutral_transients=transient,face_component_conservation_error=conservative_error,
                pure_phase_transfer_errors=pure_errors,total_population_preservation_error=complement_error,
                scope='Zero-flow neutral recurrence and algebraic face identities; nonneutral hydrodynamic accuracy is not implied.')
    histories=[]
    for boundary in ('bounce','balance','flux'):
        phi=np.zeros((8,64));phi[:,16:48]=1;bulk=phi[0].copy();history=[]
        for time in range(1,args.steps+1):
            color=2*bulk-1;g=(np.roll(color,-1)-np.roll(color,1))/2
            flux=np.where(abs(g)>1e-9,np.sign(g),0)*bulk*(1-bulk)
            bulk=bulk+a*(np.roll(bulk,1)+np.roll(bulk,-1)-2*bulk)+.7*b*(np.roll(flux,1)-np.roll(flux,-1))
            phi=step(phi,.7,boundary)
            if time&(time-1)==0 or time==args.steps:
                history.append(dict(step=time,max_bulk_transient_error=float(abs(phi-bulk).max())))
        histories.append(dict(boundary=boundary,beta=.7,height=8,length=64,history=history))
    result['neutral_scalar_comparison']=histories
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
