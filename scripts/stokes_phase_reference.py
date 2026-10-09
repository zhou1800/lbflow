#!/usr/bin/env python3
"""Staggered Stokes / conservative phase-field references and kinetic ablations.

Periodic flow direction, neutral flat walls, constant viscosity, unit density.
The hydrodynamic solver uses no LBM collision or production force code.
An optional phase-stencil ablation retains discrete projection transport.
The sharpening equation is not a free-energy/chemical-potential model.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time

import numpy as np

from front_mobility import match_beta, WIDTH_REFERENCE
from phase_recurrence import VELOCITIES, WEIGHTS

LOADED_SOURCE_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


class StokesSlit:
    """MAC velocity from a Fourier/biharmonic stream function, including mean flux."""
    def __init__(self,ny,nx=None,height=16.,nu=1/6,slip=0.,fast=False):
        self.ny=ny;self.nx=nx or 4*ny;self.h=height/ny;self.nu=nu;self.fast=fast
        h=self.h;m=ny-1
        frequencies=np.arange(self.nx//2+1)
        lam=-4*np.sin(np.pi*frequencies/self.nx)**2/h**2
        wall=2/(1+2*slip/h)/h**4
        if fast:
            # The same Dirichlet Laplacian plus two boundary rank-one terms.
            modes=np.arange(1,ny)
            ly=-4*np.sin(np.pi*modes/(2*ny))**2/h**2
            self.sine_weights=1/(ly[:,None]+lam[None])**2
            first=np.sqrt(2/ny)*np.sin(np.pi*modes/ny)
            self.wall_modes=[np.sqrt(2)*first*(modes%2==parity) for parity in (1,0)]
            self.wall_denominators=[1/wall+(q[:,None]**2*self.sine_weights).sum(axis=0) for q in self.wall_modes]
        else:
            lap=np.diag(np.full(m,-2/h**2))+np.diag(np.full(m-1,1/h**2),1)+np.diag(np.full(m-1,1/h**2),-1)
            operators=lap[None]+lam[:,None,None]*np.eye(m)[None]
            biharmonic=np.einsum('kij,kjl->kil',operators,operators)
            biharmonic[:,0,0]+=wall;biharmonic[:,-1,-1]+=wall
            self.inverse=np.linalg.inv(nu*biharmonic)
        mean_lap=np.diag(np.full(ny,2/h**2))+np.diag(np.full(ny-1,-1/h**2),1)+np.diag(np.full(ny-1,-1/h**2),-1)
        ghost_ratio=(2*slip-h)/(2*slip+h)
        mean_lap[0,0]-=ghost_ratio/h**2;mean_lap[-1,-1]-=ghost_ratio/h**2
        self.mean_inverse=np.linalg.inv(nu*mean_lap)

    def solve(self,fx,fy):
        """fx at x faces/cell y; fy at y faces/cell x, with zero wall velocity."""
        h=self.h;curl=(fy[1:-1]-np.roll(fy[1:-1],1,axis=1))/h-(fx[1:]-fx[:-1])/h
        transformed=np.fft.rfft(curl,axis=1)
        transformed[:,0]=0  # Mean flow is solved independently below.
        if self.fast:
            from scipy.fft import dst
            interior=dst(transformed,type=1,axis=0,norm='ortho')*self.sine_weights
            for q,denominator in zip(self.wall_modes,self.wall_denominators):
                coefficient=(q[:,None]*interior).sum(axis=0)/denominator
                interior-=self.sine_weights*q[:,None]*coefficient[None]
            interior=dst(interior,type=1,axis=0,norm='ortho')/self.nu
        else:interior=np.einsum('kij,jk->ik',self.inverse,transformed,optimize=False)
        psi=np.zeros((self.ny+1,self.nx));psi[1:-1]=np.fft.irfft(interior,n=self.nx,axis=1)
        u=np.diff(psi,axis=0)/h
        u+=self.mean_inverse@fx.mean(axis=1)[:,None]
        v=-(np.roll(psi,-1,axis=1)-psi)/h
        return u,v

    def divergence(self,u,v):
        return (np.roll(u,-1,axis=1)-u)/self.h+np.diff(v,axis=0)/self.h


def planar_bridge(nx,beta):
    """Independent symmetric zero-flux face recurrence, two interfaces."""
    high=(1+(np.sqrt(1+beta*beta)-1)/beta)/2
    low=[1-high]
    for _ in range(nx//4-1):
        a=low[-1];c=(1-beta)*a+beta*a*a
        low.append(2*c/(1+beta+np.sqrt((1+beta)**2-4*beta*c)))
    low=np.array(low);profile=np.concatenate((low[::-1],1-low,1-low[::-1],low))
    if len(profile)!=nx:raise ValueError('Bridge length must be divisible by four')
    return profile


def logistic_bridge(nx,h,b):
    x=(np.arange(nx)+.5)*h;length=nx*h
    distance=length/4-abs(x-length/2)
    return (1+np.tanh(b*distance))/2


def match_logistic_width(width,h):
    # Same measured crossing convention, independent analytic tanh profile.
    lo,hi=.1,2.
    x=(np.arange(100)+.5)*h
    for _ in range(56):
        b=(lo+hi)/2;fraction=(1-np.tanh(b*x))/2
        measured=2*np.interp(.05,fraction[::-1],x[::-1])
        if measured>width:lo=b
        else:hi=b
    return (lo+hi)/2


def centered_gradient(phi,h):
    gx=(np.roll(phi,-1,axis=1)-np.roll(phi,1,axis=1))/(2*h)
    extended=np.vstack((phi[0],phi,phi[-1]))
    gy=(extended[2:]-extended[:-2])/(2*h)
    return gx,gy


def normals(phi,h,cutoff=1e-9):
    gx,gy=centered_gradient(phi,h);norm=np.hypot(gx,gy)
    nx=np.divide(gx,norm,out=np.zeros_like(gx),where=norm>cutoff)
    ny=np.divide(gy,norm,out=np.zeros_like(gy),where=norm>cutoff)
    return nx,ny,norm


def distance_field(phi,b,epsilon=1e-12):
    # Only reconstruction is regularized; the evolving fraction is not clipped.
    interior=np.maximum(epsilon,np.minimum(1-epsilon,phi))
    return np.log(interior/(1-interior))/(2*b)


def capillary_faces(phi,h,sigma,b=None):
    geometry=phi if b is None else distance_field(phi,b)
    nx,ny,norm=normals(geometry,h)
    nx_x=(np.roll(nx,-1,axis=1)-np.roll(nx,1,axis=1))/(2*h)
    ny_x=(np.roll(ny,-1,axis=1)-np.roll(ny,1,axis=1))/(2*h)
    ex=np.vstack((nx[0],nx,nx[-1]));ey=np.vstack((-ny[0],ny,-ny[-1]))
    nx_y=(ex[2:]-ex[:-2])/(2*h);ny_y=(ey[2:]-ey[:-2])/(2*h)
    curvature=(nx*nx-1)*nx_x+nx*ny*(nx_y+ny_x)+(ny*ny-1)*ny_y
    curvature=np.where(norm>1e-9,curvature,0.)
    kx=(curvature+np.roll(curvature,1,axis=1))/2
    fx=sigma*kx*(phi-np.roll(phi,1,axis=1))/h
    fy=np.zeros((phi.shape[0]+1,phi.shape[1]))
    fy[1:-1]=sigma*(curvature[1:]+curvature[:-1])/2*np.diff(phi,axis=0)/h
    return fx,fy,curvature


def phase_rhs(phi,u,v,h,b,diffusion=1/6,transformed=False):
    """Conservative centered advective flux and sharpening/diffusion face flux."""
    geometry=distance_field(phi,b) if transformed else phi
    left=np.roll(phi,1,axis=1);gx=(geometry-np.roll(geometry,1,axis=1))/h
    _,gy=centered_gradient(geometry,h)
    tangent=(gy+np.roll(gy,1,axis=1))/2;norm=np.hypot(gx,tangent)
    normal_x=np.divide(gx,norm,out=np.zeros_like(gx),where=norm>1e-9)
    product=phi*(1-phi)
    if transformed:
        jump=geometry-np.roll(geometry,1,axis=1)
        chain=np.divide(phi-left,jump,out=b*(product+np.roll(product,1,axis=1)),where=abs(jump)>1e-13)
        qx=diffusion*chain*(gx-normal_x)
    else:qx=diffusion*(gx-b*(product+np.roll(product,1,axis=1))*normal_x)
    fluxx=u*(phi+left)/2-qx
    gx_center,_=centered_gradient(geometry,h);gy_face=np.diff(geometry,axis=0)/h
    tangent=(gx_center[1:]+gx_center[:-1])/2;norm=np.hypot(gy_face,tangent)
    normal_y=np.divide(gy_face,norm,out=np.zeros_like(gy_face),where=norm>1e-9)
    fluxy=np.zeros_like(v)
    if transformed:
        jump=np.diff(geometry,axis=0)
        chain=np.divide(np.diff(phi,axis=0),jump,out=b*(product[1:]+product[:-1]),where=abs(jump)>1e-13)
        qy=diffusion*chain*(gy_face-normal_y)
    else:qy=diffusion*(gy_face-b*(product[1:]+product[:-1])*normal_y)
    fluxy[1:-1]=v[1:-1]*(phi[1:]+phi[:-1])/2-qy
    return -(np.roll(fluxx,-1,axis=1)-fluxx)/h-np.diff(fluxy,axis=0)/h


def projection_rest(phi,beta,cutoff=1e-9):
    """Zero-flow phase stencil for a kinetic ablation, not an independent PDE."""
    extra=np.vstack((phi[0],phi,phi[-1]));g=np.zeros((2,)+phi.shape)
    for (ex,ey),w in zip(VELOCITIES,WEIGHTS):
        neighbor=np.roll(extra,-ex,axis=1)[1+ey:1+ey+len(phi)]
        g[0]+=6*w*ex*neighbor;g[1]+=6*w*ey*neighbor
    norm=np.linalg.norm(g,axis=0);normal=np.divide(g,norm[None],out=np.zeros_like(g),where=norm[None]>cutoff)
    outgoing=np.array([w*(phi+beta*phi*(1-phi)*(ex*normal[0]+ey*normal[1])) for (ex,ey),w in zip(VELOCITIES,WEIGHTS)])
    result=np.zeros_like(phi)
    for q,(ex,ey) in enumerate(VELOCITIES):
        incoming=np.roll(outgoing[q],(ey,ex),axis=(0,1))
        if ey:
            mirror=int(np.where(np.all(VELOCITIES==(ex,-ey),axis=1))[0][0])
            row=0 if ey>0 else -1
            incoming[row]=np.roll(outgoing[mirror,row],ex)
        result+=incoming
    return result


def transverse_relaxation_rhs(phi,h,strength):
    """Selective neutral-wall biharmonic diagnostic; not a physical phase law."""
    if not np.isfinite(h) or h<=0 or not np.isfinite(strength) or not 0<=strength<=1:
        raise ValueError('Use positive grid spacing and diagnostic strength in [0,1]')
    extended=np.vstack((phi[0],phi,phi[-1]))
    laplacian=(extended[2:]-2*phi+extended[:-2])/h**2
    extended=np.vstack((laplacian[0],laplacian,laplacian[-1]))
    biharmonic=(extended[2:]-2*laplacian+extended[:-2])/h**2
    return -strength*h*h*biharmonic/24


def observation(phi,u,v,h,clock):
    nx=phi.shape[1];length=nx*h
    z=np.sum(phi*np.exp(2j*np.pi*(np.arange(nx)+.5)/nx)[None])
    centroid=float(np.angle(z)*length/(2*np.pi))
    return dict(time=clock,mean_speed=float(u.mean()),centroid_modulo=centroid,
        max_speed=float(max(abs(u).max(),abs(v).max())),phase_mass=float(phi.sum()*h*h),
        phase_range=[float(phi.min()),float(phi.max())])


def evolve(ny,force,end=2000.,dt_factor=1.,sigma=.005,width=WIDTH_REFERENCE,save=None,write=250.,transformed=False,slip=0.,phase_stencil='fv',initial_profile='auto',transverse_relaxation=0.,fast_stokes=False):
    if ny<4 or ny%2 or not all(np.isfinite(x) for x in (force,end,dt_factor,sigma,width,write,slip,transverse_relaxation)):
        raise ValueError('Use an even height >=4 and finite parameters')
    if min(end,dt_factor,width,write)<=0 or min(sigma,slip)<0:
        raise ValueError('Duration, time factor, width and output interval must be positive; sigma and slip nonnegative')
    if not 0<=transverse_relaxation<=1 or (transverse_relaxation and (phase_stencil!='fv' or not transformed)):
        raise ValueError('Transverse diagnostic requires transformed FV and strength in [0,1]')
    solver=StokesSlit(ny,slip=slip,fast=fast_stokes);h=solver.h;beta=match_beta(width/h);b=beta/h
    if transformed and phase_stencil=='fv':b=match_logistic_width(width,h)
    profile=logistic_bridge(solver.nx,h,b) if transformed and phase_stencil=='fv' else planar_bridge(solver.nx,beta)
    if initial_profile=='projection':profile=planar_bridge(solver.nx,beta)
    elif initial_profile=='logistic':profile=logistic_bridge(solver.nx,h,match_logistic_width(width,h))
    phi=np.broadcast_to(profile,(ny,solver.nx)).copy()
    dt=dt_factor*h*h;steps=int(np.ceil(end/dt));dt=end/steps
    history=[];mass0=phi.sum()*h*h;started=time.monotonic();next_output=0.
    def rhs(field):
        fx,fy,k=capillary_faces(field,h,sigma,b if transformed else None);fx+=force
        u,v=solver.solve(fx,fy)
        if phase_stencil=='projection':
            rate=phase_rhs(field,u,v,h,b,diffusion=0.)+(projection_rest(field,beta)-field)/h**2
        else:rate=phase_rhs(field,u,v,h,b,transformed=transformed)
        if transverse_relaxation:rate+=transverse_relaxation_rhs(field,h,transverse_relaxation)
        return rate,u,v
    for step in range(steps+1):
        rate,u,v=rhs(phi);clock=step*dt
        if step==steps or clock>=next_output-1e-9:
            row=observation(phi,u,v,h,clock)
            row['max_divergence']=float(abs(solver.divergence(u,v)).max())
            row['relative_mass_change']=float((row['phase_mass']-mass0)/mass0)
            history.append(row);next_output+=write
            if save:np.savez_compressed(save/f'field.{step:010d}.npz',phi=phi,u=u,v=v,time=clock,h=h)
            print(json.dumps(dict(ny=ny,force=force,time=clock,mean_speed=row['mean_speed'],seconds=round(time.monotonic()-started))),flush=True)
        if step==steps:break
        trial=phi+dt*rate;next_rate,_,_=rhs(trial)
        phi=phi+dt*(rate+next_rate)/2
        if not np.isfinite(phi).all() or phi.min() < -2e-11 or phi.max()>1+2e-11:
            raise ArithmeticError(f'Invalid phase fraction at step {step+1}: {phi.min()}, {phi.max()}')
    positions=np.unwrap([row['centroid_modulo']*2*np.pi/(solver.nx*h) for row in history])*solver.nx*h/(2*np.pi)
    for row,x in zip(history,positions):row['centroid_unwrapped']=float(x)
    first,last=history[len(history)//2],history[-1]
    speed=(last['centroid_unwrapped']-first['centroid_unwrapped'])/(last['time']-first['time']) if end else 0.
    return dict(parameters=dict(ny=ny,nx=solver.nx,h=h,dt=dt,end=end,force=force,sigma=sigma,nu=1/6,diffusion=1/6,beta=None if transformed and phase_stencil=='fv' else beta,b=b,width_target=width,angle=90,transformed=transformed,slip=slip,phase_stencil=phase_stencil,initial_profile=initial_profile,transverse_relaxation=transverse_relaxation,fast_stokes=fast_stokes),
        history=history,late_centroid_speed=speed,seconds=time.monotonic()-started,
        source_sha256=LOADED_SOURCE_SHA256,
        scope=('Independent hydrodynamics with a retained discrete phase stencil; kinetic ablation, not independent continuum validation.' if phase_stencil=='projection' else 'Selective transverse fourth-order phase intervention with theoretical coefficient; numerical mechanism diagnostic, not a physical phase law or the full projection stencil.' if transverse_relaxation else 'Independent periodic neutral-wall creeping-flow/finite-volume phase reference. No physical slip calibration or nonneutral wall equivalence implied.'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ny',type=int,default=16);parser.add_argument('--force',type=float,default=1e-7)
    parser.add_argument('--end',type=float,default=2000.);parser.add_argument('--dt-factor',type=float,default=1.)
    parser.add_argument('--sigma',type=float,default=.005);parser.add_argument('--root',required=True)
    parser.add_argument('--width',type=float,default=WIDTH_REFERENCE,help='Sampled physical W09; independently specified, never fitted to flow.')
    parser.add_argument('--transformed',action='store_true')
    parser.add_argument('--slip',type=float,default=0.)
    parser.add_argument('--phase-stencil',choices=['fv','projection'],default='fv')
    parser.add_argument('--initial-profile',choices=['auto','projection','logistic'],default='auto')
    parser.add_argument('--transverse-relaxation',type=float,default=0.,help='Selective numerical diagnostic strength [0,1]; coefficient h^2/24, default off.')
    parser.add_argument('--fast-stokes',action='store_true',help='Equivalent sine-basis/rank-two inversion; default keeps the original dense solve.')
    args=parser.parse_args();root=Path(args.root);root.mkdir(parents=True,exist_ok=False)
    source=root/'source';source.mkdir()
    for name in ('stokes_phase_reference.py','front_mobility.py','phase_recurrence.py'):
        shutil.copyfile(Path(__file__).parent/name,source/name)
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir()}
    if hashes['stokes_phase_reference.py']!=LOADED_SOURCE_SHA256:
        raise RuntimeError('Reference source changed between import and freezing')
    result=evolve(args.ny,args.force,args.end,args.dt_factor,args.sigma,width=args.width,save=root,transformed=args.transformed,slip=args.slip,phase_stencil=args.phase_stencil,initial_profile=args.initial_profile,transverse_relaxation=args.transverse_relaxation,fast_stokes=args.fast_stokes)
    result['source_files_sha256']=hashes
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
