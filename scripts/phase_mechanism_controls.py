#!/usr/bin/env python3
"""Planar-profile resistance predictions and independently defined operator checks.

The reduced model keeps a translated planar normal profile. Its capillary
response comes from the reference Stokes solve and does not use front data.
Failure of that reduction is reported; it is not fitted into agreement.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from front_mobility import WIDTH_REFERENCE
from stokes_phase_reference import (StokesSlit,capillary_faces,phase_rhs,
    logistic_bridge,match_logistic_width,transverse_relaxation_rhs)


def operator_checks():
    rng=np.random.default_rng(1800);rows=[]
    for H in (16,32,48):
        h=16/H;y=(np.arange(H)+.5)*h
        for m in (1,2,3):
            mode=np.cos(m*np.pi*y/16)[:,None]*np.ones((1,4*H))
            rate=transverse_relaxation_rhs(mode,h,1)
            gamma=h*h/24*(4/h**2*np.sin(m*np.pi/(2*H))**2)**2
            error=float(abs(rate+gamma*mode).max())
            if error>3e-13:raise AssertionError('Added operator does not match its independent Fourier eigenvalue')
            rows.append(dict(height=H,mode=m,semidiscrete_rate=gamma,max_operator_error=error))
        phi=rng.random((H,4*H));parts=[transverse_relaxation_rhs(phi,h,1,k) for k in ('uniform','wall','bulk')]
        split=float(abs(parts[0]-parts[1]-parts[2]).max())
        mass=max(float(abs(p.sum()*h*h)) for p in parts);work=[float((phi*p).sum()*h*h) for p in parts]
        if split>1e-13 or mass>3e-12 or max(work)>=0:raise AssertionError('Localized operator violates its conservation/negative-work identity')
        planar=np.broadcast_to(logistic_bridge(4*H,h,match_logistic_width(WIDTH_REFERENCE,h)),phi.shape)
        if any(np.any(transverse_relaxation_rhs(planar,h,1,k)) for k in ('uniform','wall','bulk')):
            raise AssertionError('Selective operator changes the independent planar profile')
        rows.append(dict(height=H,max_partition_error=split,max_mass_rate=mass,quadratic_work=work))
    return rows


def reduced(H,slip,sigma=.005,force=1e-7,epsilon=1e-4):
    h=16/H;nx=4*H;b=match_logistic_width(WIDTH_REFERENCE,h);solver=StokesSlit(H,slip=slip,fast=True)
    x=(np.arange(nx)+.5)*h;p=logistic_bridge(nx,h,b)
    derivative=-np.sign(x-32)*2*b*p*(1-p);weight=h*derivative[:nx//2].sum()
    phi=np.broadcast_to(p,(H,nx)).copy()
    def project(rate):return -h*rate[:,:nx//2].sum(axis=1)/weight
    def capillary(eta):
        delta=-eta[:,None]*derivative[None]
        plus=capillary_faces(phi+epsilon*delta,h,sigma,b);minus=capillary_faces(phi-epsilon*delta,h,sigma,b)
        fx=(plus[0]-minus[0])/(2*epsilon);fy=(plus[1]-minus[1])/(2*epsilon)
        u,v=solver.solve(fx,fy)
        return project(phase_rhs(phi,u,v,h,b,diffusion=0.)),u,v
    columns=[capillary(e)[0] for e in np.eye(H)];matrix=np.array(columns).T
    rigid=float(abs(matrix@np.ones(H)).max())
    if rigid>1e-8:raise AssertionError('Reduced capillary operator does not preserve rigid displacement')
    # Remove only the verified roundoff/finite-difference residue of that symmetry.
    matrix-=matrix.sum(axis=1)[:,None]/H
    u0,v0=solver.solve(np.full((H,nx),force),np.zeros((H+1,nx)))
    drive=project(phase_rhs(phi,u0,v0,h,b,diffusion=0.))
    rows=[]
    for a in (0.,.5,1.):
        numerical=np.array([transverse_relaxation_rhs(e[:,None],h,a)[:,0] for e in np.eye(H)]).T
        op=matrix+numerical;system=np.zeros((H+1,H+1));system[:H,:H]=op;system[:H,-1]=-1;system[-1,:H]=1
        rhs=np.r_[-drive,0.];answer=np.linalg.solve(system,rhs);eta=answer[:H];speed=answer[-1]
        _,u,v=capillary(eta)
        residual=float(abs(system@answer-rhs).max())
        _,_,null_vectors=np.linalg.svd(op.T)
        left=null_vectors[-1];left/=left.sum()
        left_residual=float(abs(left@op).max())
        mobility_error=float(abs(float(left@drive)-speed)/abs(speed))
        eigenvalues=np.linalg.eigvals(op)
        positive=float(eigenvalues.real.max())
        if left_residual>1e-11 or mobility_error>1e-9 or positive>1e-10:
            raise AssertionError('Reduced translation/null-vector consistency or stability failed')
        rows.append(dict(height=H,slip=slip,sigma=sigma,force=force,strength=a,
            predicted_phase_speed=float(speed),predicted_mean_speed=float((u+u0).mean()),
            max_displacement=float(abs(eta).max()),linear_system_residual=residual,
            rigid_mode_residual_before_projection=rigid,epsilon=epsilon,shape=eta.tolist(),
            left_null_residual=left_residual,left_null_phase_speed_relative_error=mobility_error,
            maximum_eigenvalue_real_part=positive,left_null_vector=left.tolist()))
    rates=[]
    for m in (1,2,3):
        mode=np.cos(m*np.pi*(np.arange(H)+.5)/H)
        cap_rate=float(-(mode@(matrix@mode))/(mode@mode))
        num_rate=h*h/24*(4/h**2*np.sin(m*np.pi/(2*H))**2)**2
        rates.append(dict(height=H,slip=slip,sigma=sigma,mode=m,
            initial_projected_capillary_rate=cap_rate,numerical_rate=num_rate,ratio=num_rate/cap_rate,
            scope='Initial projected response on a held planar profile; not a full coupled eigenvalue or universal error criterion.'))
    return rows,rates


def calculate():
    rows=[];rates=[]
    for H,slip,sigma in [(h,s,.005) for h in (16,32,48) for s in (1.,4.)]+[(h,1.,s) for h in (16,32) for s in (.0025,.01)]:
        r,k=reduced(H,slip,sigma);rows.extend(r);rates.extend(k)
    repeat,_=reduced(32,1.,epsilon=2e-4)
    return dict(operator_checks=operator_checks(),reduced_predictions=rows,initial_rate_diagnostics=rates,
        finite_difference_repeat=repeat,no_front_measurements_used=True,
        scope='Selective-operator verification and an explicitly approximate held-normal-profile resistance reduction. Its accuracy is measured separately, never assumed.')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',required=True)
    args=parser.parse_args();root=Path(args.root);root.mkdir(parents=True,exist_ok=False)
    source=root/'source';source.mkdir();hashes={}
    for p in Path(__file__).parent.glob('*.py'):
        shutil.copy2(p,source/p.name);hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    result=calculate();result['source_sha256']=hashes
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['reduced_predictions']),flush=True)


if __name__=='__main__':main()
