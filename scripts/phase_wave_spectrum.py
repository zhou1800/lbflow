#!/usr/bin/env python3
"""Independent linear spectrum of neutral projection phase-shape transport."""
import argparse
import json
from pathlib import Path

import numpy as np

from front_mobility import match_beta, WIDTH_REFERENCE
from phase_recurrence import bulk_profile, coefficients


def matrix(profile,beta,wavenumber,cutoff=1e-9):
    """Exact zero-flow stencil symbol; variation of the unit normal is included."""
    n=len(profile);eye=np.eye(n);left=np.roll(eye,1,axis=0);right=np.roll(eye,-1,axis=0)
    g=(np.roll(profile,-1)-np.roll(profile,1))/2
    # Production normal uses the color gradient of chi=2*phi-1.
    active=2*abs(g)>cutoff;normal=np.where(active,np.sign(g),0.)
    a,b=1/9,1/36;c=np.cos(wavenumber);s=np.sin(wavenumber)
    diffusion=(4/9+2*a*c)*eye+(a+2*b*c)*(left+right)
    coefficient=(1-2*profile)*normal
    normal_part=beta*(a+2*b*c)*(np.roll(coefficient[:,None]*eye,1,axis=0)-np.roll(coefficient[:,None]*eye,-1,axis=0))
    smooth=(2/3)*eye+(left+right)/6
    fraction=np.divide(profile*(1-profile),abs(g),out=np.zeros_like(g),where=active)
    tangent=fraction[:,None]*smooth
    tangent_part=2*beta*s*s*(a*tangent+b*(np.roll(tangent,1,axis=0)+np.roll(tangent,-1,axis=0)))
    return diffusion+normal_part+tangent_part


def expansion_controls():
    """Diagnostic no-cutoff equilibrium checks of the formal 1/24 limit."""
    from stokes_phase_reference import planar_bridge, projection_rest
    rows=[]
    for height in (32,64,128):
        h=16/height;beta=match_beta(WIDTH_REFERENCE/h)
        profile=planar_bridge(4*height,beta)
        stationary=float(abs(projection_rest(np.broadcast_to(profile,(4,len(profile))),beta,cutoff=0.)-profile[None]).max())
        values,vectors=np.linalg.eig(matrix(profile,beta,8*np.pi/height,cutoff=0.))
        derivative=(np.roll(profile,-1)-np.roll(profile,1))/2
        overlap=abs(vectors.conj().T@derivative)/np.linalg.norm(vectors,axis=0)/np.linalg.norm(derivative)
        value=values[int(np.argmax(overlap))]
        if abs(value.imag)>1e-12 or not 0<value.real<1 or stationary>3e-14:
            raise AssertionError('Diagnostic no-cutoff equilibrium or phase mode failed')
        rate=float(-np.log(value.real)/h**2);coefficient=rate/(h*h*(8*np.pi/16)**4)
        rows.append(dict(height=height,mode=8,cutoff=0.,profile='Analytic zero-flux face recurrence',
            max_stationary_residual=stationary,physical_relaxation_rate=rate,
            coefficient=coefficient,relative_coefficient_difference=coefficient*24-1))
    errors=[abs(r['relative_coefficient_difference']) for r in rows]
    orders=[float(np.log2(a/b)) for a,b in zip(errors,errors[1:])]
    if min(orders)<1.8 or errors[-1]>.006:raise AssertionError('Formal coefficient trend is not reproduced')
    return dict(controls=rows,coefficient_error_orders=orders,
        scope='Exact matrix on well-separated periodic planar equilibria with diagnostic cutoff zero; mathematical expansion control, not a production LBM trajectory or a universal-interface validation.')


def spectrum(height,mode=1,profile_steps=8000,width=WIDTH_REFERENCE):
    from stokes_phase_reference import projection_rest
    if height<16 or height%2 or not 1<=mode<height or profile_steps<1:
        raise ValueError('Use an even height >=16, a resolved transverse mode and positive profile duration')
    if not np.isfinite(width) or width<=0:raise ValueError('Use a finite positive physical width')
    dx=16/height;length=4*height;beta=match_beta(width/dx)
    ratio=coefficients()['bulk_segregation']/coefficients()['bulk_diffusion']
    profile=bulk_profile(length,beta/ratio,steps=profile_steps);wave=mode*np.pi/height
    stationary=float(abs(projection_rest(np.broadcast_to(profile,(4,length)),beta)-profile[None]).max())
    if stationary>3e-14:
        raise ValueError('Planar profile is not stationary; increase --profile-steps before interpreting its spectrum')
    values,vectors=np.linalg.eig(matrix(profile,beta,wave))
    derivative=(np.roll(profile,-1)-np.roll(profile,1))/2
    overlap=abs(np.conj(vectors).T@derivative)/np.linalg.norm(vectors,axis=0)/np.linalg.norm(derivative)
    chosen=int(np.argmax(overlap));value=values[chosen]
    if abs(value.imag)>1e-12 or value.real<=0 or value.real>=1:raise AssertionError('Expected a stable real translational phase mode')
    rate=float(-np.log(value.real)/dx**2)
    return dict(height=height,mode=mode,beta=beta,wavenumber_lattice=wave,
        eigenvalue=float(value.real),physical_relaxation_rate=rate,
        rate_divided_by_dx_squared=rate/dx**2,translation_mode_overlap=float(overlap[chosen]),
        profile=profile.tolist(),eigenvector=np.real(vectors[:,chosen]).tolist(),
        profile_steps=profile_steps,max_stationary_residual=stationary,
        spectral_radius=float(abs(values).max()))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    parser.add_argument('--heights',type=int,nargs='+',default=[16,32,48]);parser.add_argument('--modes',type=int,nargs='+',default=[1])
    parser.add_argument('--profile-steps',type=int,default=8000)
    parser.add_argument('--width',type=float,default=WIDTH_REFERENCE)
    args=parser.parse_args();rows=[spectrum(h,m,args.profile_steps,args.width) for h in args.heights for m in args.modes]
    result=dict(measurements=rows,protocol='Analytic transverse Fourier symbol of the exact projection phase stencil at an independently relaxed planar profile. The unit normal is differentiated; the production color-gradient cutoff is retained; no relaxation rate is fitted.',
        scope='Zero-flow, zero-sigma, neutral reflected scalar operator. No moving-contact-line constitutive law or full hydrodynamic model implied.')
    if args.width!=WIDTH_REFERENCE:result['physical_width']=args.width
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ('height','mode','physical_relaxation_rate','translation_mode_overlap')} for r in rows],indent=2))


if __name__=='__main__':main()
