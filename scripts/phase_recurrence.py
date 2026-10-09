#!/usr/bin/env python3
"""Independent unit-density, zero-velocity D2Q9 phase-transport diagnostic.

This isolates recoloring and impermeable reflection from fluid dynamics.
It tests a discrete stationary-profile condition, not an interface PDE fit.
"""
import argparse
import json
from pathlib import Path

import numpy as np

VELOCITIES=np.array([(0,0),(1,0),(0,1),(-1,0),(0,-1),
                     (1,1),(-1,1),(-1,-1),(1,-1)])
WEIGHTS=np.array([4/9]+[1/9]*4+[1/36]*4)
OPPOSITE=np.array([0,3,4,1,2,7,8,5,6])


def coefficients():
    a,b=1/9,1/36
    bulk_a,bulk_b=a+2*b,a+np.sqrt(2)*b
    wall_a=a+b
    factor=(wall_a*bulk_b/bulk_a-b/np.sqrt(2))/a
    return dict(coordinate_weight=a,diagonal_weight=b,bulk_diffusion=bulk_a,
                bulk_segregation=bulk_b,wall_diffusion=wall_a,
                legacy_wall_segregation=a+b/np.sqrt(2),tangential_factor=float(factor),
                corrected_wall_segregation=a*factor+b/np.sqrt(2),
                wall_residual_ratio=wall_a/bulk_a)


def bulk_profile(length,beta,steps=8000):
    phi=np.zeros(length);phi[length//4:3*length//4]=1
    c=coefficients();a,b=c['bulk_diffusion'],c['bulk_segregation']
    for _ in range(steps):
        color=2*phi-1
        gradient=(np.roll(color,-1)-np.roll(color,1))/2
        n=np.where(abs(gradient)>1e-9,np.sign(gradient),0.)
        flux=n*phi*(1-phi)
        phi=((1-2*a)*phi+a*(np.roll(phi,1)+np.roll(phi,-1))+
             beta*b*(np.roll(flux,1)-np.roll(flux,-1)))
    return phi


def step(phi,beta,tangential=False):
    height,length=phi.shape
    color=2*phi-1
    gradient=np.zeros((2,height,length))
    boundary=np.zeros_like(phi,dtype=bool);boundary[[0,-1]]=True
    for (dx,dy),weight in zip(VELOCITIES[1:],WEIGHTS[1:]):
        neighbor=np.roll(color,(-dy,-dx),axis=(0,1))
        valid=np.ones_like(phi,dtype=bool)
        if dy>0:valid[-1]=False
        if dy<0:valid[0]=False
        difference=np.where(boundary,neighbor-color,neighbor)
        gradient[0]+=3*weight*dx*difference*valid
        gradient[1]+=3*weight*dy*difference*valid
    # Available-stencil inverse for the renormalized wall gradient.
    gradient[0,[0,-1]]*=6/5
    gradient[1,[0,-1]]*=2
    if tangential:gradient[1,[0,-1]]=0
    norm=np.linalg.norm(gradient,axis=0)
    normal=np.divide(gradient,norm[None],out=np.zeros_like(gradient),where=norm[None]>1e-9)
    normal[0,[0,-1]]=np.where(norm[[0,-1]]>1e-9,
                            np.where(gradient[0,[0,-1]]>=0,1.,-1.),0)
    normal[1,[0,-1]]=0
    outgoing=np.zeros((9,height,length))
    for q,((dx,dy),weight) in enumerate(zip(VELOCITIES,WEIGHTS)):
        projection=(dx*normal[0]+dy*normal[1])/np.hypot(dx,dy) if q else np.zeros_like(phi)
        if tangential and q and dy==0:projection[[0,-1]]*=coefficients()['tangential_factor']
        outgoing[q]=weight*(phi+beta*phi*(1-phi)*projection)
    incoming=[]
    for q,(dx,dy) in enumerate(VELOCITIES):
        shifted=np.roll(outgoing[q],(dy,dx),axis=(0,1))
        if dy>0:shifted[0]=outgoing[OPPOSITE[q],0]
        if dy<0:shifted[-1]=outgoing[OPPOSITE[q],-1]
        incoming.append(shifted)
    return np.sum(incoming,axis=0)


def equilibrium_seed(target,length,height,beta,axis='x',tangential=True,time=1):
    """Write independently reconstructed outgoing populations, not a VTI fit."""
    phi=bulk_profile(length,beta)
    color=2*phi-1
    gradient=(np.roll(color,-1)-np.roll(color,1))/2
    normal=np.where(abs(gradient)>1e-9,np.sign(gradient),0.)
    populations=np.zeros((2,9,height+2,length))
    for q,((dx,dy),weight) in enumerate(zip(VELOCITIES,WEIGHTS)):
        projection=dx*normal/np.hypot(dx,dy) if q else np.zeros(length)
        f0=weight*(phi+beta*phi*(1-phi)*projection)
        populations[0,q,1:-1]=f0
        populations[1,q,1:-1]=weight-f0
        if tangential and q and dy==0:
            wall=weight*(phi+beta*phi*(1-phi)*projection*coefficients()['tangential_factor'])
            populations[0,q,[1,-2]]=wall
            populations[1,q,[1,-2]]=weight-wall
    if axis=='y':
        # Transpose the spatial field and the velocity directions together.
        swap=[int(np.where(np.all(VELOCITIES==v[::-1],axis=1))[0][0]) for v in VELOCITIES]
        populations=populations[:,swap].transpose(0,1,3,2)
    elif axis!='x':raise ValueError('The wall tangent must be x or y')
    target=Path(target);target.parent.mkdir(parents=True,exist_ok=True)
    Path(str(target)+'.rst').write_bytes(populations.astype('=f8').tobytes())
    ny,nx=populations.shape[-2:]
    Path(str(target)+'.rst.meta').write_text(
        f'DRLAB_OUTGOING_V1\n{time} 2 2 9\n{ny} {nx} \nNEXT_STREAM {time+1}\n'
        f'INITIAL_COLLISION_INCLUDED {time+1}\nCONVERGENCE_HISTORY_RESET_ON_LOAD\n')
    return phi


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--beta',type=float,default=.7)
    parser.add_argument('--height',type=int,default=16)
    parser.add_argument('--steps',type=int,default=8000)
    parser.add_argument('--baseline-receipt',help='Compare the phase-only diagnostic to a recorded zero-sigma solver case')
    args=parser.parse_args()
    length=4*args.height
    bulk=bulk_profile(length,args.beta,args.steps)
    lifted=np.broadcast_to(bulk,(args.height,length)).copy()
    residuals={}
    fields={}
    for correction in (False,True):
        mode='tangential' if correction else 'legacy'
        residuals[mode]=float(np.max(abs(step(lifted,args.beta,correction)-lifted)))
        phi=np.zeros_like(lifted);phi[:,length//4:3*length//4]=1
        for _ in range(args.steps):phi=step(phi,args.beta,correction)
        fields[mode]=phi
    result=dict(protocol='Independent zero-velocity, unit-total-density phase recurrence; renormalized neutral-wall projection; no hydrodynamic solver or force.',
                coefficients=coefficients(),beta=args.beta,height=args.height,steps=args.steps,
                lifted_bulk_fixedpoint_residuals=residuals,
                wall_color_distortions={k:float(np.max(abs(2*(v[0]-v[args.height//2])))) for k,v in fields.items()},
                final_phase_profiles={k:v.tolist() for k,v in fields.items()})
    if args.baseline_receipt:
        from benchmarks import geometry_fields
        receipt=json.loads(Path(args.baseline_receipt).read_text())
        assert receipt['parameters']['sigma']==0
        arrays,rho,chi,velocity=geometry_fields(receipt)
        observed=(chi[0,1:-1]+1)/2
        result['max_phase_difference_from_full_solver']=float(np.max(abs(observed-fields['legacy'])))
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('final_phase_profiles',)}))


if __name__=='__main__':main()
