#!/usr/bin/env python3
"""Independent wall fixed points, exact period-two seeds and replay checks."""
import argparse
import json
from pathlib import Path
import uuid

import numpy as np

from benchmarks import Study,compare,geometry_fields,spec
from phase_recurrence import equilibrium_seed,VELOCITIES,WEIGHTS
from wall_consistency import cases


def checkerboard_seed(target,model,amplitude=1e-6,time=1):
    height,length=4,8
    velocity=amplitude*(-1.)**np.arange(length)
    populations=np.zeros((1 if model=='scsp' else 2,9,height,length))
    for q,((dx,dy),weight) in enumerate(zip(VELOCITIES,WEIGHTS)):
        eu=dx*velocity
        populations[0,q]=weight*(1+3*eu+4.5*eu**2-1.5*velocity**2)
    Path(str(target)+'.rst').write_bytes(populations.astype('=f8').tobytes())
    Path(str(target)+'.rst.meta').write_text(
        f'DRLAB_OUTGOING_V1\n{time} {len(populations)} 2 9\n{height} {length} \n'
        f'NEXT_STREAM {time+1}\nINITIAL_COLLISION_INCLUDED {time+1}\nCONVERGENCE_HISTORY_RESET_ON_LOAD\n')
    return velocity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',required=True);p.add_argument('--root',required=True)
    p.add_argument('--mpiexec',default='mpiexec');p.add_argument('--numproc-flag',default='-n')
    args=p.parse_args()
    study=Study(args.solver,Path(args.root)/('wall-checks-'+uuid.uuid4().hex[:10]),args.mpiexec,args.numproc_flag)
    fixed=[]
    for beta in (.35,.7,.9):
        for length in (64,65):
            for form in ('normal','gradient'):
                cfg=next(cases([16],[beta],['tangential'],[.005],['csf']))
                cfg.update(name=f'fixed-beta{beta:g}-L{length}-{form}',size=(length,18,1),
                           steps=501,write=250,restart=250,curvature_stencil=form)
                seed=study.root/f'seed-beta{beta:g}-L{length}'
                phi=equilibrium_seed(seed,length,16,beta)
                r=study.run(cfg,start=1,checkpoint=seed)
                _,rho,chi,v=geometry_fields(r)
                speed=float(np.linalg.norm(v,axis=-1).max())
                error=float(abs(chi[0,1:-1]-(2*phi-1)).max())
                if speed>2e-14 or error>2e-11:raise AssertionError('Manufactured wall fixed point was not preserved')
                fixed.append(dict(name=r['name'],max_speed=speed,color_error=error))
                if beta==.7 and length==64 and form=='gradient':
                    reference=r;reference_cfg=cfg
    parallel=study.run(dict(reference_cfg,name='fixed-parallel'),ranks=2,start=1,
                       checkpoint=study.root/'seed-beta0.7-L64')
    comparisons=[compare(reference,parallel)]
    checkpoint=Path(reference['output_files'][-1]).parent/'field.0000000250'
    for ranks in (1,2):
        resumed=study.run(dict(reference_cfg,name=f'fixed-resumed-{ranks}'),ranks,250,checkpoint)
        comparisons.append(compare(reference,resumed))
    rotated=next(cases([16],[.7],['tangential'],[.005],['csf'],axes=['y']))
    rotated.update(name='fixed-rotated-SRT',eq='SRT',steps=500,write=250,curvature_stencil='gradient')
    seed=study.root/'seed-rotated';phi=equilibrium_seed(seed,64,16,.7,axis='y')
    r=study.run(rotated,start=1,checkpoint=seed);_,rho,chi,v=geometry_fields(r)
    speed=float(np.linalg.norm(v,axis=-1).max());error=float(abs(chi[0,:,1:-1]-(2*phi[:,None]-1)).max())
    if speed>2e-14 or error>2e-11:raise AssertionError('Rotated fixed point was not preserved')
    fixed.append(dict(name=r['name'],max_speed=speed,color_error=error))
    convergence=[]
    for model in ('scsp','color'):
        for collision in ('SRT','MRT'):
            seed=study.root/f'checkerboard-{model}';u=checkerboard_seed(seed,model)
            cfg=spec(f'checkerboard-{model}-{collision}',nx=8,ny=4,nz=1,dimensions=2,
                     model=model,eq=collision,tau=(1.,1.) if model=='color' else (1.,),
                     steps=16,write=2,report_dir=1,stop_error=1e-6)
            r=study.run(cfg,start=1,checkpoint=seed)
            _,rho,chi,v=geometry_fields(r)
            mode_error=float(abs(v[0,:,:,0]+u).max())
            if r['timing']['completed_steps']!=15 or mode_error>2e-14:
                raise AssertionError('Period-two mode falsely converged or lost its exact propagation')
            log=(Path(r['output_files'][-1]).parent/'run.log').read_text()
            if 'adjacentVelocityResidual' not in log:raise AssertionError('Adjacent-step convergence was not exercised')
            settled=study.run(dict(cfg,name=f'settled-{model}-{collision}'))
            if settled['timing']['completed_steps']>=16:raise AssertionError('Stationary flow no longer converges')
            convergence.append(dict(name=r['name'],mode_error=mode_error,
                                    oscillatory_steps=r['timing']['completed_steps'],settled_steps=settled['timing']['completed_steps']))
    guards=[]
    base=next(cases([16],[.7],['tangential'],[.005],['csf']))
    base.update(steps=2,write=1)
    for name,kw in [('bad-recoloring',dict(overrides={'wallRecoloring':'"unknown"'})),
                    ('bad-curvature',dict(overrides={'curvatureStencil':'"unknown"'})),
                    ('gradient-with-perturbation',dict(capillary_model='perturbation',curvature_stencil='gradient')),
                    ('unsupported-gradient',dict(wall_gradient='legacy')),
                    ('unsupported-3d',dict(dimensions=3,size=(64,18,4))),
                    ('wall-corner',dict(label=lambda x,y,z:255 if y in (0,17) or (x==8 and y==8) else 0))]:
        r=study.run(dict(base,name=name,expect_failure=True,**kw))
        guards.append(dict(name=name,exit=r['exit']))
    summary=dict(fixed_points=fixed,comparisons=comparisons,convergence=convergence,guards=guards,
                 solver_sha256=study.solver_sha,
                 scope='Flat neutral D2Q9 fixed points and short-time propagation/replay; no dynamic or irregular-wall accuracy claim.')
    (study.root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
