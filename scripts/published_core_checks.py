#!/usr/bin/env python3
"""Independent recoloring, momentum, wall and replay checks for core controls."""
import argparse
import json
from pathlib import Path
import subprocess
import uuid

import numpy as np

from benchmarks import Study, compare, geometry_fields, sha, spec
from capillary_response import cases
from phase_recurrence import VELOCITIES, WEIGHTS, bulk_profile, coefficients
from wall_equilibrium_checks import checkerboard_seed
from wall_flux_checks import bulk_step


def projection_seed(path, length=65, height=16):
    beta = .7*coefficients()['bulk_segregation']/coefficients()['bulk_diffusion']
    phi = bulk_profile(length, .7)
    gradient = (np.roll(2*phi-1, -1)-np.roll(2*phi-1, 1))/2
    normal = np.where(abs(gradient)>1e-9, np.sign(gradient), 0.)
    f = np.zeros((2, 9, height+2, length))
    for q, ((dx, dy), w) in enumerate(zip(VELOCITIES, WEIGHTS)):
        phase = w*(phi+beta*phi*(1-phi)*dx*normal)
        f[0, q, 1:-1] = phase
        f[1, q, 1:-1] = w-phase
    Path(str(path)+'.rst').write_bytes(f.astype('=f8').tobytes())
    Path(str(path)+'.rst.meta').write_text(
        f'DRLAB_OUTGOING_V1\n1 2 2 9\n{height+2} {length} \n'
        'NEXT_STREAM 2\nINITIAL_COLLISION_INCLUDED 2\nCONVERGENCE_HISTORY_RESET_ON_LOAD\n')
    return phi, beta


def recoloring_check(receipt):
    cfg = receipt['parameters']; nx, ny, _ = cfg['size']
    folder = Path(receipt['output_files'][-1]).parent
    f = np.fromfile(folder/'field.0000000001.rst', dtype='=f8').reshape(2, 9, ny, nx)
    _, rho, chi, _ = geometry_fields(receipt)
    rho, chi = rho[0], chi[0]
    gradient = np.zeros((2, ny, nx))
    for (dx, dy), w in zip(VELOCITIES, WEIGHTS):
        gradient[0] += 3*w*dx*np.roll(chi, (-dy, -dx), axis=(0, 1))
        gradient[1] += 3*w*dy*np.roll(chi, (-dy, -dx), axis=(0, 1))
    norm = np.linalg.norm(gradient, axis=0)
    normal = np.divide(gradient, norm[None], out=np.zeros_like(gradient), where=norm[None]>1e-9)
    projection = np.einsum('qd,dyx->qyx', VELOCITIES, normal)
    if cfg['recoloring_model']=='cosine':
        lengths = np.linalg.norm(VELOCITIES, axis=1)
        projection = np.divide(projection, lengths[:, None, None],
                               out=np.zeros_like(projection), where=lengths[:, None, None]>0)
    phi = (1+chi)/2
    expected = phi[None]*f.sum(axis=0)+cfg['beta']*rho[None]*phi[None]*(1-phi[None])*WEIGHTS[:, None, None]*projection
    error = float(abs(f[0]-expected).max())
    if error>3e-14: raise AssertionError('Actual outgoing recoloring disagrees with the independent moment formula')
    return dict(name=receipt['name'], max_population_error=error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True); parser.add_argument('--audit-driver', required=True)
    parser.add_argument('--root', required=True); parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--numproc-flag', default='-n')
    args = parser.parse_args()
    study = Study(args.solver, Path(args.root)/('published-core-'+uuid.uuid4().hex[:10]), args.mpiexec, args.numproc_flag)
    recoloring, transients, fixed, forcing, parity, guards = [], [], [], [], [], []
    for model in ('cosine', 'projection'):
        cfg = spec('population-'+model, nx=32, ny=16, nz=1, dimensions=2,
                   steps=1, write=1, sigma=0., capillary_model='csf',
                   recoloring_model=model, mrt_spectrum='liu2015',
                   label=lambda x,y,z: 0 if (x+2*y)%32<16 else 100,
                   indices=['0 "VD" "N/A" "N/A" "0"', '100 "VD" "N/A" "N/A" "1"'])
        recoloring.append(recoloring_check(study.run(cfg)))
    for spectrum in ('legacy', 'liu2015'):
        cfg = next(cases([16], [90], ['ghost'], capillary_models=['csf']))
        beta = .7*coefficients()['bulk_segregation']/coefficients()['bulk_diffusion']
        cfg.update(name='transient-'+spectrum, size=(65, 18, 1), beta=beta,
                   steps=64, write=1, phase_boundary='flux', recoloring_model='projection', mrt_spectrum=spectrum)
        r = study.run(cfg); phi = np.zeros(65); phi[16:48]=1; error=speed=0.
        for t in range(65):
            if t: phi=bulk_step(phi, .7)  # Analytic beta match makes the full recurrence identical.
            _, rho, chi, v = geometry_fields(r, t)
            error=max(error, float(abs((1+chi[0, 1:-1])/2-phi).max()))
            speed=max(speed, float(np.linalg.norm(v, axis=-1).max()))
        if error>3e-12 or speed>2e-14: raise AssertionError('Published projection lost matched neutral transient')
        transients.append(dict(name=r['name'], max_color_error=error, max_speed=speed))
        seed = study.root/('projection-seed-'+spectrum); phi, beta=projection_seed(seed)
        cfg.update(name='bounce-fixed-'+spectrum, steps=129, write=64, phase_boundary='bounce')
        r=study.run(cfg, start=1, checkpoint=seed); _,rho,chi,v=geometry_fields(r)
        error=float(abs(chi[0,1:-1]-(2*phi-1)).max());speed=float(np.linalg.norm(v,axis=-1).max())
        if error>2e-11 or speed>2e-14: raise AssertionError('Published projection already has a neutral bounce-back fixed point')
        fixed.append(dict(name=r['name'], color_error=error, max_speed=speed))
        cfg=spec('force-'+spectrum, nx=8, ny=4, nz=1, dimensions=2,
                 tau=(.6,1.4), steps=100, write=50, capillary_model='csf',
                 recoloring_model='projection', mrt_spectrum=spectrum,
                 label=lambda x,y,z:75, indices=['75 "GV" "0.000001" "1" "1"'])
        r=study.run(cfg);v=geometry_fields(r)[3];error=float(abs(v[...,0]-1e-4).max())
        if error>1e-12 or abs(v[...,1:]).max()>1e-12: raise AssertionError('Momentum rates and half-force convention do not integrate the prescribed force')
        forcing.append(dict(name=r['name'], error=error))
        if spectrum=='liu2015':
            inp=Path(r['output_files'][-1]).parent/'input.prm'
            proc=subprocess.run([args.audit_driver,str(inp)],env=study.env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
            if proc.returncode: raise AssertionError('Published MRT audit failed')
            data=json.loads(next(line[len('MOMENT_AUDIT '):] for line in proc.stdout.splitlines() if line.startswith('MOMENT_AUDIT ')))
            rates=np.array(data['rates']);expected=np.array([1,1.63,1.14,1,1.92,1,1.92,1/1.2,1/1.2])
            if abs(rates-expected).max()>2e-14: raise AssertionError('Published MRT rates differ from specified spectrum')
            (study.root/'operators.json').write_text(json.dumps(data,indent=2)+'\n')
            audit=dict(max_rate_error=float(abs(rates-expected).max()), driver_sha256=sha(args.audit_driver), input_sha256=sha(inp))
        seed=study.root/('checkerboard-'+spectrum);u=checkerboard_seed(seed,'color')
        cfg.update(name='checkerboard-'+spectrum, tau=(1.,1.), report_dir=1, stop_error=1e-6,
                   steps=16, write=2, indices=['0 "VD" "N/A" "N/A" "0"'], label=lambda x,y,z:0)
        r=study.run(cfg,start=1,checkpoint=seed);v=geometry_fields(r)[3]
        error=float(abs(v[0,:,:,0]+u).max())
        if error>2e-14 or r['timing']['completed_steps']!=15: raise AssertionError('Published MRT falsely removes or converges the exact alternating equilibrium mode')
        parity.append(dict(name=r['name'], error=error, completed_steps=r['timing']['completed_steps']))
    cfg=next(cases([16],[60],['ghost'],capillary_models=['csf']))
    cfg.update(name='nonneutral-core', recoloring_model='projection', mrt_spectrum='liu2015',
               beta=beta, phase_boundary='flux', tau=(.6,1.4), viscosity_interpolation='harmonic',
               steps=128, write=64, restart=64)
    ref=study.run(cfg);par=study.run(dict(cfg,name='nonneutral-core-parallel'),2)
    comparisons=[compare(ref,par)];populations=[]
    checkpoint=Path(ref['output_files'][-1]).parent/'field.0000000064'
    for rank in (1,2):
        r=study.run(dict(cfg,name='nonneutral-core-resumed-'+str(rank)),rank,64,checkpoint)
        comparisons.append(compare(ref,r))
        populations.append(r)
    populations.append(par)
    final=Path(ref['output_files'][-1]).parent/'field.0000000128.rst'
    for r in populations:
        other=Path(r['output_files'][-1]).parent/'field.0000000128.rst'
        if final.read_bytes()!=other.read_bytes():raise AssertionError('Core MPI/restart outgoing populations differ')
    cfg.update(steps=2,write=1)
    invalid=[('bad-recoloring',dict(recoloring_model='unknown')),
             ('bad-spectrum',dict(mrt_spectrum='unknown')),
             ('srt-spectrum',dict(eq='SRT')),
             ('perturbation-spectrum',dict(capillary_model='perturbation')),
             ('wrong-weights',dict(alpha=.4)),
             ('wrong-density',dict(density=(1.,2.))),
             ('three-dimensions',dict(dimensions=3,size=(64,18,4))),
             ('duplicate-multiplier',dict(wall_recoloring='tangential')),
             ('single-component',dict(model='scsp',tau=(1.,),density=(1.,)))]
    for name,change in invalid:
        r=study.run(dict(cfg,name=name,expect_failure=True,**change));guards.append(dict(name=name,exit=r['exit']))
    summary=dict(recoloring=recoloring,matched_transients=transients,published_fixed_points=fixed,
                 forcing=forcing,parity=parity,mrt_audit=audit,comparisons=comparisons,
                 bitwise_population_comparisons=len(populations),guards=guards,solver_sha256=study.solver_sha)
    (study.root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
