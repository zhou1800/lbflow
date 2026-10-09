#!/usr/bin/env python3
"""Separate alternating momentum, capillary forcing and maintained circulation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from benchmarks import Study, fields, geometry_fields, sha
from capillary_response import cases, measure as pressure_measure
from staggered_diagnostics import measure as adjacent_measure


def physical_projection(receipt, time):
    _, rho, chi, velocity = geometry_fields(receipt, time)
    z, y, x = np.indices(rho.shape); fluid = rho > 0
    return np.array([np.mean((rho*velocity[..., d]*(-1.)**coord)[fluid])
                     for d, coord in enumerate((x, y))])


def force_audit(study, receipt, driver, time):
    original = Path(receipt['output_files'][-1]).parent
    folder = study.root/'audit'/receipt['name']/str(time); folder.mkdir(parents=True)
    for filename in ('geometry.raw', f'field.{time:010d}.rst', f'field.{time:010d}.rst.meta'):
        shutil.copyfile(original/filename, folder/filename)
    text = (original/'input.prm').read_text().replace(str(original), str(folder))
    text = text.replace(f'startTime {receipt["start"]};', f'startTime {time};')
    inp = folder/'input.prm'; inp.write_text(text)
    process = subprocess.run([str(driver), str(inp), '--capillary-state'], cwd=folder,
                             env=study.env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, timeout=30)
    (folder/'audit.log').write_text(process.stdout)
    if process.returncode: raise RuntimeError('Production capillary-state audit failed')
    data = json.loads(next(line[len('CAPILLARY_STATE '):] for line in process.stdout.splitlines()
                           if line.startswith('CAPILLARY_STATE ')))
    shape = data['size']; rho = np.array(data['density']).reshape(shape); fluid = rho > 0
    y, x = np.indices(shape)
    force = np.array(data['force']).reshape((2,)+tuple(shape))[::-1]
    projection = np.array([np.mean((force[d]*(-1.)**coord)[fluid]) for d, coord in enumerate((x, y))])
    populations = np.fromfile(folder/f'field.{time:010d}.rst', dtype='=f8').reshape((2, 9)+tuple(shape))
    from phase_recurrence import VELOCITIES
    momentum = np.einsum('qyx,qd->dyx', populations.sum(axis=0), VELOCITIES)
    outgoing = np.array([np.mean((momentum[d]*(-1.)**coord)[fluid]) for d, coord in enumerate((x, y))])
    return dict(step=time, force_projection=projection.tolist(), outgoing_projection=outgoing.tolist(),
                input_sha256=sha(inp), checkpoint_sha256=sha(folder/f'field.{time:010d}.rst'))


def adjacent(study, seed, name, sigma=None, steps=8):
    cfg = dict(seed['parameters']); cut = cfg['steps']; folder = Path(seed['output_files'][-1]).parent
    cfg.update(name=name, steps=cut+steps, write=1, restart=1,
               geometry_file=str(folder/'geometry.raw'), report_dir=0)
    if sigma is not None: cfg['sigma'] = sigma
    return study.run(cfg, start=cut, checkpoint=folder/f'field.{cut:010d}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True); parser.add_argument('--audit-driver', required=True)
    parser.add_argument('--root', required=True); parser.add_argument('--seed', required=True)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args(); study = Study(args.solver, args.root)
    source_driver = Path(args.audit_driver).resolve()
    driver = study.root/'.bin'/sha(source_driver)/source_driver.name
    driver.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(source_driver, driver); driver.chmod(source_driver.stat().st_mode)
    seed = json.loads(Path(args.seed).read_text()); short = adjacent(study, seed, 'held-force-adjacent')
    cut = seed['parameters']['steps']; points=[]
    for time in range(cut+1, cut+9):
        point = force_audit(study, short, driver, time)
        point['physical_projection'] = physical_projection(short, time).tolist(); points.append(point)
    residuals=[]; collision_residuals=[]
    for left, right in zip(points, points[1:]):
        p, next_p = np.array(left['physical_projection']), np.array(right['physical_projection'])
        g, next_g = np.array(left['force_projection']), np.array(right['force_projection'])
        residuals.append((next_p+p-.5*(next_g-g)).tolist())
        collision_residuals.append((np.array(left['outgoing_projection'])-p-.5*g).tolist())
    if np.max(abs(np.array(residuals))) > 3e-12: raise AssertionError('Alternating momentum/force relation failed')
    # A changed-physics continuation diagnoses the free total-fluid modes.
    cfg = dict(seed['parameters']); folder = Path(seed['output_files'][-1]).parent
    cfg.update(name='force-removed-relaxation', sigma=0., steps=cut+2048, write=256,
               restart=cut+2048, geometry_file=str(folder/'geometry.raw'), report_dir=0)
    removed = study.run(cfg, start=cut, checkpoint=folder/f'field.{cut:010d}')
    removal_history=[]
    for filename in removed['output_files']:
        time=int(Path(filename).name.split('.')[1]); _,rho,chi,v=geometry_fields(removed,time)
        removal_history.append(dict(step=time, max_speed=float(np.linalg.norm(v,axis=-1).max()),
                                    physical_projection=physical_projection(removed,time).tolist()))
    removed_adjacent=adjacent(study,removed,'force-removed-adjacent',steps=2)
    removal_diagnostics=adjacent_measure(removed,removed_adjacent)
    configs=[]
    for length,collision,shift in [(64,'MRT',0),(65,'MRT',0),(64,'SRT',0),(65,'SRT',0),(64,'MRT',1)]:
        cfg=next(cases([16],[60],['ghost'],capillary_models=['csf']))
        cfg.update(name=f'parity-L{length}-{collision}-shift{shift}',size=(length,18,1),eq=collision,
                   phase_boundary='flux',steps=32000,write=8000,restart=32000,
                   label=lambda x,y,z,s=shift,l=length:255 if y in (0,17) else (0 if 16<=((x-s)%l)<48 else 100))
        configs.append(cfg)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool: runs=list(pool.map(study.run,configs))
    measurements=[]
    for run in runs:
        replay=adjacent(study,run,run['parameters']['name']+'-adjacent',steps=2)
        measurements.append(dict(pressure=pressure_measure(run),adjacent=adjacent_measure(run,replay)))
    untranslated=next(r for r in runs if r['parameters']['size'][0]==64 and r['parameters']['eq']=='MRT' and 'shift0' in r['name'])
    translated=next(r for r in runs if 'shift1' in r['name'])
    left,right=fields(untranslated),fields(translated)
    translation_errors={key:float(abs(np.roll(left[key],1,axis=2)-right[key]).max()) for key in left}
    if max(translation_errors.values()) > 2e-12:raise AssertionError('Translated periodic geometry changes the solution')
    result=dict(protocol='Even-grid force/alternating-momentum identity; explicit force-removal diagnostic; even/odd domain and exact translation controls.',
                seed_receipt_sha256=sha(args.seed),solver_sha256=study.solver_sha,audit_driver_sha256=sha(driver),
                force_identity_points=points,identity_residuals=residuals,collision_residuals=collision_residuals,
                force_removed_history=removal_history,force_removed_adjacent=removal_diagnostics,
                parity=measurements,translation_errors=translation_errors,
                scope='Diagnostic changed physics and small odd-domain geometry difference; no damping/filtering or general wetting accuracy claim.')
    (study.root/'parity-forcing.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(max_identity_residual=float(np.max(abs(np.array(residuals)))),
                         force_removed_pair_speed=removal_diagnostics['adjacent_pair_average_speed'],translation_errors=translation_errors)))


if __name__=='__main__':main()
