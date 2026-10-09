#!/usr/bin/env python3
"""Audit recorded capillary forces with an independently verified Stokes solve.

Two adjacent physical states are retained, without averaging inside the
production solver. The same-phase-field force comparison is a diagnostic,
not an independent validation of interface motion.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil

import numpy as np

from benchmarks import Study, geometry_fields, sha
from parity_forcing import adjacent, force_audit
from stokes_phase_reference import StokesSlit, capillary_faces, match_logistic_width
from front_mobility import WIDTH_REFERENCE


def remeasure(receipt, root):
    cfg=receipt['parameters'];height=cfg['height'];h=cfg['dx']
    solver=StokesSlit(height,slip=cfg['physical_slip'])
    points=[]
    for step in (receipt['start']+1,receipt['start']+2):
        folder=Path(root)/'audit'/receipt['name']/str(step)
        data=json.loads(next(line[len('CAPILLARY_STATE '):]
            for line in (folder/'audit.log').read_text().splitlines()
            if line.startswith('CAPILLARY_STATE ')))
        shape=tuple(data['size'])
        density=np.array(data['density']).reshape(shape)
        color=np.array(data['color']).reshape(shape)
        _,rho,chi,velocity=geometry_fields(receipt,step);fluid=density>0
        state_error=max(float(abs(density-rho[0]).max()),float(abs(color[fluid]-chi[0][fluid]).max()))
        if state_error>3e-12:raise AssertionError('Audit scalar state differs from its physical output')
        # Audit stores coordinate components in the historical reversed order.
        force=np.array(data['force']).reshape((2,)+shape)[::-1,1:-1]/h**3
        fx=(force[0]+np.roll(force[0],1,axis=1))/2+cfg['physical_force']
        fy=np.zeros((height+1,4*height));fy[1:-1]=(force[1,1:]+force[1,:-1])/2
        u,v=solver.solve(fx,fy)
        phi=(1+chi[0,1:-1])/2
        b=match_logistic_width(WIDTH_REFERENCE,h)
        reference_x,reference_y,_=capillary_faces(phi,h,cfg['physical_sigma'],b)
        reference_x+=cfg['physical_force']
        ur,vr=solver.solve(reference_x,reference_y)
        raw=float(velocity[0,1:-1,:,0].mean()/h)
        files=('input.prm','geometry.raw',f'field.{step:010d}.rst',f'field.{step:010d}.rst.meta','audit.log')
        points.append(dict(step=step,raw_mean_speed=raw,
            actual_force_stokes_speed=float(u.mean()),same_phi_reference_force_speed=float(ur.mean()),
            relative_actual_force_speed_difference=float(u.mean()/raw-1),
            relative_same_phi_force_speed_difference=float(ur.mean()/raw-1),
            max_audit_state_error=state_error,
            hashes={name:sha(folder/name) for name in files}))
    return dict(height=height,name=receipt['name'],points=points)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True);parser.add_argument('--audit-driver',required=True)
    parser.add_argument('--root',required=True);parser.add_argument('--seeds',nargs='+',required=True)
    args=parser.parse_args();study=Study(args.solver,args.root)
    original=Path(args.audit_driver).resolve();driver=study.root/'.bin'/sha(original)/original.name
    driver.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(original,driver);driver.chmod(original.stat().st_mode)
    def run(path):
        seed=json.loads(Path(path).read_text())
        receipt=adjacent(study,seed,f'force-stokes-H{seed["parameters"]["height"]}',steps=2)
        for step in (receipt['start']+1,receipt['start']+2):force_audit(study,receipt,driver,step)
        result=remeasure(receipt,study.root)
        result.update(seed_receipt_sha256=sha(path),replay_receipt_name=receipt['name'])
        return result
    with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,args.seeds))
    data=dict(measurements=rows,solver_sha256=study.solver_sha,audit_driver_sha256=sha(driver),
        analysis_source_sha256=sha(__file__),
        protocol='Neutral periodic slit, prescribed physical Navier slip, actual audited force converted to MAC faces; independent same-phase-field force comparison.',
        scope='Creeping-flow hydraulic and force diagnostics. No full interface-motion validation or exact discrete equivalence implied.')
    (study.root/'force-stokes.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ('height','name')}|{'relative_actual_force_speed_difference':r['points'][0]['relative_actual_force_speed_difference']} for r in rows]))


if __name__=='__main__':main()
