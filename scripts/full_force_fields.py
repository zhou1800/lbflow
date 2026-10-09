#!/usr/bin/env python3
"""Compare complete recorded force/velocity fields, retaining their gauges.

Capillary-force norms and the velocity they drive are different observables:
a gradient-force difference may change pressure without driving Stokes flow.
The MAC-to-cell average is specified, not fitted to the LBM fields.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from benchmarks import geometry_fields, sha
from collect_capillary_controls import verify
from collect_phase_kinetics import receipt_at_root
from stokes_phase_reference import StokesSlit, capillary_faces, match_logistic_width
from front_mobility import WIDTH_REFERENCE


def norm_difference(left,right):
    denominator=sum(float(np.sum(a*a)) for a in left)
    error=sum(float(np.sum((a-b)**2)) for a,b in zip(left,right))
    return float(np.sqrt(error/denominator)) if denominator else None


def collect(phase_root,output):
    output.mkdir(parents=True,exist_ok=False);rows=[]
    campaign=phase_root/'hydraulic-audit-canonical'
    for path in sorted((campaign/'results').glob('*.json')):
        receipt=receipt_at_root(json.loads(path.read_text()),phase_root);verify(receipt)
        p=receipt['parameters'];height=p['height'];h=p['dx'];solver=StokesSlit(height,slip=p['physical_slip'])
        for step in (receipt['start']+1,receipt['start']+2):
            log=campaign/'audit'/receipt['name']/str(step)/'audit.log'
            text=log.read_text();state=json.loads(next(line[16:] for line in text.splitlines() if line.startswith('CAPILLARY_STATE ')))
            shape=tuple(state['size']);density=np.array(state['density']).reshape(shape)
            color=np.array(state['color']).reshape(shape);_,rho,chi,velocity=geometry_fields(receipt,step)
            state_error=max(float(abs(density-rho[0]).max()),float(abs(color[density>0]-chi[0][density>0]).max()))
            if state_error>3e-12:raise AssertionError('Audited state differs from physical output')
            force=np.array(state['force']).reshape((2,)+shape)[::-1,1:-1]/h**3
            fx=(force[0]+np.roll(force[0],1,axis=1))/2
            fy=np.zeros((height+1,4*height));fy[1:-1]=(force[1,1:]+force[1,:-1])/2
            phi=(1+chi[0,1:-1])/2;b=match_logistic_width(WIDTH_REFERENCE,h)
            rx,ry,_=capillary_faces(phi,h,p['physical_sigma'],b)
            u,v=solver.solve(fx+p['physical_force'],fy);ur,vr=solver.solve(rx+p['physical_force'],ry)
            raw=np.moveaxis(velocity[0,1:-1],-1,0)/h
            actual=np.array([(u+np.roll(u,-1,axis=1))/2,(v[1:]+v[:-1])/2])
            reference=np.array([(ur+np.roll(ur,-1,axis=1))/2,(vr[1:]+vr[:-1])/2])
            row=dict(height=height,step=step,max_state_error=state_error,
                relative_actual_force_cell_velocity_error=norm_difference(raw,actual),
                relative_reference_force_cell_velocity_error=norm_difference(raw,reference),
                relative_capillary_face_force_difference=norm_difference((fx,fy),(rx,ry)),
                relative_force_driven_mac_velocity_difference=norm_difference((u,v),(ur,vr)),
                max_actual_cell_velocity_difference=float(abs(raw-actual).max()),
                max_reference_cell_velocity_difference=float(abs(raw-reference).max()),
                lbm_mean_speed=float(raw[0].mean()),actual_force_stokes_mean=float(u.mean()),reference_force_stokes_mean=float(ur.mean()),
                audit_log_sha256=sha(log),receipt_sha256=sha(path),solver_sha256=receipt['solver_sha256'])
            rows.append(row)
            np.savez_compressed(output/f'H{height}-step{step}.npz',lbm_velocity=raw,actual_force_velocity=actual,
                reference_force_velocity=reference,actual_fx=fx,actual_fy=fy,reference_fx=rx,reference_fy=ry,h=h)
    if len(rows)!=6:raise AssertionError('Expected two audited states at each of three grids')
    result=dict(measurements=rows,comparison='Cell LBM versus specified two-face MAC averages; forces compared on MAC faces. Gradient-force gauges are retained.',
        scope='Hydraulic fields on supplied recorded phases; not independent validation of their evolution.')
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(rows),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--phase-root',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();collect(Path(args.phase_root).resolve(),Path(args.output).resolve())


if __name__=='__main__':main()
