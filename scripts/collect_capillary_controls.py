#!/usr/bin/env python3
"""Reanalyze and checksum the separate closed-capillary experiment campaigns."""
import argparse
import json
from pathlib import Path

import numpy as np
from benchmarks import geometry_fields,measure_capillary,sha
from capillary_response import measure
from porous_media import closed_history


def verify(receipt):
    folder=Path(receipt['output_files'][-1]).parent
    assert receipt['exit']==0 and receipt['finished']
    assert sha(folder/'input.prm')==receipt['input_sha256']
    assert sha(folder/'geometry.raw')==receipt['geometry_sha256']
    executable=next(Path(part) for part in receipt['command'] if Path(part).name=='LBFlow')
    assert sha(executable)==receipt['solver_sha256']


def slit_roots(roots):
    measurements=[]
    for root in roots:
        for path in sorted((Path(root)/'results').glob('*.json')):
            r=json.loads(path.read_text())
            verify(r)
            m=measure(r)
            if 'seed_steps' in r['parameters']:raise ValueError('Unexpected seed protocol')
            if r['start']:
                m['seed_steps']=r['start']
                old=next(h for h in m['history'] if h['step']==r['start'])
                m['seed_pressure_coefficient']=old['pressure_coefficient']
                m['hold_pressure_coefficient_change']=m['history'][-1]['pressure_coefficient']-old['pressure_coefficient']
            m['receipt_sha256']=sha(path)
            measurements.append(m)
    return measurements


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial-roots',nargs='+',required=True)
    parser.add_argument('--hold-roots',nargs='+',required=True)
    parser.add_argument('--bulk-initial-root',required=True)
    parser.add_argument('--bulk-root',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    initial=slit_roots(args.initial_roots)
    holds=slit_roots(args.hold_roots)
    bulk=json.loads((Path(args.bulk_root)/'bulk-and-recolor.json').read_text())
    bulk_receipts={}
    for root in (args.bulk_initial_root,args.bulk_root):
        for path in (Path(root)/'results').glob('*.json'):
            receipt=json.loads(path.read_text())
            verify(receipt)
            bulk_receipts[receipt['name']]=receipt
    for m in bulk['measurements']:
        path=Path(args.bulk_root)/'results'/f'{m["name"]}.json'
        if not path.exists():path=Path(args.bulk_initial_root)/'results'/f'{m["name"]}.json'
        r=json.loads(path.read_text())
        verify(r)
        _,rho,chi,v=geometry_fields(r)
        fresh=measure_capillary(r) if m['name'].startswith('laplace') else dict(
            max_spurious_speed=float(np.linalg.norm(v,axis=-1).max()),
            max_wall_color_difference_from_center=float(np.max(abs(chi[0,1]-chi[0,8]))))
        fresh.update(**closed_history(r))
        # Independently reproduce the reported scalar/profile diagnostics.
        for key in ['pressure_jump','relative_laplace_error','max_spurious_speed',
                    'max_wall_color_difference_from_center']:
            if key in fresh and fresh[key]!=m[key]:raise AssertionError('Bulk reanalysis mismatch')
        m['receipt_sha256']=sha(path)
    output=dict(protocol='Direct closed-slit pressure, static neutral-wall flow, signed 60/120 controls, three widths and longer holds; bulk-drop and zero-sigma controls.',
                scope='Equal density and viscosity, tau=1, beta=0.7. Joint grid/sharp-interface sequence. No entry-threshold, trapping law, fixed-width convergence, high-density or curved-3D-wetting claim.',
                planar_width_09=4.765235021842287,viscosity=1/6,
                measurements=initial,long_holds=holds,bulk_controls=bulk['measurements'],
                initial_drop_measurements=bulk['initial_drop_measurements'],
                physical_executions=len(initial)+len(holds)+len(bulk_receipts),
                verification='All included final measurements reanalyzed from raw fields; input, geometry and frozen executable hashes verified.')
    target=Path(args.output)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(dict(initial=len(initial),holds=len(holds),bulk=len(bulk['measurements']),
                          physical_executions=output['physical_executions'])))


if __name__=='__main__':main()
