#!/usr/bin/env python3
"""Extend completed slit controls and retain their entire measured histories."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

from benchmarks import Study,sha
from capillary_response import measure


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True)
    parser.add_argument('--initial-roots',nargs='+',required=True)
    parser.add_argument('--root',required=True)
    parser.add_argument('--heights',nargs='+',type=int,default=[16,32,48])
    parser.add_argument('--angles',nargs='+',type=float,default=[60])
    parser.add_argument('--walls',nargs='+',default=['renormalized'])
    parser.add_argument('--betas',nargs='+',type=float)
    parser.add_argument('--phase-boundaries',nargs='+',choices=['bounce','flux'],default=['bounce','flux'])
    parser.add_argument('--factor',type=int,default=2)
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    if args.factor<2:parser.error('The continuation factor must be at least two')
    seeds=[]
    identities=set()
    for root in args.initial_roots:
        for path in sorted((Path(root)/'results').glob('*.json')):
            receipt=json.loads(path.read_text())
            if 'parameters' not in receipt:continue  # Separate measurement JSON is not a run receipt.
            p=receipt['parameters']
            if p.get('height') not in args.heights or p.get('angle') not in args.angles or p.get('wall_gradient') not in args.walls:
                continue
            if args.betas is not None and p['beta'] not in args.betas:continue
            if p.get('phase_boundary','bounce') not in args.phase_boundaries:continue
            identity=(p['height'],p['angle'],p['wall_gradient'],p.get('capillary_model','perturbation'),
                      p['beta'],p.get('phase_boundary','bounce'),p.get('wall_recoloring','legacy'),
                      p.get('curvature_stencil','normal'),p['eq'],tuple(p['tau']),tuple(p['size']),p['sigma'],
                      p.get('recoloring_model','cosine'),p.get('mrt_spectrum','legacy'))
            if identity in identities:raise ValueError('Duplicate seed configuration')
            if receipt['exit']!=0 or not receipt['finished']:raise ValueError('Incomplete seed')
            identities.add(identity)
            seeds.append((path,receipt))
    if not seeds:parser.error('No completed matching seeds')
    study=Study(args.solver,args.root)

    def extend(seed):
        path,original=seed
        p=original['parameters']
        folder=Path(original['output_files'][-1]).parent
        cut=p['steps']
        new=dict(p,name=p['name']+'-long-hold',steps=args.factor*cut,restart=args.factor*cut,
                 geometry_file=str(folder/'geometry.raw'))
        result=study.run(new,original['ranks'],cut,folder/f'field.{cut:010d}')
        result['output_files']=original['output_files']+result['output_files']
        result['seed_receipt_sha256']=sha(path)
        result['seed_solver_sha256']=original['solver_sha256']
        (study.root/'results'/f'{result["name"]}.json').write_text(json.dumps(result,indent=2)+'\n')
        m=measure(result)
        m['seed_steps']=cut
        m['seed_pressure_coefficient']=measure(original)['history'][-1]['pressure_coefficient']
        m['hold_pressure_coefficient_change']=m['history'][-1]['pressure_coefficient']-m['seed_pressure_coefficient']
        return m

    with ThreadPoolExecutor(max_workers=args.jobs) as pool: measurements=list(pool.map(extend,seeds))
    (study.root/'capillary-response.json').write_text(json.dumps(dict(protocol='Restart the same outgoing state with unchanged physical parameters; include seed fields in the closed-mass and pressure audit.',measurements=measurements),indent=2)+'\n')
    for m in measurements:
        print(json.dumps(dict(name=m['name'],pressure_error=m['relative_pressure_error'],
                              hold_coefficient_change=m['hold_pressure_coefficient_change'],
                              late_color=m['history'][-1]['color_rms_change'])))


if __name__=='__main__':main()
