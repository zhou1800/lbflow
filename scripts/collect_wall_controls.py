#!/usr/bin/env python3
"""Recompute wall-study summaries and verify recorded simulation inputs."""
import argparse
import json
from pathlib import Path

import numpy as np
from benchmarks import geometry_fields,measure_capillary,sha
from phase_recurrence import bulk_profile
from capillary_response import measure as meniscus_measure
from collect_capillary_controls import verify
from staggered_diagnostics import measure as staggered_measure
from wall_consistency import measure as wall_measure


def receipts(roots):
    for root in roots:
        for path in sorted((Path(root)/'results').glob('*.json')):
            r=json.loads(path.read_text());verify(r)
            if r['start'] and 'checkpoint_sha256' in r:
                folder=Path(r['output_files'][-1]).parent
                assert sha(folder/f'field.{r["start"]:010d}.rst')==r['checkpoint_sha256']
                assert sha(folder/f'field.{r["start"]:010d}.rst.meta')==r['checkpoint_metadata_sha256']
            yield path,r


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--neutral-roots',nargs='+',required=True)
    p.add_argument('--capillary-roots',nargs='+',required=True)
    p.add_argument('--bulk-root',required=True)
    p.add_argument('--staggered-roots',nargs='+',required=True)
    p.add_argument('--recurrence',required=True);p.add_argument('--fixed-summary',required=True)
    p.add_argument('--compatibility',required=True);p.add_argument('--output',required=True)
    p.add_argument('--development-receipt',help='Retain the failed odd-grid fixed point before the neutral-mask fix')
    args=p.parse_args();neutral=[];domains=[];capillary=[];bulk=[];provenance=[]
    def record(path,r,m):
        m['receipt_sha256']=sha(path)
        provenance.append(dict(name=r['name'],receipt_sha256=sha(path),input_sha256=r['input_sha256'],
                               geometry_sha256=r['geometry_sha256'],solver_sha256=r['solver_sha256'],
                               fields=[dict(name=Path(f).name,sha256=sha(f)) for f in r['output_files']]))
        return m
    for path,r in receipts(args.neutral_roots):
        if 'radius' in r['parameters']:domains.append(record(path,r,measure_capillary(r)))
        else:neutral.append(record(path,r,wall_measure(r)))
    for path,r in receipts(args.capillary_roots):capillary.append(record(path,r,meniscus_measure(r)))
    for path,r in receipts([args.bulk_root]):
        m=measure_capillary(r);m['parameters']=r['parameters'];bulk.append(record(path,r,m))
    staggered=[]
    for root in args.staggered_roots:
        for path,r in receipts([root]):
            provenance.append(dict(name=r['name'],receipt_sha256=sha(path),
                                   solver_sha256=r['solver_sha256'],checkpoint_sha256=r['checkpoint_sha256']))
        summary=json.loads((Path(root)/'staggered.json').read_text())
        for original_m in summary['measurements']:
            resumed_name=original_m['parameters']['name']+'-adjacent-1rank'
            r=json.loads((Path(root)/'results'/f'{resumed_name}.json').read_text())
            seed_command=original_m.get('seed_receipt')
            if not seed_command:raise ValueError('Staggered summary must identify its seed receipt')
            original=json.loads(Path(seed_command).read_text());verify(original)
            fresh=staggered_measure(original,r)
            for key in ('raw_speed','adjacent_pair_average_speed','adjacent_velocity_rms','two_step_velocity_rms'):
                assert fresh[key]==original_m[key]
            fresh['seed_receipt_sha256']=sha(seed_command);staggered.append(fresh)
    fixed_path=Path(args.fixed_summary);fixed=json.loads(fixed_path.read_text())
    fixed_cases=fixed_path.parent/'results'
    fixed_receipts=[]
    for path in sorted(fixed_cases.glob('fixed-*.json')):
        r=json.loads(path.read_text());verify(r);fixed_receipts.append(dict(name=r['name'],receipt_sha256=sha(path),
                       checkpoint_sha256=r.get('checkpoint_sha256'),solver_sha256=r['solver_sha256']))
    data=dict(protocol='Independent scalar recurrence, fixed-point seeds, neutral/rotated/parity controls, adjacent physical fields, bulk Laplace and exploratory nonneutral holds.',
              scope='Flat Cartesian D2Q9 neutral fixed point, equal reference densities and tau=1; no general wetting, dynamic contact-line, high-density or trapping claim.',
              recurrence=json.loads(Path(args.recurrence).read_text()),neutral=neutral,nonneutral=capillary,
              matched_bulk=bulk,bulk_domain_controls=domains,staggered=staggered,
              fixed_point_checks=fixed,compatibility=json.loads(Path(args.compatibility).read_text()),
              fixed_receipts=fixed_receipts,provenance=provenance,
              verification='Recomputed physical summaries from raw fields; input, geometry, frozen executable and supplied restart checksums verified. Staggered summaries recomputed from their seed receipts.')
    if args.development_receipt:
        r=json.loads(Path(args.development_receipt).read_text());verify(r)
        _,rho,chi,v=geometry_fields(r);cfg=r['parameters']
        reference=2*bulk_profile(cfg['size'][0],cfg['beta'])-1
        data['development_failure']=dict(case=r['name'],solver_sha256=r['solver_sha256'],
            input_sha256=r['input_sha256'],geometry_sha256=r['geometry_sha256'],
            checkpoint_sha256=r['checkpoint_sha256'],max_speed=float(np.linalg.norm(v,axis=-1).max()),
            max_color_error=float(abs(chi[0,1:-1]-reference).max()),
            max_wall_color_difference=float(abs(chi[0,1]-chi[0,cfg['height']//2+1]).max()),
            source_patch='run-kernels/wall-pre-neutral-norm.patch',
            scope='Development failure before neutral tangent-magnitude/cutoff consistency; not a final validated equilibrium.')
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(neutral=len(neutral),nonneutral=len(capillary),bulk=len(bulk),staggered=len(staggered),
                         fixed=len(fixed['fixed_points']))))


if __name__=='__main__':main()
