#!/usr/bin/env python3
"""Recompute the wall-flux/location study and retain raw-input provenance."""
import argparse
import json
from pathlib import Path

import numpy as np

from benchmarks import geometry_fields, sha
from capillary_response import measure
from collect_capillary_controls import verify
from phase_recurrence import bulk_profile
from porous_media import closed_history
from staggered_diagnostics import measure as adjacent_measure


def neutral_measure(receipt):
    cfg = receipt['parameters']; height = cfg['height']
    _, rho, chi, velocity = geometry_fields(receipt)
    phi = (chi[0, 1:-1]+1)/2
    reference = bulk_profile(cfg['size'][0], cfg['beta'], cfg['steps'])
    center = .5*(phi[height//2-1]+phi[height//2])
    audit = closed_history(receipt)
    return dict(name=receipt['name'], height=height, beta=cfg['beta'], sigma=cfg['sigma'],
                wall_gradient=cfg['wall_gradient'], phase_boundary=cfg.get('phase_boundary', 'bounce'),
                max_speed=float(np.linalg.norm(velocity, axis=-1).max()),
                max_transverse_color_difference=float(abs(2*(phi-center)).max()),
                max_bulk_phase_error=float(abs(phi-reference).max()),
                wall_phase_profile=phi[0].tolist(), center_phase_profile=center.tolist(),
                mass=audit['max_relative_mass_drift'], history=audit['history'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roots', nargs='+', required=True)
    parser.add_argument('--independent', required=True)
    parser.add_argument('--checks', required=True)
    parser.add_argument('--comparisons', required=True)
    parser.add_argument('--adjacent-root', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    neutral, initial, holds, provenance = [], [], [], []
    for root_name in args.roots:
        root = Path(root_name)
        for path in sorted((root/'results').glob('*.json')):
            receipt = json.loads(path.read_text())
            if 'parameters' not in receipt: continue
            verify(receipt)
            cfg = receipt['parameters']; folder = Path(receipt['output_files'][-1]).parent
            if receipt['start']:
                for key, suffix in [('checkpoint_sha256', '.rst'), ('checkpoint_metadata_sha256', '.rst.meta')]:
                    if sha(folder/f'field.{receipt["start"]:010d}{suffix}') != receipt[key]:
                        raise AssertionError('Continuation seed checksum changed')
            result = neutral_measure(receipt) if cfg['angle'] == 90 else measure(receipt)
            result.update(solver_sha256=receipt['solver_sha256'], input_sha256=receipt['input_sha256'],
                          geometry_sha256=receipt['geometry_sha256'], receipt_sha256=sha(path))
            if cfg['angle'] == 90: neutral.append(result)
            elif receipt['start']:
                result['seed_steps'] = receipt['start']
                seed = next(item for item in result['history'] if item['step'] == receipt['start'])
                result['hold_pressure_coefficient_change'] = result['history'][-1]['pressure_coefficient']-seed['pressure_coefficient']
                holds.append(result)
            else: initial.append(result)
            provenance.append(dict(root=root.name, case=receipt['name'], receipt_sha256=sha(path),
                solver_sha256=receipt['solver_sha256'], input_sha256=receipt['input_sha256'],
                geometry_sha256=receipt['geometry_sha256'], start=receipt['start'],
                checkpoint_sha256=receipt.get('checkpoint_sha256'),
                checkpoint_metadata_sha256=receipt.get('checkpoint_metadata_sha256'),
                seed_receipt_sha256=receipt.get('seed_receipt_sha256'),
                fields=[dict(path=f'{Path(filename).parent.parent.parent.name}/cases/{Path(filename).parent.name}/{Path(filename).name}',
                             sha256=sha(filename)) for filename in receipt['output_files']]))
    adjacent_root = Path(args.adjacent_root)
    adjacent = []
    for item in json.loads((adjacent_root/'staggered.json').read_text())['measurements']:
        seed_path = Path(item['seed_receipt']); seed = json.loads(seed_path.read_text()); verify(seed)
        resumed_path = adjacent_root/'results'/f'{seed["parameters"]["name"]}-adjacent-1rank.json'
        resumed = json.loads(resumed_path.read_text()); verify(resumed)
        fresh = adjacent_measure(seed, resumed)
        for key in ('raw_speed', 'adjacent_pair_average_speed', 'adjacent_velocity_rms', 'two_step_velocity_rms'):
            if fresh[key] != item[key]: raise AssertionError('Adjacent raw-field reanalysis changed')
        fresh['seed_receipt_sha256'] = sha(seed_path); fresh['receipt_sha256'] = sha(resumed_path)
        adjacent.append(fresh)
    data = dict(protocol='Closed equal-density D2Q9 flat slits: conservative relative-phase transfer, center/tangent/halfway/ghost color gradients, signed pressure and four-duration holds. H16/32/48 use joint grid/sharp-interface scaling.',
                scope='Neutral complete transient and fixed-point controls; bounded static nonneutral pressure comparison. No general wetting, open-reservoir, curved-wall, thermodynamic or trapping claim.',
                independent=json.loads(Path(args.independent).read_text()),
                contracts=json.loads(Path(args.checks).read_text()),
                prior_pilot_confirmation=json.loads(Path(args.comparisons).read_text()),
                neutral=neutral, nonneutral=initial, long_holds=holds, adjacent_fields=adjacent, provenance=provenance,
                physical_executions=len(provenance),
                verification='Every included measurement recomputed from raw fields; input, geometry, immutable executable and continuation seed checksums verified. Historical hashes retained by source kernel.')
    output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(dict(neutral=len(neutral), initial=len(initial), holds=len(holds), physical_executions=len(provenance))))


if __name__ == '__main__': main()
