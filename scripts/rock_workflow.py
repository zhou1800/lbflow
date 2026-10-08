#!/usr/bin/env python3
"""Run a bounded digital-rock application with mass, late-time and MPI audits.

An evolving two-phase flux is reported as an advective flux. It is never
automatically promoted to relative permeability or a measured rock property.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np

from benchmarks import Study, compare, geometry_fields, sha, spec
from porous_media import closed_history


def rock_case(dataset, name, model='scsp', tau=(1.,), force=2e-6,
              phase=0, steps=4000, angle=60, two_phase=False):
    if steps < 4 or steps % 4:
        raise ValueError('The base run duration must be a positive multiple of four')
    dataset = Path(dataset).resolve()
    manifest = json.loads((dataset/'manifest.json').read_text())
    refinement = manifest['refinement']
    geometry = 'two-phase.raw' if two_phase else 'single-phase.raw'
    path = dataset/geometry
    if sha(path) != manifest['files'][geometry]['sha256']:
        raise ValueError('Prepared image has changed since its provenance was recorded')
    lattice_force = force/refinement**3
    nx, ny, nz = manifest['simulation_size_xyz']
    return spec(name, nx=nx, ny=ny, nz=nz, model=model, tau=tau,
        steps=steps*refinement**2, write=(steps//4)*refinement**2,
        restart=(steps//2)*refinement**2, force=lattice_force,
        geometry_file=str(path), sigma=.002/refinement,
        source_manifest=manifest, contact_angle=angle,
        two_phase=two_phase, timeout=7200,
        indices=[f'75 "GV" "{lattice_force:.17g}" "1" "{phase}"']+
                ([f'76 "GV" "{lattice_force:.17g}" "1" "1"'] if two_phase else [])+
                [f'255 "SD" "{angle}" "N/A" "N/A"'])


def measure(result):
    p = result['parameters']
    audit = closed_history(result)
    before, final = audit['history'][-2:]
    a = geometry_fields(result, before['step'])
    b = geometry_fields(result)
    mask = b[1] > 0
    velocity_difference = np.linalg.norm((b[3]-a[3])[mask])
    velocity_scale = max(np.linalg.norm(b[3][mask]), 1e-30)
    changes = dict(relative_velocity_l2=float(velocity_difference/velocity_scale),
                   darcy_velocity_absolute=(np.array(final['darcy_velocity'])-
                       np.array(before['darcy_velocity'])).tolist())
    if a[2] is not None:
        delta = (b[2]-a[2])[mask]
        changes.update(color_rms=float(np.sqrt(np.mean(delta**2))),
                       color_max=float(np.max(np.abs(delta))))
    record = dict(name=result['name'], model=p['model'], tau=p['tau'],
        force=p['force'], sigma=p['sigma'], contact_angle=p['contact_angle'],
        two_phase=p['two_phase'], dataset=p['source_manifest'], **audit,
        late_field_changes=changes, timing=result['timing'],
        solver_sha256=result['solver_sha256'], geometry_sha256=result['geometry_sha256'])
    if not p['two_phase'] and p['force']:
        phase = 0 if p['model'] == 'scsp' else int(p['indices'][0].split()[-1].strip('"'))
        nu = (p['tau'][phase]-.5)/3
        k = nu*final['darcy_velocity'][0]/p['force']
        record.update(K_lattice=k, K_m2=k*p['source_manifest']['simulation_voxel_size_m']**2,
                      permeability_scope='Single-phase body-force estimate on the specified '
                                         'periodically repeated crop, equal unit density.')
    else:
        record['permeability_scope'] = 'No relative permeability inferred from this transient run.'
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--mode', choices=('single', 'two', 'refinement', 'contracts'), default='single')
    parser.add_argument('--steps', type=int)
    args = parser.parse_args()
    study = Study(args.solver, args.root)
    # Independent cases may run simultaneously. Let the scheduler distribute
    # functional tests; dedicated timings use explicit per-rank CPU affinity.
    study.env.setdefault('OMPI_MCA_hwloc_base_binding_policy', 'none')
    steps = args.steps or (16000 if args.mode == 'two' else 8000 if args.mode == 'refinement' else 4000)
    cases = []
    if args.mode == 'single':
        for force in (1e-6, 2e-6, 4e-6):
            cases.append(rock_case(args.dataset, f'rock-scsp-F{force}', force=force, steps=steps))
        # Pure-phase reductions exercise both local viscosity branches on a rock.
        for phase, tau in ((0, .6), (1, 1.4)):
            phase_steps = max(steps, 12000 if tau == .6 else 4000)
            cases.extend((rock_case(args.dataset, f'rock-scsp-tau{tau}', tau=(tau,), steps=phase_steps),
                rock_case(args.dataset, f'rock-color-pure{phase}-tau{tau}', model='color',
                          tau=(.6, 1.4), phase=phase, steps=phase_steps)))
    elif args.mode == 'two':
        cases = [rock_case(args.dataset, f'rock-two-theta{angle}-F{force}', model='color',
            tau=(1., 1.), force=force, angle=angle, two_phase=True, steps=steps)
            for angle, force in ((60, 0.), (60, 2e-6), (120, 2e-6))]
    elif args.mode == 'refinement':
        cases = [rock_case(args.dataset, 'rock-refined-scsp', steps=4000)]
        cases += [rock_case(args.dataset, 'rock-refined-two-F2e-06', model='color',
                  tau=(1., 1.), force=2e-6, two_phase=True, steps=steps)]
    else:
        case = rock_case(args.dataset, 'rock-contract-continuous', model='color',
                        tau=(1., 1.), two_phase=True, steps=400)
        continuous = study.run(case)
        parallel = study.run(dict(case, name='rock-contract-parallel'), 2)
        cut = case['restart']
        checkpoint = Path(continuous['output_files'][-1]).parent/f'field.{cut:010d}'
        resumed = study.run(dict(case, name='rock-contract-resumed'), 2, cut, checkpoint)
        summary = dict(comparisons=[compare(continuous, parallel), compare(continuous, resumed)],
                       measurement=measure(continuous))
        (study.root/'rock-contracts.json').write_text(json.dumps(summary, indent=2)+'\n')
        return
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda case: study.run(case, args.ranks), cases))
    measurements = [measure(result) for result in results]
    comparisons = []
    if args.mode == 'single':
        for phase, tau in ((0, .6), (1, 1.4)):
            scsp = next(r for r in results if r['parameters']['name'] == f'rock-scsp-tau{tau}')
            color = next(r for r in results if r['parameters']['name'] == f'rock-color-pure{phase}-tau{tau}')
            left, right = geometry_fields(scsp), geometry_fields(color)
            density_error = float(np.max(np.abs(left[1]-right[1])))
            raw_velocity_error = float(np.max(np.abs(left[3]-right[3])))
            # SCSP uses the legacy reference-density (incompressible) momentum
            # convention; color fluid divides by local total density. Compare
            # like velocity definitions instead of demanding raw equality.
            normalized = left[3]/np.where(left[1] > 0, left[1], 1.)[..., None]
            velocity_error = float(np.max(np.abs(normalized-right[3])))
            if max(density_error, velocity_error) > 2e-11:
                raise AssertionError('Single-phase reduction on the rock disagrees')
            comparisons.append(dict(phase=phase, tau=tau,
                max_density_difference=density_error, max_velocity_difference=velocity_error,
                raw_velocity_difference=raw_velocity_error,
                comparison_definition='SCSP velocity divided by SCSP density (unit reference '
                                      'density), compared with color-fluid mixture velocity.'))
    (study.root/f'rock-{args.mode}.json').write_text(json.dumps(
        dict(measurements=measurements, pure_phase_comparisons=comparisons), indent=2)+'\n')


if __name__ == '__main__':
    main()
