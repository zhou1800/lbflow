#!/usr/bin/env python3
"""Check wall transport against an independent recurrence and replay contracts."""
import argparse
import json
from itertools import product
import math
from pathlib import Path
import shutil
import subprocess
import uuid

import numpy as np

from benchmarks import Study, compare, fields, geometry_fields, sha, spec
from capillary_response import cases
from phase_recurrence import coefficients, equilibrium_seed, VELOCITIES, WEIGHTS
from porous_media import closed_history


def bulk_step(phi, beta):
    """One-dimensional recurrence, independent of wall/population indexing."""
    gradient = (np.roll(2*phi-1, -1)-np.roll(2*phi-1, 1))/2
    segregating = np.where(abs(gradient) > 1e-9, np.sign(gradient), 0)*phi*(1-phi)
    a, b = coefficients()['bulk_diffusion'], coefficients()['bulk_segregation']
    return phi + a*(np.roll(phi, 1)+np.roll(phi, -1)-2*phi) + beta*b*(np.roll(segregating, 1)-np.roll(segregating, -1))


def ghost_audit(study, result, driver):
    """Compare actual outgoing-state gradients with a virtual full-grid stencil."""
    original = Path(result['output_files'][-1]).parent
    folder = study.root/'ghost-gradient-audit'; folder.mkdir()
    for filename in ('geometry.raw', 'field.0000000128.rst', 'field.0000000128.rst.meta'):
        shutil.copyfile(original/filename, folder/filename)
    text = (original/'input.prm').read_text().replace(str(original), str(folder))
    text = text.replace('startTime 0;', 'startTime 128;')
    inp = folder/'input.prm'; inp.write_text(text)
    completed = subprocess.run([driver, str(inp), '--capillary-state'], env=study.env,
                               cwd=folder, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, timeout=30)
    (folder/'audit.log').write_text(completed.stdout)
    if completed.returncode: raise AssertionError('Ghost state audit failed')
    data = json.loads(next(line[len('CAPILLARY_STATE '):] for line in completed.stdout.splitlines()
                           if line.startswith('CAPILLARY_STATE ')))
    color = np.array(data['color']).reshape(data['size'])
    extended = color.copy(); cot = 1/math.tan(math.radians(result['parameters']['angle']))
    for ghost, first, second in [(0, 1, 2), (-1, -2, -3)]:
        tangent = [(.5*(np.roll(color[row], -1)-np.roll(color[row], 1))) for row in (first, second)]
        extended[ghost] = color[first]-cot*abs(1.5*tangent[0]-.5*tangent[1])
    gradient = np.zeros((2,)+color.shape)
    for (dx, dy), weight in zip(VELOCITIES, WEIGHTS):
        neighbor = np.roll(extended, (-dy, -dx), axis=(0, 1))
        gradient[0] += 3*weight*dy*neighbor
        gradient[1] += 3*weight*dx*neighbor
    actual = np.array(data['gradient']).reshape(gradient.shape)
    error = float(abs(actual[:, 1:-1]-gradient[:, 1:-1]).max())
    if error > 2e-14: raise AssertionError('Published ghost continuation differs from the independent grid calculation')
    return dict(max_gradient_error=error, input_sha256=sha(inp), audit_driver_sha256=sha(driver),
                note='Outgoing-state gradients only; not a comparison with incoming-state velocities.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--audit-driver', required=True)
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--numproc-flag', default='-n')
    args = parser.parse_args()
    study = Study(args.solver, Path(args.root)/('wall-flux-'+uuid.uuid4().hex[:10]),
                  args.mpiexec, args.numproc_flag)
    transients, fixed, pure, guards, surface_transients = [], [], [], [], []
    for beta in (.35, .7, .9):
        for length in (64, 65):
            cfg = next(cases([16], [90], ['renormalized'], beta, ['csf']))
            # Test every early step, including odd times; no even-time aliasing.
            cfg.update(name=f'transient-beta{beta:g}-L{length}', size=(length, 18, 1),
                       sigma=0., steps=64, write=1, restart=32, phase_boundary='flux')
            result = study.run(cfg)
            phi = np.zeros(length); phi[16:48] = 1
            error = speed = 0.
            for time in range(65):
                if time: phi = bulk_step(phi, beta)
                _, rho, chi, velocity = geometry_fields(result, time)
                error = max(error, float(abs((chi[0, 1:-1]+1)/2-phi).max()))
                speed = max(speed, float(np.linalg.norm(velocity, axis=-1).max()))
            if error > 3e-12 or speed > 2e-14:
                raise AssertionError('Neutral wall does not preserve the complete bulk transient')
            transients.append(dict(name=result['name'], max_transient_error=error, max_speed=speed))
    for gradient, collision, beta in product(('halfway', 'ghost'), ('SRT', 'MRT'), (.35, .7, .9)):
        cfg = next(cases([16], [90], [gradient], beta, ['csf']))
        cfg.update(name=f'surface-transient-{gradient}-{collision}-beta{beta:g}', eq=collision,
                   size=(65, 18, 1), steps=64, write=1, phase_boundary='flux')
        result = study.run(cfg)
        speed = max(float(np.linalg.norm(geometry_fields(result, time)[3], axis=-1).max())
                    for time in range(65))
        if speed > 2e-14:
            raise AssertionError('A neutral surface-force transient generates resolved flow')
        surface_transients.append(dict(name=result['name'], max_speed_every_step=speed))
    for gradient in ('renormalized', 'tangent', 'halfway', 'ghost'):
        cfg = next(cases([16], [90], [gradient], capillary_models=['csf']))
        cfg.update(name='fixed-'+gradient, phase_boundary='flux', steps=129, write=64,
                   restart=64)
        seed = study.root/'equilibrium'
        phi = equilibrium_seed(seed, 65, 16, .7, tangential=False)
        cfg['size'] = (65, 18, 1)
        result = study.run(cfg, start=1, checkpoint=seed)
        _, rho, chi, velocity = geometry_fields(result)
        error = float(abs(chi[0, 1:-1]-(2*phi-1)).max())
        speed = float(np.linalg.norm(velocity, axis=-1).max())
        if error > 2e-11 or speed > 2e-14:
            raise AssertionError('Nonzero-sigma wall equilibrium was not preserved')
        fixed.append(dict(name=result['name'], color_error=error, max_speed=speed))
    # Nonneutral transient with wall-normal MPI halo exchanges and cross-rank restarts.
    cfg = next(cases([16], [60], ['halfway'], capillary_models=['csf']))
    cfg.update(name='moving-interface', phase_boundary='flux', steps=128, write=64, restart=64)
    reference = study.run(cfg)
    parallel = study.run(dict(cfg, name='moving-interface-parallel'), ranks=2)
    comparisons = [compare(reference, parallel)]
    def population_comparison(reference, other):
        a = Path(reference['output_files'][-1]).parent/'field.0000000128.rst'
        b = Path(other['output_files'][-1]).parent/'field.0000000128.rst'
        same = a.read_bytes() == b.read_bytes()
        if not same: raise AssertionError('Parallel/restarted outgoing populations differ')
        return dict(left=reference['name'], right=other['name'], bitwise_populations=same)
    population_checks = [population_comparison(reference, parallel)]
    checkpoint = Path(reference['output_files'][-1]).parent/'field.0000000064'
    for ranks in (1, 2):
        resumed = study.run(dict(cfg, name=f'moving-interface-restart-{ranks}'), ranks, 64, checkpoint)
        comparisons.append(compare(reference, resumed))
        population_checks.append(population_comparison(reference, resumed))
    rotated = dict(cfg, name='moving-interface-rotated-SRT', eq='SRT', size=(18, 64, 1),
                   label=lambda x, y, z: 255 if x in (0, 17) else (0 if 16 <= y < 48 else 100))
    unrotated = study.run(dict(cfg, name='moving-interface-SRT', eq='SRT'))
    result = study.run(rotated)
    left, right = fields(unrotated), fields(result)
    rotation_errors = {}
    for name in left:
        array = right[name].transpose(0, 2, 1, 3)
        if name == 'v': array = array[..., [1, 0]]
        rotation_errors[name] = float(abs(left[name]-array).max())
    if max(rotation_errors.values()) > 3e-12:
        raise AssertionError('Rotated wall flux differs beyond roundoff')
    mass = closed_history(reference)['max_relative_mass_drift']
    if max(mass.values()) > 1e-12: raise AssertionError('Wall flux loses component mass')
    ghost_cfg = dict(cfg, name='moving-interface-ghost', wall_gradient='ghost')
    ghost_reference = study.run(ghost_cfg)
    ghost_parallel = study.run(dict(ghost_cfg, name='moving-interface-ghost-parallel'), ranks=2)
    comparisons.append(compare(ghost_reference, ghost_parallel))
    ghost_checkpoint = Path(ghost_reference['output_files'][-1]).parent/'field.0000000064'
    ghost_resumed = study.run(dict(ghost_cfg, name='moving-interface-ghost-resumed'), 2, 64, ghost_checkpoint)
    comparisons.append(compare(ghost_reference, ghost_resumed))
    population_checks.extend([population_comparison(ghost_reference, ghost_parallel),
                              population_comparison(ghost_reference, ghost_resumed)])
    audit = ghost_audit(study, ghost_reference, str(Path(args.audit_driver).resolve()))
    for phase in (0, 1):
        for collision in ('SRT', 'MRT'):
            cfg = spec(f'pure-{phase}-{collision}', nx=12, ny=10, nz=1, dimensions=2,
                       eq=collision, tau=(.6, 1.4), steps=100, write=50,
                       wall_gradient='halfway', capillary_model='csf',
                       label=lambda x, y, z: 255 if y in (0, 9) else 75,
                       indices=[f'75 "GV" "0.000001" "1" "{phase}"', '255 "SD" "60" "N/A" "N/A"'])
            baseline = study.run(cfg)
            flux = study.run(dict(cfg, name=cfg['name']+'-flux', phase_boundary='flux'))
            pure.append(compare(baseline, flux))
    base = next(cases([16], [60], ['halfway'], capillary_models=['csf']))
    base.update(phase_boundary='flux', steps=2, write=1)
    invalid = [('bad-phase-boundary', dict(phase_boundary='unknown')),
               ('duplicate-recoloring', dict(wall_recoloring='tangential')),
               ('unsupported-gradient', dict(wall_gradient='legacy')),
               ('unsupported-weights', dict(alpha=.4)),
               ('unsupported-density', dict(density=(1., 2.))),
               ('unsupported-3d', dict(dimensions=3, size=(64, 18, 4))),
               ('open-reservoir', dict(indices=['0 "PB" "1" "0" "0"', '255 "SD" "60" "N/A" "N/A"'])),
               ('singular-angle', dict(indices=['0 "VD" "N/A" "N/A" "0"', '255 "SD" "0" "N/A" "N/A"'])),
               ('thin-fluid-layer', dict(size=(64, 4, 1), label=lambda x, y, z: 255 if y in (0, 3) else 0)),
               ('ghost-thin-layer', dict(wall_gradient='ghost', size=(64, 4, 1), label=lambda x, y, z: 255 if y in (0, 3) else 0)),
               ('ghost-singular-angle', dict(wall_gradient='ghost', indices=['0 "VD" "N/A" "N/A" "0"', '100 "VD" "N/A" "N/A" "1"', '255 "SD" "180" "N/A" "N/A"'])),
               ('wall-corner', dict(label=lambda x, y, z: 255 if y in (0, 17) or (x == 8 and y == 8) else 0))]
    for name, options in invalid:
        result = study.run(dict(base, name=name, expect_failure=True, **options))
        guards.append(dict(name=name, exit=result['exit']))
    summary = dict(neutral_transients=transients, surface_transients=surface_transients,
                   fixed_points=fixed, comparisons=comparisons, population_checks=population_checks,
                   rotation_errors=rotation_errors, mass=mass, pure_phases=pure, guards=guards,
                   ghost_gradient_audit=audit,
                   solver_sha256=study.solver_sha,
                   scope='Flat D2Q9 transient/fixed-point identities, pure-phase no-slip, short nonneutral replay; not a wetting accuracy result.')
    (study.root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__': main()
