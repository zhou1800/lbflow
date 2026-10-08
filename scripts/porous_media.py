#!/usr/bin/env python3
"""Independent pore-flow measurements and analytic two-phase pore benchmarks.

Phase flux here means the phase-fraction-weighted *advective* Darcy flux of
the common mixture velocity. It is not a measurement of the recoloring flux.
All volume averages include solid voxels, whose velocities are zero.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path

import numpy as np

from benchmarks import Study, geometry_fields, spec


def pore_measurement(result, at=None):
    arrays, rho, chi, velocity = geometry_fields(result, at)
    pore = rho > 0
    if not pore.any():
        raise AssertionError('No fluid voxels in the output')
    speed = np.linalg.norm(velocity, axis=-1)
    measurement = dict(darcy_velocity=velocity.mean(axis=(0, 1, 2)).tolist(),
                       porosity=float(pore.mean()),
                       total_mass=float(rho.sum()),
                       density_range=[float(rho[pore].min()), float(rho[pore].max())],
                       max_speed=float(speed.max()),
                       max_mach=float(speed.max()*math.sqrt(3)))
    if chi is not None:
        phi = np.divide(arrays['r0'][..., 0], rho, out=np.zeros_like(rho), where=pore)
        fractions = (phi, pore.astype(float)-phi)
        measurement.update(
            phase_masses=[float(arrays[f'r{i}'].sum()) for i in (0, 1)],
            phase_saturations=[float(fraction[pore].mean()) for fraction in fractions],
            phase_advective_darcy_flux=[(fraction[..., None]*velocity).mean(
                axis=(0, 1, 2)).tolist() for fraction in fractions],
            diffuse_pore_fraction=float((np.abs(chi[pore]) < .9).mean()))
        flux = np.array(measurement['phase_advective_darcy_flux']).sum(axis=0)
        if not np.allclose(flux, measurement['darcy_velocity'], rtol=0, atol=1e-14):
            raise AssertionError('Phase fluxes do not sum to mixture Darcy velocity')
    return measurement


def closed_history(result):
    history = []
    for name in result['output_files']:
        at = int(Path(name).name.split('.')[1])
        history.append(dict(step=at, **pore_measurement(result, at)))
    initial = history[0]
    if initial['step'] != 0:
        raise ValueError('Closed mass audit needs the initial snapshot')
    keys = ['total_mass'] + (['phase_masses'] if 'phase_masses' in initial else [])
    mass_errors = {}
    for key in keys:
        reference = np.asarray(initial[key])
        errors = [np.abs(np.asarray(item[key])-reference)/np.maximum(reference, 1.)
                  for item in history]
        mass_errors[key] = float(np.max(errors))
        if mass_errors[key] > 5e-9:
            raise AssertionError(f'{result["name"]}: closed {key} drift {mass_errors[key]}')
    return dict(history=history, max_relative_mass_drift=mass_errors)


def bundle_cases(scales=(1, 2, 3), dimensions=(2, 3)):
    cases = []
    for ndim in dimensions:
        for scale in scales:
            h, ny = 8*scale, 36*scale
            widths = (h, 2*h)
            force = 1e-6/scale**3
            for taus in ((1.,), (.6, 1.4), (1.4, .6)):
                model = 'scsp' if len(taus) == 1 else 'color'
                viscosities = [(tau-.5)/3 for tau in taus]
                decay = max(widths[i]**2/(viscosities[i % len(taus)]*math.pi**2)
                            for i in (0, 1))
                steps = math.ceil(18*decay/1000)*1000

                def label(x, y, z, h=h, scale=scale, model=model):
                    if 4*scale <= y < 4*scale+h:
                        return 75
                    if 8*scale+h <= y < 8*scale+3*h:
                        return 75 if model == 'scsp' else 76
                    return 255

                cases.append(spec(
                    f'bundle-d{ndim}-s{scale}-{model}-tau'+'-'.join(map(str, taus)),
                    nx=4*scale, ny=ny, nz=1 if ndim == 2 else 4*scale,
                    dimensions=ndim, model=model, tau=taus, steps=steps, write=steps//2,
                    timeout=3600,
                    force=force, scale=scale, widths=widths, label=label,
                    indices=[f'75 "GV" "{force:.17g}" "1" "0"']+
                            ([f'76 "GV" "{force:.17g}" "1" "1"'] if model == 'color' else [])+
                            ['255 "SD" "90" "N/A" "N/A"']))
    return cases


def bundle_study(study, jobs=1, scales=(1, 2, 3), dimensions=(2, 3)):
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        results = list(pool.map(study.run, bundle_cases(scales, dimensions)))
    absolute = {}
    measured = []
    for result in results:
        p = result['parameters']
        audit = closed_history(result)
        final = audit['history'][-1]
        ny = p['size'][1]
        widths = p['widths']
        # A midpoint sample of u(y)=F*y*(H-y)/(2*mu) has this exact sum.
        weights = np.array([h**3/12+h/24 for h in widths])/ny
        continuous = sum(h**3/12 for h in widths)/ny
        expected_k = float(weights.sum())
        record = dict(name=result['name'], dimensions=p['dimensions'], scale=p['scale'],
                      model=p['model'], tau=p['tau'], expected_K_midpoint=expected_k,
                      expected_K_continuum=continuous, **audit,
                      timing=result['timing'], solver_sha256=result['solver_sha256'])
        if p['model'] == 'scsp':
            nu = (p['tau'][0]-.5)/3
            k = nu*final['darcy_velocity'][0]/p['force']
            if abs(k/expected_k-1) > 2e-6:
                raise AssertionError(f'{result["name"]}: analytic permeability disagreement')
            absolute[(p['dimensions'], p['scale'])] = k
            record.update(measured_K=k, relative_midpoint_error=abs(k/expected_k-1),
                          relative_continuum_error=abs(k/continuous-1))
        else:
            q = np.array(final['phase_advective_darcy_flux'])[:, 0]
            nu = (np.array(p['tau'])-.5)/3
            expected_q = p['force']*weights/nu
            if np.max(np.abs(q/expected_q-1)) > 2e-6:
                raise AssertionError(f'{result["name"]}: incorrect phase flux/viscosity')
            record.update(phase_flux=q.tolist(), analytic_phase_flux=expected_q.tolist(),
                          relative_phase_flux_errors=np.abs(q/expected_q-1).tolist(),
                          reference_relative_permeability=(weights/weights.sum()).tolist())
        measured.append(record)
    for record in measured:
        if record['model'] == 'color':
            baseline = absolute[(record['dimensions'], record['scale'])]
            nu = (np.array(record['tau'])-.5)/3
            force = 1e-6/record['scale']**3
            record['measured_relative_permeability'] = (
                nu*np.array(record['phase_flux'])/(force*baseline)).tolist()
            record['single_phase_K_baseline'] = baseline
    summary = dict(protocol='Two isolated parallel slit pores; periodic x (and z in 3D), '
                            'halfway no-slip walls; equal unit density; no fluid-fluid interface.',
                   flux_definition='Volume average of phase fraction times mixture velocity; '
                                   'advective flux, including solid voxels in the averaging volume.',
                   measurements=measured)
    (study.root/'bundle.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


def entry_cases(scales=(1, 2, 4), walls=('legacy', 'renormalized'), dimensions=(2,)):
    cases = []
    if len(dimensions) > 1:
        return [case for ndim in dimensions for case in entry_cases(scales, walls, (ndim,))]
    ndim = dimensions[0]
    for scale in scales:
        h, nx = 16*scale, 64*scale+1
        sigma = .005/scale
        # Contact angle is measured through the initially resident phase 1.
        pc = (2 if ndim == 2 else 4)*sigma*math.cos(math.radians(60))/h
        for wall in walls:
            for ratio in (.5, 1.5):
                dp = ratio*pc

                def label(x, y, z, h=h, nx=nx, ndim=ndim):
                    return (255 if y in (0, h+1) or (ndim == 3 and z in (0, h+1)) else 50 if x == 0 else
                            200 if x == nx-1 else 0 if x < h else 100)

                cases.append(spec(
                    f'entry-'+('d3-' if ndim == 3 else '')+f's{scale}-{wall}-ratio{ratio}',
                    nx=nx, ny=h+2, nz=1 if ndim == 2 else h+2,
                    dimensions=ndim, steps=8000*scale**2, write=2000*scale**2,
                    timeout=3600,
                    tau=(1., 1.), sigma=sigma, scale=scale, height=h, pressure_drop=dp,
                    sharp_entry_pressure=pc, pressure_ratio=ratio, wall_gradient=wall,
                    label=label, indices=[f'50 "PB" "{1+3*dp:.17g}" "1" "0"',
                        f'200 "PB" "1" "{3 if ndim == 2 else 2}" "1"', '0 "VD" "N/A" "N/A" "0"',
                        '100 "VD" "N/A" "N/A" "1"', '255 "SD" "60" "N/A" "N/A"']))
    return cases


def entry_study(study, jobs=1, scales=(1, 2, 4), walls=('legacy', 'renormalized'), dimensions=(2,)):
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        results = list(pool.map(study.run, entry_cases(scales, walls, dimensions)))
    measured = []
    for result in results:
        p = result['parameters']
        history = []
        for path in result['output_files']:
            step = int(Path(path).name.split('.')[1])
            a, rho, chi, velocity = geometry_fields(result, step)
            phi = np.divide(a['r0'][..., 0], rho, out=np.zeros_like(rho), where=rho > 0)
            interior = (0 if p['dimensions'] == 2 else slice(1, -1), slice(1, -1), slice(1, -1))
            history.append(dict(step=step,
                phase0_equivalent_length=float(phi[interior].sum()/p['height']**(p['dimensions']-1)),
                pore_mean_ux=float(velocity[interior+(0,)].mean()),
                **pore_measurement(result, step)))
        late = history[-1]['phase0_equivalent_length']-history[-2]['phase0_equivalent_length']
        if (p['pressure_ratio']-.999)*late <= 0:
            raise AssertionError(f'{result["name"]}: incorrect late displacement direction')
        # Diagnostic inference assumes a developed slit profile. Contact-line
        # dissipation and diffuse walls can violate it; it is not a static Pc measurement.
        mu = 1/6
        if p['dimensions'] == 2:
            conductance = (p['height']**2+.5)/(12*mu*(p['size'][0]-1))
        else:
            odd = np.arange(1, 400, 2, dtype=float)
            square_coefficient = (1-192/math.pi**5*np.sum(np.tanh(odd*math.pi/2)/odd**5))/12
            conductance = square_coefficient*p['height']**2/(mu*(p['size'][0]-1))
        apparent_pc = p['pressure_drop']-history[-1]['pore_mean_ux']/conductance
        measured.append(dict(name=result['name'], dimensions=p['dimensions'], scale=p['scale'], height=p['height'],
            wall_gradient=p['wall_gradient'], sigma=p['sigma'], pressure_ratio=p['pressure_ratio'],
            sharp_entry_pressure=p['sharp_entry_pressure'], pressure_drop=p['pressure_drop'],
            history=history, late_length_change=late, apparent_entry_pressure=apparent_pc,
            apparent_pressure_relative_error=apparent_pc/p['sharp_entry_pressure']-1,
            solver_sha256=result['solver_sha256'], timing=result['timing']))
    summary = dict(protocol='Equal viscosity, 60-degree slit (2D) or square tube (3D); open pressure boundaries, '
                            'phase 0 on the left, phase 1 on the right; diffusive refinement.',
                   mass_scope='Open system: mass change is permitted and is not a closed mass audit.',
                   apparent_pressure_scope='Hydraulic inference from mean velocity; excludes no '
                                           'contact-line correction and is not a calibrated static Pc.',
                   measurements=measured)
    (study.root/'entry.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', choices=('bundle', 'entry'))
    parser.add_argument('--solver', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--scales', type=int, nargs='+')
    parser.add_argument('--dimensions', type=int, nargs='+', choices=(2, 3))
    parser.add_argument('--walls', nargs='+', choices=('legacy', 'renormalized'),
                        default=('legacy', 'renormalized'))
    args = parser.parse_args()
    study = Study(args.solver, args.root)
    if args.study == 'bundle':
        bundle_study(study, args.jobs, args.scales or (1, 2, 3), args.dimensions or (2, 3))
    else:
        entry_study(study, args.jobs, args.scales or (1, 2, 4), args.walls, args.dimensions or (2,))


if __name__ == '__main__':
    main()
