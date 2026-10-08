#!/usr/bin/env python3
"""Separate diffuse viscosity bias from lattice error in equal-density slit flow.

Continuum references follow d(mu du/dy)/dy = -F with no-slip halfway walls.
The planar-width prediction uses a tanh profile, with no flow-data fitting.
This is a controlled constitutive study, not a new viscosity model.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from benchmarks import Study, geometry_fields, sha, spec
from porous_media import closed_history


def sharp_reference(y, height, viscosities, force):
    interface = height / 2
    mu0, mu1 = viscosities
    moments = [(interface**(n+1)/mu0 +
                (height**(n+1)-interface**(n+1))/mu1)/(n+1) for n in range(3)]
    constant = force * moments[1] / moments[0]
    first = (constant*y-force*y*y/2)/mu0
    second = ((constant*interface-force*interface**2/2)/mu0 +
              (constant*(y-interface)-force*(y*y-interface**2)/2)/mu1)
    return np.where(y <= interface, first, second), force*(
        moments[2]-moments[1]**2/moments[0])


def diffuse_reference(y, grid, viscosity, force):
    inverse = 1/viscosity
    moments = [float(np.trapezoid(grid**n*inverse, grid)) for n in range(3)]
    constant = force*moments[1]/moments[0]
    derivative = (constant-force*grid)*inverse
    integral = np.concatenate(([0.], np.cumsum(
        (derivative[1:]+derivative[:-1])*np.diff(grid)/2)))
    return dict(profile=np.interp(y, grid, integral),
                volume_flux=force*(moments[2]-moments[1]**2/moments[0]),
                moments=moments, stress_intercept=constant,
                endpoint_residual=float(integral[-1]))


def viscosity(phi, viscosities, closure):
    mu0, mu1 = viscosities
    if closure == 'arithmetic':
        return phi*mu0+(1-phi)*mu1
    return 1/(phi/mu0+(1-phi)/mu1)


def thin_interface_coefficients(ratio):
    contrast = (1-ratio)/(1+ratio)
    denominator = 4-3*contrast**2
    return dict(arithmetic_linear=3*contrast**3*math.log(ratio)/denominator,
                harmonic_quadratic=2*math.pi**2*contrast**2/denominator,
                small_parameter='xi/H for a centered tanh interface; continuous volume flux')


def check_continuum():
    height, force, mu = 16., 1e-6, 1/6
    grid = np.linspace(0, height, 16385)
    y = np.arange(16)+.5
    reference = diffuse_reference(y, grid, np.full_like(grid, mu), force)
    expected = force*y*(height-y)/(2*mu)
    profile_error = float(np.max(np.abs(reference['profile']-expected)))
    flux_error = abs(reference['volume_flux']/(force*height**3/(12*mu))-1)
    # Independent quadrature of the infinite-domain tanh excess compliance.
    xi, mu0, mu1 = 2., 1/30, .3
    position = np.linspace(-30*xi, 30*xi, 131073)
    phi = (1-np.tanh(position/xi))/2
    sharp = np.where(position < 0, 1/mu0, 1/mu1)
    # Assign the symmetric value at the jump for trapezoidal integration.
    sharp[len(sharp)//2] = (1/mu0+1/mu1)/2
    measured = float(np.trapezoid(1/viscosity(phi, (mu0, mu1), 'arithmetic')-sharp,
                                  position))
    exact = xi/2*(1/mu1-1/mu0)*math.log(mu1/mu0)
    harmonic = float(np.trapezoid(1/viscosity(phi, (mu0, mu1), 'harmonic')-sharp,
                                  position))
    if profile_error > 1e-12 or flux_error > 1e-7 or abs(measured-exact) > 1e-8 or abs(harmonic) > 1e-8:
        raise AssertionError('Independent continuum controls failed')
    # Check the analytic thin-interface orders independently by quadrature.
    fine = np.linspace(0, 1, 262145)
    c, d = (1/mu0+1/mu1)/2, (1/mu1-1/mu0)/2
    sharp_flux = c/12-d*d/(16*c)
    coefficients = thin_interface_coefficients(mu1/mu0)
    limits = []
    for epsilon in (.01, .005, .0025):
        fraction = (1-np.tanh((fine-.5)/epsilon))/2
        values = {}
        for closure, order in (('arithmetic', 1), ('harmonic', 2)):
            flow = diffuse_reference(np.array([.5]), fine, viscosity(fraction, (mu0, mu1),
                                                                  closure), 1.)['volume_flux']
            values[closure] = float((flow/sharp_flux-1)/epsilon**order)
        limits.append(dict(xi_over_height=epsilon, coefficients=values))
    for closure, label in (('arithmetic', 'arithmetic_linear'), ('harmonic', 'harmonic_quadratic')):
        errors = [abs(item['coefficients'][closure]-coefficients[label]) for item in limits]
        if not errors[2] < errors[1] < errors[0]:
            raise AssertionError('Continuum thin-interface limits do not approach their derivation')
    return dict(constant_viscosity_profile_max_error=profile_error,
                constant_viscosity_flux_relative_error=flux_error,
                tanh_compliance_quadrature=measured, tanh_compliance_formula=exact,
                harmonic_excess_compliance=harmonic,
                analytic_thin_interface_coefficients=coefficients,
                continuum_quadrature_limits=limits)


def planar(study, betas):
    records = []
    for beta in betas:
        case = spec(f'transport-planar-beta{beta}', nx=128, ny=4, nz=1,
                    dimensions=2, tau=(.6, 1.4), beta=beta, steps=8000, write=4000,
                    label=lambda x, y, z: 75 if x < 64 else 76,
                    indices=['75 "VD" "N/A" "N/A" "0"',
                             '76 "VD" "N/A" "N/A" "1"'])
        result = study.run(case)
        _, rho, chi, velocity = geometry_fields(result)
        line = chi[0, 0, 32:96]
        positions = []
        for level in (.9, -.9):
            intersections = np.flatnonzero((line[:-1] >= level) & (line[1:] < level))
            if len(intersections) != 1:
                raise AssertionError('Planar reference lacks a unique interface')
            i = intersections[0]
            positions.append(float(i+(level-line[i])/(line[i+1]-line[i])))
        previous = geometry_fields(result, 4000)[2]
        width = positions[1]-positions[0]
        records.append(dict(beta=beta, width_09=width, xi_tanh=width/(2*np.arctanh(.9)),
                            max_speed=float(np.linalg.norm(velocity, axis=-1).max()),
                            late_color_max_change=float(np.max(np.abs(chi-previous))),
                            density_range=[float(rho.min()), float(rho.max())],
                            solver_sha256=result['solver_sha256'],
                            input_sha256=result['input_sha256'],
                            geometry_sha256=result['geometry_sha256'],
                            **closed_history(result)))
    return records


def cases(betas, heights):
    for beta in betas:
        for height in heights:
            scale = height/16
            force = 1e-6/scale**3
            steps = math.ceil(20*height**2/((.6-.5)/3*math.pi**2)/1000)*1000
            for closure in ('arithmetic', 'harmonic'):
                yield spec(f'transport-H{height}-beta{beta}-{closure}', nx=4, ny=height+2,
                           nz=1, dimensions=2, tau=(.6, 1.4), steps=steps,
                           write=steps//2, beta=beta, sigma=.005/scale,
                           force=force, height=height, viscosity_interpolation=closure,
                           timeout=1800,
                           label=lambda x, y, z, h=height: 255 if y in (0, h+1) else (
                               75 if y <= h//2 else 76),
                           indices=[f'75 "GV" "{force:.17g}" "1" "0"',
                                    f'76 "GV" "{force:.17g}" "1" "1"',
                                    '255 "SD" "90" "N/A" "N/A"'])


def measure(result, calibration):
    p = result['parameters']
    arrays, rho, chi, velocity = geometry_fields(result)
    height, closure = p['height'], p['viscosity_interpolation']
    viscosities = [(tau-.5)/3 for tau in p['tau']]
    y = np.arange(height)+.5
    phi = (arrays['r0'][..., 0]/np.where(rho > 0, rho, 1)).mean(axis=(0, 2))[1:-1]
    observed = velocity[..., 0].mean(axis=(0, 2))[1:-1]
    grid = np.linspace(0, height, height*256+1)
    sharp_profile, sharp_flux = sharp_reference(y, height, viscosities, p['force'])
    actual = diffuse_reference(y, grid, viscosity(np.interp(grid, y, phi), viscosities,
                                                  closure), p['force'])
    xi = calibration['xi_tanh']
    predicted_phi = (1-np.tanh((grid-height/2)/xi))/2
    predicted = diffuse_reference(y, grid, viscosity(predicted_phi, viscosities, closure),
                                   p['force'])
    previous = geometry_fields(result, p['steps']-p['write'])
    error = lambda profile: float(np.linalg.norm(observed-profile)/np.linalg.norm(profile))
    change = float(np.linalg.norm(velocity-previous[3])/np.linalg.norm(velocity))
    mass = closed_history(result)['max_relative_mass_drift']
    # Midpoint sampling is shared by the velocity-profile comparisons.
    sampled_sharp = float(sharp_profile.sum())
    record = dict(name=result['name'], height=height, beta=p['beta'],
                  closure=closure, tau=p['tau'], force=p['force'], sigma=p['sigma'],
                  steps=p['steps'], width_09=calibration['width_09'],
                  width_over_height=calibration['width_09']/height,
                  profile=observed.tolist(), phase0_profile=phi.tolist(),
                  sharp_profile=sharp_profile.tolist(),
                  measured_diffuse_profile=actual['profile'].tolist(),
                  planar_tanh_predicted_profile=predicted['profile'].tolist(),
                  sharp_profile_relative_error=error(sharp_profile),
                  measured_diffuse_profile_relative_error=error(actual['profile']),
                  tanh_prediction_profile_relative_error=error(predicted['profile']),
                  measured_flux_shift=float(observed.sum()/sampled_sharp-1),
                  measured_diffuse_flux_shift=float(actual['profile'].sum()/sampled_sharp-1),
                  tanh_predicted_flux_shift=float(predicted['profile'].sum()/sampled_sharp-1),
                  sharp_continuous_volume_flux=sharp_flux,
                  diffuse_continuous_volume_flux=actual['volume_flux'],
                  planar_tanh_continuous_volume_flux=predicted['volume_flux'],
                  diffuse_moments=actual['moments'],
                  late_velocity_relative_change=change,
                  late_color_max_change=float(np.max(np.abs(chi-previous[2]))),
                  max_speed=float(np.linalg.norm(velocity, axis=-1).max()),
                  density_range=[float(rho[rho > 0].min()), float(rho.max())],
                  max_relative_mass_drift=mass, solver_sha256=result['solver_sha256'],
                  input_sha256=result['input_sha256'], geometry_sha256=result['geometry_sha256'])
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--betas', nargs='+', type=float, default=[.35, .7, .9])
    parser.add_argument('--heights', nargs='+', type=int, default=[16, 32, 64])
    args = parser.parse_args()
    if any(h < 8 or h % 2 for h in args.heights) or any(b <= 0 or b >= 1 for b in args.betas):
        parser.error('Use even heights >= 8 and 0 < beta < 1')
    checks = check_continuum()
    study = Study(args.solver, args.root)
    calibration = planar(study, args.betas)
    by_beta = {m['beta']: m for m in calibration}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(study.run, cases(args.betas, args.heights)))
    records = [measure(result, by_beta[result['parameters']['beta']]) for result in results]
    root = Path(__file__).resolve().parents[1]
    summary = dict(protocol='Equal-density D2Q9 MRT layered slit; viscosity ratio 9:1; '
                   'no slip; separately measured planar interface; no flow-data fitting.',
                   scaling='Fixed relaxation times, F ~ H^-3 and sigma ~ H^-1. Fixed beta '
                   'shrinks the physical interface with refinement. This is a coupled '
                   'sharp-interface/grid sequence, not a fixed-interface convergence study.',
                   scope='Controlled viscous transport only. No capillary-entry, trapping, '
                   'universal confinement threshold or novel viscosity closure is demonstrated.',
                   source_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root,
                                                           text=True).strip(),
                   analysis_script_sha256=sha(__file__), checks=checks,
                   planar_calibration=calibration, measurements=records)
    (study.root/'interfacial-transport.json').write_text(json.dumps(summary, indent=2)+'\n')
    for m in records:
        print(json.dumps({k: m[k] for k in ('name', 'width_over_height', 'measured_flux_shift',
                         'tanh_predicted_flux_shift', 'measured_diffuse_profile_relative_error',
                         'late_velocity_relative_change')}))


if __name__ == '__main__':
    main()
