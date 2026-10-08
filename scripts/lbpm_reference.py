#!/usr/bin/env python3
"""Run unmodified LBPM as an independent single-phase pore-flow reference.

The external GPL-3.0 project is not bundled or relicensed as MIT. Verify and
build the pinned source separately, then supply its actual executable here.
No performance comparison is made: I/O, data layout and timers differ.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from benchmarks import sha


LBPM_REVISION = '6d686d354e5b8140841d3601e4c8c0e4e4b77e48'


def reference_case(binary, folder, pore, tau=1., force=2e-6, steps=6000):
    folder.mkdir(parents=True, exist_ok=False)
    raw = folder/'geometry.raw'
    pore.astype(np.uint8).tofile(raw)
    # Native rank-local input permits exact periodic ghost labels, rather than
    # inferring them from a pre-exchange decomposition file.
    np.pad(pore, 1, mode='wrap').astype(np.int8).tofile(folder/'ID.00000')
    nz, ny, nx = pore.shape
    # Unit voxel length: permeability is calculated independently from velocity
    # and the exact mask, rather than trusting LBPM's converted CSV columns.
    text = f'''MRT {{
 tau = {tau:.17g}
 F = {force:.17g}, 0, 0
 timestepMax = {steps}
 analysis_interval = {steps//6}
 tolerance = -1.0
 Restart = false
}}
Domain {{
 nproc = 1, 1, 1
 n = {nx}, {ny}, {nz}
 N = {nx}, {ny}, {nz}
 voxel_length = 1.0
 BC = 0
 ReadType = "8bit"
 ReadValues = 0, 1
 WriteValues = 0, 1
}}
Visualization {{
 write_silo = false
}}
'''
    input_file = folder/'input.db'
    input_file.write_text(text)
    start = time.monotonic()
    process = subprocess.run([str(binary), str(input_file)], cwd=folder,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1200)
    (folder/'run.log').write_text(process.stdout)
    process.check_returncode()
    rows = np.loadtxt(folder/'Permeability.csv', skiprows=1)
    if rows.shape[0] != 6 or rows[-1, 0] != steps or not np.isfinite(rows).all():
        raise AssertionError('LBPM did not produce the requested complete finite history')
    # Check that the external solver actually loaded the same geometry,
    # including the periodic padding used by its sparse layout.
    loaded = np.fromfile(folder/'ID.00000', dtype=np.int8).reshape((nz+2, ny+2, nx+2))
    if not np.array_equal(loaded > 0, np.pad(pore, 1, mode='wrap')):
        raise AssertionError('LBPM input mask or periodic padding differs')
    porosity = float(pore.mean())
    nu = (tau-.5)/3
    history = [dict(step=int(row[0]), pore_average_velocity=row[9:12].tolist(),
                    darcy_velocity=(porosity*row[9:12]).tolist(),
                    K_lattice=float(nu*porosity*row[9]/force)) for row in rows]
    result = dict(name=folder.name, tau=tau, force=force, porosity=porosity,
        size_xyz=[nx, ny, nz], source_revision=LBPM_REVISION,
        solver_sha256=sha(binary), geometry_sha256=sha(raw), input_sha256=sha(input_file),
        periodic_padded_mask_verified=True, exit=process.returncode,
        convergence_control='LBPM tolerance=-1 disables its nonnegative flow-change stopping condition.',
        process_seconds=time.monotonic()-start, history=history,
        csv_scope='Velocity columns printed to eight significant digits; K recalculated '
                  'from exact porosity, viscosity and force. Converted K columns unused.')
    (folder/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(name=folder.name, K_lattice=history[-1]['K_lattice'], exit=0)), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--source-checkout', required=True)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--refined-dataset', help='Optional subdivision of the same original crop')
    args = parser.parse_args()
    checkout = Path(args.source_checkout).resolve()
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    modified = subprocess.check_output(['git', '-C', str(checkout), 'diff', '--name-only'], text=True).strip()
    if revision != LBPM_REVISION or modified:
        raise ValueError('Supply the unmodified pinned LBPM v2026.04 source')
    root = Path(args.root).resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = root/'lbpm_permeability_simulator'
    shutil.copyfile(args.solver, binary)
    binary.chmod(0o755)
    measurements = []
    slit = np.ones((4, 18, 8), dtype=bool)
    slit[:, (0, -1), :] = False
    for tau in (.6, 1., 1.4):
        result = reference_case(binary, root/f'lbpm-slit-tau{tau}', slit, tau=tau, steps=18000)
        expected = (16**3/12+16/24)/18
        error = abs(result['history'][-1]['K_lattice']/expected-1)
        if error > 2e-6:
            raise AssertionError('Independent LBPM slit reference disagrees with Poiseuille flow')
        result.update(kind='slit', expected_K_midpoint=expected, relative_analytic_error=error)
        measurements.append(result)
    dataset = Path(args.dataset).resolve()
    manifest = json.loads((dataset/'manifest.json').read_text())
    path = dataset/'single-phase.raw'
    if sha(path) != manifest['files'][path.name]['sha256']:
        raise ValueError('Dataset checksum mismatch')
    nx, ny, nz = manifest['simulation_size_xyz']
    pore = np.fromfile(path, dtype=np.uint8).reshape(nz, ny, nx) == 75
    for tau, force in ((1., 1e-6), (1., 2e-6), (1., 4e-6), (.6, 2e-6), (1.4, 2e-6)):
        result = reference_case(binary, root/f'lbpm-rock-tau{tau}-F{force}', pore,
                                tau=tau, force=force, steps=18000 if tau == .6 else 6000)
        result.update(kind='rock', K_m2=result['history'][-1]['K_lattice']*
                      manifest['simulation_voxel_size_m']**2)
        measurements.append(result)
    summary = dict(software='LBPM v2026.04', source_revision=revision,
        source_url='https://github.com/OPM/LBPM/tree/'+revision,
        license='GPL-3.0; external reference, not distributed with LBFlow',
        comparison_scope='Independent single-phase permeability, not a matched multiphase '
                         'or runtime comparison. The three slit controls verify force, '
                         'viscosity and volume-average conventions first.',
        dataset=manifest, measurements=measurements)
    if args.refined_dataset:
        refined = Path(args.refined_dataset).resolve()
        metadata = json.loads((refined/'manifest.json').read_text())
        if any(metadata[key] != manifest[key] for key in ('crop_origin_xyz', 'crop_size_xyz')):
            raise ValueError('Refined reference must use the same original crop')
        path = refined/'single-phase.raw'
        if sha(path) != metadata['files'][path.name]['sha256']:
            raise ValueError('Refined image checksum mismatch')
        nx, ny, nz = metadata['simulation_size_xyz']
        pore = np.fromfile(path, dtype=np.uint8).reshape(nz, ny, nx) == 75
        scale = metadata['refinement']
        result = reference_case(binary, root/'lbpm-rock-refined', pore,
            force=2e-6/scale**3, steps=6000*scale**2)
        result.update(kind='rock-refined', K_m2=result['history'][-1]['K_lattice']*
                      metadata['simulation_voxel_size_m']**2)
        measurements.append(result)
        summary['refined_dataset'] = metadata
    (root/'lbpm-reference.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
