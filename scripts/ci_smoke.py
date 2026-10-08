#!/usr/bin/env python3
"""Small executable regressions for forcing, MPI, restart and failure contracts."""
import argparse
import json
from pathlib import Path
import struct
import uuid

import numpy as np
from benchmarks import Study, compare, force_cases, geometry_fields, spec


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--numproc-flag', default='-n')
    args = parser.parse_args()
    root = Path(args.root) / ('contracts-' + uuid.uuid4().hex[:10])
    study = Study(args.solver, root, args.mpiexec, args.numproc_flag)
    force_results = []
    for case in force_cases():
        result = study.run(case)
        timing = result['timing']
        if not timing or timing['completed_steps'] != case['steps']:
            raise AssertionError('Timing summary did not count the completed updates')
        if timing['fluid_nodes'] != np.prod(case['size']) or timing['ranks'] != 1:
            raise AssertionError('Timing summary has an incorrect lattice/rank count')
        durations = [value for key, value in timing.items() if key.endswith('_seconds_max')]
        if not all(np.isfinite(value) and value >= 0 for value in durations):
            raise AssertionError('Invalid timing duration')
        if timing['numerical_seconds_max'] > timing['run_seconds_max']:
            raise AssertionError('Numerical duration exceeds the enclosing run duration')
        velocity = geometry_fields(result)[3]
        error = float(np.max(np.abs(velocity[..., 0] - 1e-4)))
        transverse = float(np.max(np.abs(velocity[..., 1:])))
        if max(error, transverse) > 1e-12:
            raise AssertionError(f'{case["name"]}: incorrect force integration')
        force_results.append({'name': result['name'], 'absolute_error': error})

    drop = spec('drop-continuous', nx=12, ny=12, nz=12, steps=40,
                write=20, restart=20,
                label=lambda x, y, z: 100 if (x-5.5)**2+(y-5.5)**2+(z-5.5)**2 <= 9 else 0,
                indices=['0 "VD" "N/A" "N/A" "0"', '100 "VD" "N/A" "N/A" "1"'])
    continuous = study.run(drop)
    parallel = study.run(dict(drop, name='drop-parallel'), ranks=2)
    comparisons = [compare(continuous, parallel)]
    checkpoint = Path(continuous['output_files'][-1]).parent / 'field.0000000020'
    for ranks in (1, 2):
        resumed = study.run(dict(drop, name=f'drop-resumed-to{ranks}'), ranks, 20, checkpoint)
        if resumed['timing']['completed_steps'] != 20 or resumed['timing']['ranks'] != ranks:
            raise AssertionError('Incorrect restart/rank timing accounting')
        comparisons.append(compare(continuous, resumed))
    if not all(item['bitwise_fields'] for item in comparisons):
        raise AssertionError('MPI/restart fields are not bitwise identical')

    base = force_cases()[-1]
    failures = []
    for name, overrides in [('invalid-tau', {'rexlationTime': '.5 .5'}),
                            ('invalid-lattice-speed', {'latticeSpeed': 2})]:
        result = study.run(dict(base, name=name, overrides=overrides, expect_failure=True))
        failures.append({'name': name, 'exit': result['exit']})
    for name in ('missing-stage', 'short-populations', 'nonfinite-population'):
        copied = root / name
        populations = Path(str(checkpoint)+'.rst').read_bytes()
        metadata = Path(str(checkpoint)+'.rst.meta').read_text()
        if name == 'missing-stage':
            metadata = 'NO_VALID_STAGE\n'
        elif name == 'short-populations':
            populations = populations[:-8]
        else:
            populations = struct.pack('<d', float('nan')) + populations[8:]
        Path(str(copied)+'.rst').write_bytes(populations)
        Path(str(copied)+'.rst.meta').write_text(metadata)
        result = study.run(dict(drop, name=name, expect_failure=True),
                           2 if name == 'nonfinite-population' else 1, 20, copied)
        failures.append({'name': name, 'exit': result['exit']})
    (root/'summary.json').write_text(json.dumps(
        {'forcing': force_results, 'comparisons': comparisons, 'failures': failures}, indent=2))
    print(f'Passed {len(force_results)} forcing cases, {len(comparisons)} field comparisons '
          f'and {len(failures)} failure contracts. Evidence: {root}')


if __name__ == '__main__':
    main()
