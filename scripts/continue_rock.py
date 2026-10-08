#!/usr/bin/env python3
"""Continue a rock trajectory from a complete outgoing checkpoint.

A deliberately stopped parent run can supply a complete checkpoint; its exit
status is retained. Optional one-update comparisons check both global fields
and outgoing population bytes before changing MPI rank count.
"""
import argparse
import json
from pathlib import Path

from benchmarks import Study, compare, sha
from rock_workflow import measure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--seed-run', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--start', type=int, required=True)
    parser.add_argument('--end', type=int, required=True)
    parser.add_argument('--ranks', type=int, required=True)
    parser.add_argument('--check-ranks', type=int, nargs='+')
    args = parser.parse_args()
    seed = json.loads(Path(args.seed_run).read_text())
    if sha(args.solver) != seed['solver_sha256']:
        raise ValueError('Continue with the exact frozen executable used for the seed')
    if args.end <= args.start:
        raise ValueError('The end must follow the checkpoint time')
    folder = Path(seed['output_files'][0]).parent
    checkpoint = folder/f'field.{args.start:010d}'
    if not all(Path(str(checkpoint)+suffix).exists() for suffix in ('.rst', '.rst.meta')):
        raise ValueError('A complete population checkpoint and sidecar are required')
    geometry = folder/'geometry.raw'
    if sha(geometry) != seed['geometry_sha256']:
        raise ValueError('Seed geometry changed')
    # The production loader verifies dimensions, outgoing stage, size and
    # finite populations. Input physics is inherited without substitutions.
    case = dict(seed['parameters'], geometry_file=str(geometry), steps=args.end,
                name='rock-continued', timeout=7200)
    study = Study(args.solver, args.root)
    study.env.setdefault('OMPI_MCA_hwloc_base_binding_policy', 'none')
    checks, trials = [], []
    for ranks in args.check_ranks or ():
        trial = study.run(dict(case, name=f'rock-one-update-to{ranks}',
            steps=args.start+1, write=1, report_interval=1, restart=args.start+1),
            ranks, args.start, checkpoint)
        trials.append(trial)
    for trial in trials[1:]:
        comparison = compare(trials[0], trial)
        left = Path(trials[0]['output_files'][-1]).with_suffix('.rst')
        right = Path(trial['output_files'][-1]).with_suffix('.rst')
        comparison['outgoing_population_bytes_identical'] = sha(left) == sha(right)
        if not comparison['outgoing_population_bytes_identical']:
            raise AssertionError('Rank change altered outgoing populations')
        checks.append(comparison)
    result = study.run(case, args.ranks, args.start, checkpoint)
    merged = dict(result, parameters=dict(case), output_files=[
        path for path in seed['output_files'] if int(Path(path).name.split('.')[1]) <= args.start
    ]+result['output_files'])
    measurement = measure(merged)
    measurement['execution_segments'] = [dict(name=seed['name'], ranks=seed['ranks'],
        start=seed['start'], retained_until=args.start, exit=seed['exit'], finished=seed['finished']),
        dict(name=result['name'], ranks=result['ranks'], start=args.start, end=args.end,
             exit=result['exit'], finished=result['finished'])]
    summary = dict(measurement=measurement, rank_change_checks=checks,
        protocol='Complete outgoing checkpoint; exact executable and inherited physics; '
                 'closed-history audit combines retained initial snapshots and resumed output.',
        checkpoint_sha256=sha(Path(str(checkpoint)+'.rst')))
    (study.root/'rock-continuation.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
