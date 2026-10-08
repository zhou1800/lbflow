#!/usr/bin/env python3
"""Sequential workstation timings with pinned ranks and repeated fixed inputs.

The in-solver MPI_Wtime numerical category includes exchange, boundaries,
macroscopic reconstruction, collision and finite checks. Field/diagnostic and
checkpoint I/O are reported separately. This is not a GPU or cluster benchmark.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np

from benchmarks import Study, sha
from rock_workflow import rock_case


def tagged_lines(path, tag):
    return [json.loads(line[len(tag)+1:]) for line in Path(path).read_text().splitlines()
            if line.startswith(tag+' ')]


def hardware(cpus):
    result = dict(platform=platform.platform(), python=platform.python_version(),
                  cpu_affinity=sorted(os.sched_getaffinity(0)), benchmark_cpus=cpus,
                  load_average=os.getloadavg(), meminfo=Path('/proc/meminfo').read_text(),
                  lscpu=subprocess.check_output(['lscpu'], text=True))
    result['cpu_frequency_metadata'] = {}
    for cpu in cpus:
        path = Path(f'/sys/devices/system/cpu/cpu{cpu}/cpufreq')
        result['cpu_frequency_metadata'][str(cpu)] = {
            key: (path/key).read_text().strip() for key in
            ('cpuinfo_max_freq', 'scaling_governor', 'scaling_driver') if (path/key).exists()}
    result['mpi_version'] = subprocess.check_output(['mpiexec', '--version'], text=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--cpus', type=int, nargs='+', required=True)
    parser.add_argument('--ranks', type=int, nargs='+', default=(1, 2, 4, 8))
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--steps', type=int, default=400)
    args = parser.parse_args()
    if max(args.ranks) > len(args.cpus) or args.repetitions < 3:
        raise ValueError('Use enough CPUs and at least three measured repetitions')
    before = hardware(args.cpus)
    wrapper = Path(__file__).with_name('pin_rank.py').resolve()
    study = Study(args.solver, args.root, launch_prefix=(sys.executable, str(wrapper)))
    study.env['LBFLOW_BENCH_CPUS'] = ','.join(map(str, args.cpus))
    study.env['OMPI_MCA_hwloc_base_binding_policy'] = 'none'
    measurements, warmups = [], []
    # Interleave rank counts by repetition to reduce monotonic thermal/load bias.
    for repetition in range(-1, args.repetitions):
        for ranks in args.ranks:
            case = rock_case(args.dataset, f'performance-r{ranks}-repeat{repetition}',
                model='color', tau=(1., 1.), force=2e-6, two_phase=True, steps=args.steps)
            # Fixed step count; reporting cannot stop the benchmark early.
            case.update(write=case['steps'], report_interval=case['steps'],
                        restart=case['steps'])
            load_before = list(os.getloadavg())
            result = study.run(case, ranks)
            timing = result['timing']
            if timing['completed_steps'] != case['steps'] or timing['ranks'] != ranks:
                raise AssertionError('Benchmark did not complete its fixed update count')
            folder = Path(result['output_files'][-1]).parent
            affinities = tagged_lines(folder/'run.log', 'LBFlowAffinity')
            memory = tagged_lines(folder/'run.log', 'LBFlowResources')
            if len(affinities) != ranks or len(memory) != ranks:
                raise AssertionError('Missing per-rank affinity/resource receipts')
            if sorted(a['cpus'][0] for a in affinities) != sorted(args.cpus[:ranks]):
                raise AssertionError('Unexpected rank affinity')
            updates = timing['fluid_nodes']*timing['completed_steps']
            record = dict(name=result['name'], ranks=ranks, repetition=repetition,
                timing=timing, process_seconds=result['seconds'], affinities=affinities,
                max_rss_kib_per_rank=memory, load_average_before=load_before,
                load_average_after=list(os.getloadavg()),
                fluid_MLUPS=updates/timing['numerical_seconds_max']/1e6,
                geometry_sha256=result['geometry_sha256'], solver_sha256=result['solver_sha256'])
            (warmups if repetition < 0 else measurements).append(record)
    # Full global fields must agree irrespective of rank count or repetition.
    from benchmarks import compare
    files = list((study.root/'results').glob('*.json'))
    runs = [json.loads(path.read_text()) for path in files]
    reference = next(r for r in runs if r['name'].startswith('performance-r1-repeat0-'))
    comparisons = [compare(reference, r) for r in runs if r['name'] != reference['name']]
    aggregates = []
    for ranks in args.ranks:
        records = [m for m in measurements if m['ranks'] == ranks]
        times = [r['timing']['numerical_seconds_max'] for r in records]
        throughput = [r['fluid_MLUPS'] for r in records]
        aggregates.append(dict(ranks=ranks, numerical_seconds_median=float(np.median(times)),
            numerical_seconds_range=[min(times), max(times)],
            fluid_MLUPS_median=float(np.median(throughput)), fluid_MLUPS_range=[min(throughput), max(throughput)],
            process_seconds_median=float(np.median([r['process_seconds'] for r in records])),
            fields_diagnostics_seconds_median=float(np.median([
                r['timing']['fields_diagnostics_seconds_max'] for r in records])),
            checkpoint_seconds_median=float(np.median([
                r['timing']['checkpoint_seconds_max'] for r in records])),
            max_rank_rss_mib=max(m['max_rss_kib'] for r in records for m in r['max_rss_kib_per_rank'])/1024))
    baseline = aggregates[0]['numerical_seconds_median']
    for item in aggregates:
        item['speedup_from_first_rank_count'] = baseline/item['numerical_seconds_median']
    summary = dict(hardware_before=before, hardware_after=hardware(args.cpus),
        protocol='One discarded full warm-up per rank count, then three or more sequential '
                 'repetitions. Rank counts interleaved; fixed geometry and update count; '
                 'one logical CPU per local rank. CPU frequency is not fixed.',
        timing_scope='In-solver category maxima across ranks; their sum need not equal the '
                     'maximum run wall time. Process time includes construction and launcher '
                     'overhead. Initial collision counted as initialization, not a completed update.',
        throughput_scope='Fluid voxel updates per maximum-rank numerical time, including MPI and '
                         'finite checks. Dense solid arrays also incur streaming cost.',
        dataset=json.loads((Path(args.dataset)/'manifest.json').read_text()),
        wrapper_sha256=sha(wrapper), warmups=warmups, measurements=measurements,
        aggregates=aggregates, field_comparisons=comparisons)
    (study.root/'performance.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
