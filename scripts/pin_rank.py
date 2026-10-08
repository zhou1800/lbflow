#!/usr/bin/env python3
"""Linux-only benchmark launcher: one specified logical CPU per local MPI rank."""
import json
import os
import resource
import subprocess
import sys


def main():
    if not sys.platform.startswith('linux') or not hasattr(os, 'sched_setaffinity'):
        raise RuntimeError('The explicit-affinity benchmark requires Linux')
    cpus = [int(value) for value in os.environ['LBFLOW_BENCH_CPUS'].split(',')]
    local_rank = int(os.environ.get('OMPI_COMM_WORLD_LOCAL_RANK', '0'))
    local_size = int(os.environ.get('OMPI_COMM_WORLD_LOCAL_SIZE', '1'))
    if local_size > len(cpus) or len(set(cpus)) != len(cpus):
        raise ValueError('Specify enough distinct CPUs for the MPI ranks')
    os.sched_setaffinity(0, {cpus[local_rank]})
    print('LBFlowAffinity '+json.dumps(dict(local_rank=local_rank,
        cpus=sorted(os.sched_getaffinity(0)))), flush=True)
    process = subprocess.run(sys.argv[1:])
    print('LBFlowResources '+json.dumps(dict(local_rank=local_rank,
        max_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)), flush=True)
    return process.returncode if process.returncode >= 0 else 128-process.returncode


if __name__ == '__main__':
    sys.exit(main())
