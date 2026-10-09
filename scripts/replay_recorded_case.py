#!/usr/bin/env python3
"""Replay an archived case after relocating paths; optionally compare fields.

The case directory must contain input.prm, geometry.raw and any restart seed.
The caller chooses a solver revision; historical kernels are documented in
research/drlab-publication/run-kernels. No build or download is performed.
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

from benchmarks import compare,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',required=True);p.add_argument('--case',required=True)
    p.add_argument('--destination',required=True);p.add_argument('--ranks',type=int,default=1)
    p.add_argument('--mpiexec',default='mpiexec');p.add_argument('--numproc-flag',default='-n')
    p.add_argument('--compare-fields',action='store_true')
    args=p.parse_args();source=Path(args.case).resolve();dest=Path(args.destination).resolve()
    dest.mkdir(parents=True,exist_ok=False)
    if any(c in str(dest) for c in ('"','\n','\r')):raise ValueError('The parameter-file path cannot contain quotes or newlines')
    shutil.copyfile(source/'geometry.raw',dest/'geometry.raw')
    text=(source/'input.prm').read_text()
    for key,target in [('gridFile',dest/'geometry.raw'),('outPath',dest/'field'),('logPath',dest)]:
        text,count=re.subn(r'^'+key+r'\s+"[^"]*";',lambda m:f'{key} "{target}";',text,flags=re.M)
        if count!=1:raise ValueError(f'The case must have exactly one {key}')
    start=int(re.search(r'^startTime\s+(\d+);',text,re.M)[1])
    if start:
        for extension in ('.rst','.rst.meta'):
            shutil.copyfile(source/f'field.{start:010d}{extension}',dest/f'field.{start:010d}{extension}')
    inp=dest/'input.prm';inp.write_text(text)
    solver=Path(args.solver).resolve();command=[str(solver),str(inp)]
    if args.ranks>1:command=[args.mpiexec,args.numproc_flag,str(args.ranks)]+command
    result=subprocess.run(command,cwd=dest,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (dest/'run.log').write_text(result.stdout)
    if result.returncode or 'Simulation finished!' not in result.stdout:raise RuntimeError(f'Case failed; inspect {dest/"run.log"}')
    fields=sorted(dest.glob('field.*.vti'))
    if not fields:raise RuntimeError('Replay produced no field output')
    summary=dict(solver_sha256=sha(solver),input_sha256=sha(inp),geometry_sha256=sha(dest/'geometry.raw'),command=command)
    if args.compare_fields:
        recorded=source/fields[-1].name
        if not recorded.exists():raise ValueError('Corresponding recorded final field is absent')
        summary['comparison']=compare(dict(name='recorded',output_files=[str(recorded)]),
                                      dict(name='replayed',output_files=[str(fields[-1])]))
    (dest/'replay.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':main()
