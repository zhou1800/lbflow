#!/usr/bin/env python3
"""Continue recorded Laplace cases without relabelling early-time data as steady."""
import argparse,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from benchmarks import Study,capillary_cases,measure_capillary
p=argparse.ArgumentParser()
p.add_argument('--solver',required=True);p.add_argument('--initial-root',required=True)
p.add_argument('--root',required=True);p.add_argument('--jobs',type=int,default=4)
a=p.parse_args();study=Study(a.solver,a.root);initial=Path(a.initial_root)
def run(s):
    original=json.loads((initial/'results'/f'{s["name"]}-1rank.json').read_text())
    cut=s['steps'];folder=Path(original['output_files'][-1]).parent
    new=dict(s,name=s['name']+'-relaxed',steps=12000,write=2000,restart=12000,
             initial_output=original['output_files'][0])
    return study.run(new,1,cut,folder/f'field.{cut:010d}')
selected=[s for s in capillary_cases() if s['name'].startswith('laplace-d2-')]
with ThreadPoolExecutor(max_workers=a.jobs) as pool:results=list(pool.map(run,selected))
data=[measure_capillary(r) for r in results]
(study.root/'capillary.json').write_text(json.dumps({'measurements':data},indent=2))
for d in data:print(d['name'],d['relative_laplace_error'],d['late_pressure_relative_change'])
