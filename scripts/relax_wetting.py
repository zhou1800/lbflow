#!/usr/bin/env python3
"""Run one bounded batch of long wetting cases from valid population checkpoints."""
import argparse,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from benchmarks import Study
from wetting import cases,measure
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--solver',required=True);p.add_argument('--initial-root',required=True)
p.add_argument('--root',required=True);p.add_argument('--radius',type=int,required=True)
p.add_argument('--angle',type=int);p.add_argument('--method',choices=('legacy','renormalized'))
p.add_argument('--ranks',type=int,default=4);p.add_argument('--jobs',type=int,default=2)
p.add_argument('--end-time',type=int,help='Override the radius-based final step, e.g. 48000 for R12/theta60')
a=p.parse_args();study=Study(a.solver,a.root);initial=Path(a.initial_root)
runs=[s for s in cases() if s['radius']==a.radius and (a.angle is None or s['angle']==a.angle)
      and (a.method is None or s['wall_gradient']==a.method)
      and (s['angle']!=90 or a.radius==48)]
def run(s):
    name=s['name'];receipt=initial/'results'/f'{name}-1rank.json'
    end=a.end_time if a.end_time is not None else 24000*(s['radius']//12)**2
    if end<=s['steps']:raise ValueError('Continuation end time must exceed the initial run duration')
    new=dict(s,name=name+'-steady',steps=end,write=end//4,restart=end,timeout=1500)
    if receipt.exists():
        original=json.loads(receipt.read_text())
        if original['exit']!=0 or not original['finished']:raise RuntimeError('Cannot resume failed run')
        cut=s['steps'];folder=Path(original['output_files'][-1]).parent
        new['initial_output']=original['output_files'][0]
        result=study.run(new,a.ranks,cut,folder/f'field.{cut:010d}')
    else:
        # Interrupted macro outputs are insufficient to reconstruct populations.
        result=study.run(new,a.ranks)
    result['measurement']=measure(result)
    (study.root/'results'/f'{result["name"]}.measurement.json').write_text(json.dumps(result['measurement'],indent=2))
    print(result['name'],result['measurement']['principal']['angle_degrees'],
          result['measurement']['late_angle_change'],flush=True)
    return result['measurement']
with ThreadPoolExecutor(max_workers=a.jobs) as pool:data=list(pool.map(run,runs))
(study.root/f'wetting-R{a.radius}-theta{a.angle}-{a.method}.json').write_text(json.dumps({'measurements':data},indent=2))
