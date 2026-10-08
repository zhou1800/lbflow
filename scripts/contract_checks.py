#!/usr/bin/env python3
"""Failure contracts and legacy BB compatibility against the actual executable."""
import argparse,json,struct
from pathlib import Path
from benchmarks import Study,force_cases,channel_cases,compare
p=argparse.ArgumentParser()
p.add_argument('--solver',required=True);p.add_argument('--root',required=True)
a=p.parse_args();study=Study(a.solver,a.root)
base=force_cases()[-1]
seed=study.run(dict(base,name='checkpoint-seed'))
folder=Path(seed['output_files'][-1]).parent
source=folder/'field.0000000100'
records=[]
for name,change in [('missing-metadata','missing'),('short-populations','short'),('nan-population','nan')]:
    # Private copied checkpoints; valid source evidence remains unchanged.
    cp=study.root/(name+'-checkpoint')
    data=Path(str(source)+'.rst').read_bytes()
    if change=='short':data=data[:-8]
    if change=='nan':data=struct.pack('<d',float('nan'))+data[8:]
    Path(str(cp)+'.rst').write_bytes(data)
    Path(str(cp)+'.rst.meta').write_bytes(Path(str(source)+'.rst.meta').read_bytes())
    # Run copies the sidecar; removing it from the staged case is handled by
    # an invalid header so the same checker path is exercised without races.
    if change=='missing':Path(str(cp)+'.rst.meta').write_text('NO_VALID_STAGE\n')
    r=study.run(dict(base,name='guard-'+name,steps=200,expect_failure=True),1,100,cp)
    records.append(dict(name=r['name'],exit=r['exit']))
cp=study.root/'nan-population-checkpoint'
r=study.run(dict(base,name='guard-nan-population-parallel',steps=200,expect_failure=True),2,100,cp)
records.append(dict(name=r['name'],exit=r['exit']))
invalid=dict(base,name='guard-invalid-phase',indices=['75 "GV" "0.000001" "1" "2"'],expect_failure=True)
records.append({'name':invalid['name'],'exit':study.run(invalid)['exit']})
s=next(s for s in channel_cases() if s['model']=='scsp' and s['eq']=='MRT' and s['tau']==(.6,) and s['height']==8)
normal=study.run(dict(s,name='solid-SD-reference'))
bb=study.run(dict(s,name='solid-BB-reference',indices=[s['indices'][0],'255 "BB" "N/A" "N/A" "N/A"']))
result=dict(failure_checks=records,bb_compatibility=compare(normal,bb))
(study.root/'contracts.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
