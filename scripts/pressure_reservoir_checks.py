#!/usr/bin/env python3
"""Pressure CSF hydraulic references, component budgets and guarded scope."""
import argparse
import json
from pathlib import Path
import uuid

import numpy as np

from benchmarks import Study, compare, fields
from pressure_reservoirs import budget, channel, channel_measure


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True);parser.add_argument('--root',required=True)
    parser.add_argument('--mpiexec',default='mpiexec');parser.add_argument('--numproc-flag',default='-n')
    args=parser.parse_args();study=Study(args.solver,Path(args.root)/('pressure-reservoirs-'+uuid.uuid4().hex[:10]),args.mpiexec,args.numproc_flag)
    hydraulic=[]
    for spectrum in ('legacy','liu2015'):
        for phase,sign in ((0,-1),(1,1)):
            cfg=channel(8,spectrum,phase,sign,tau=1.)
            r=study.run(cfg);row=channel_measure(r)
            if row['relative_profile_error']>3e-5:raise AssertionError('Pressure-driven midpoint profile differs from the independent MRT reference')
            hydraulic.append(row)
    cfg=channel(16,'liu2015',1,1,tau=1.)
    cfg.update(name='pressure-front-short',tau=(.6,1.4),viscosity_interpolation='harmonic',
        sigma=.0005,steps=256,write=128,restart=128,
        label=lambda x,y,z:255 if y in (0,17) else 50 if x==0 else 51 if x==63 else 0 if x<32 else 100,
        indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"',
                 '50 "PB" "1.000015" "1" "0"','51 "PB" "0.999985" "3" "1"','255 "SD" "60" "N/A" "N/A"'])
    ref=study.run(cfg);parallel=study.run(dict(cfg,name=cfg['name']+'-parallel'),2)
    comparisons=[compare(ref,parallel)];receipts=[ref,parallel]
    checkpoint=Path(ref['output_files'][-1]).parent/'field.0000000128'
    for ranks in (1,2):
        resumed=study.run(dict(cfg,name=cfg['name']+'-resumed-'+str(ranks)),ranks,128,checkpoint)
        comparisons.append(compare(ref,resumed));receipts.append(resumed)
    native=Path(ref['output_files'][-1]).parent/'field.0000000256.rst'
    for r in receipts[1:]:
        other=Path(r['output_files'][-1]).parent/'field.0000000256.rst'
        if native.read_bytes()!=other.read_bytes():raise AssertionError('Pressure MPI/restart outgoing populations differ')
    budgets=[budget(r) for r in receipts]
    rotated=dict(cfg,name='pressure-front-rotated',size=(18,64,1),
        label=lambda x,y,z:255 if x in (0,17) else 50 if y==0 else 51 if y==63 else 0 if y<32 else 100,
        indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"',
                 '50 "PB" "1.000015" "2" "0"','51 "PB" "0.999985" "4" "1"','255 "SD" "60" "N/A" "N/A"'])
    rr=study.run(rotated,ranks=2);left=fields(ref);right=fields(rr);rotation={}
    for name in left:
        array=right[name].transpose(0,2,1,3)
        if name=='v':array=array[...,[1,0]]
        rotation[name]=float(abs(left[name]-array).max())
    if max(rotation.values())>3e-12:raise AssertionError('Rotated pressure/halo geometry changes the solution')
    budgets.append(budget(rr))
    # Each component must exchange a nonzero amount with its own reservoir.
    if min(abs(np.array(budgets[0]['history'][-1]['integrated_net_mass'])))<1e-6:
        raise AssertionError('Two-phase component reservoir exchange was not exercised')
    base=channel(8,'liu2015',1,1,tau=1.);base.update(steps=2,write=1)
    left='50 "PB" "1.000015" "1" "1"';right='51 "PB" "0.999985" "3" "1"';wall='255 "SD" "60" "N/A" "N/A"'
    invalid=[('unsupported-wall',dict(wall_gradient='halfway')),
        ('unsupported-phase',dict(phase_boundary='bounce')),
        ('velocity-reservoir',dict(indices=['75 "VD" "N/A" "N/A" "1"',left,'51 "VB" "0.00001" "3" "1"',wall])),
        ('missing-face',dict(indices=['75 "VD" "N/A" "N/A" "1"',left,wall],label=lambda x,y,z:255 if y in (0,9) else 50 if x==0 else 75)),
        ('outward-face',dict(indices=['75 "VD" "N/A" "N/A" "1"','50 "PB" "1.000015" "3" "1"',right,wall])),
        ('interior-pressure',dict(label=lambda x,y,z:255 if y in (0,9) else 50 if x in (0,8) else 51 if x==31 else 75)),
        ('incomplete-face',dict(label=lambda x,y,z:255 if y in (0,9) else 50 if x==0 and y!=4 else 51 if x==31 else 75)),
        ('mixed-collar',dict(indices=['75 "VD" "N/A" "N/A" "1"','0 "VD" "N/A" "N/A" "0"',left,right,wall],
            label=lambda x,y,z:255 if y in (0,9) else 50 if x==0 else 51 if x==31 else 0 if x==2 else 75)),
        ('body-force',dict(indices=['75 "GV" "0.000001" "1" "1"',left,right,wall])),
        ('wrong-axis',dict(indices=['75 "VD" "N/A" "N/A" "1"',left,'51 "PB" "0.999985" "4" "1"',wall]))]
    guards=[]
    for name,change in invalid:
        r=study.run(dict(base,name=name,expect_failure=True,**change));guards.append(dict(name=name,exit=r['exit']))
    # Check the outgoing collar path after loading a valid native stage.
    poisoned=study.root/'poisoned-collar';f=np.frombuffer(Path(str(checkpoint)+'.rst').read_bytes(),dtype='=f8').copy().reshape(2,9,18,64)
    f[0,0,8,1]-=.001;f[1,0,8,1]+=.001
    Path(str(poisoned)+'.rst').write_bytes(f.tobytes());Path(str(poisoned)+'.rst.meta').write_bytes(Path(str(checkpoint)+'.rst.meta').read_bytes())
    r=study.run(dict(cfg,name='outgoing-mixed-collar',expect_failure=True),1,128,poisoned)
    guards.append(dict(name=r['name'],exit=r['exit']))
    # A sub-tolerance other-phase fraction is still invalid if its actual
    # ghost gradient crosses the force cutoff at the first strip.
    microscopic=study.root/'microscopic-collar';f=np.frombuffer(Path(str(checkpoint)+'.rst').read_bytes(),dtype='=f8').copy().reshape(2,9,18,64)
    total=f[:,:,1:-1,2].sum(axis=0)
    f[0,:,1:-1,2]=total*(1-9e-10);f[1,:,1:-1,2]=total*9e-10
    Path(str(microscopic)+'.rst').write_bytes(f.tobytes());Path(str(microscopic)+'.rst.meta').write_bytes(Path(str(checkpoint)+'.rst.meta').read_bytes())
    controls=dict(cfg,name='incoming-force-collar',expect_failure=True,
                  indices=cfg['indices'][:-1]+['255 "SD" "1" "N/A" "N/A"'])
    r=study.run(controls,1,128,microscopic)
    log=(Path(r['command'][-1]).parent/'run.log').read_text()
    if 'force-free first interior strip' not in log:raise AssertionError('Actual collar-gradient guard was not exercised')
    guards.append(dict(name=r['name'],exit=r['exit']))
    summary=dict(hydraulic=hydraulic,comparisons=comparisons,bitwise_population_comparisons=3,budgets=budgets,
        guards=guards,rotation_errors=rotation,solver_sha256=study.solver_sha,scope='Flat pure-collar D2Q9 pressure CSF; quantitative contact-line/throat accuracy is a separate physical experiment.')
    (study.root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(dict(hydraulic_cases=len(hydraulic),guards=len(guards),comparisons=len(comparisons),root=str(study.root)),indent=2))


if __name__=='__main__':main()
