#!/usr/bin/env python3
"""Pure hydraulics, phase transport and replay contracts for prescribed slip."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import uuid

import numpy as np

from benchmarks import Study, compare, geometry_fields, fields
from periodic_contact import config, seed
from phase_recurrence import VELOCITIES


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True);parser.add_argument('--root',required=True)
    parser.add_argument('--mpiexec',default='mpiexec');parser.add_argument('--numproc-flag',default='-n')
    args=parser.parse_args();study=Study(args.solver,Path(args.root)/('navier-slip-'+uuid.uuid4().hex[:10]),args.mpiexec,args.numproc_flag)
    configs=[]
    for h in (8,16,32):
        for length in (.5,1.,2.):
            for sign in (-1,1):
                cfg=config(h,sign*1e-7,end=5000.)
                cfg.update(name=f'pure-H{h}-slip{length:g}-sign{sign}',overrides={'wallSlipLength':length/(16/h)},
                    label=lambda x,y,z,hh=h:255 if y in (0,hh+1) else 100)
                configs.append(cfg)
    with ThreadPoolExecutor(max_workers=3) as pool:runs=list(pool.map(study.run,configs))
    hydraulic=[]
    for r in runs:
        cfg=r['parameters'];h=cfg['height'];d=cfg['dx'];_,rho,chi,v=geometry_fields(r)
        slip=cfg['overrides']['wallSlipLength'];y=np.arange(h)+.5;f=cfg['physical_force']*d**3
        reference=f*3*(y*(h-y)+slip*h)
        actual=v[0,1:-1,:,0].mean(axis=1);error=float(abs(actual-reference).max()/abs(reference).max())
        if error>3e-6:raise AssertionError('Prescribed pure slip does not match the independent hydraulic formula')
        hydraulic.append(dict(name=r['name'],relative_profile_error=error,slip_lattice=slip))
    cfg=config(16,1e-7,end=128.);cfg.update(name='front-continuous',steps=128,write=64,restart=64,overrides={'wallSlipLength':1.})
    path=study.root/'planar-seed';beta,phi=seed(path,16);cfg['beta']=beta
    reference=study.run(cfg,start=1,checkpoint=path)
    parallel=study.run(dict(cfg,name='front-parallel'),ranks=2,start=1,checkpoint=path)
    expected_masses=np.array([phi.sum(),phi.size-phi.sum()]);mass_checks=[]
    def check_mass(receipt):
        for filename in receipt['output_files']:
            step=int(Path(filename).name.split('.')[1]);_,rho,chi,v=geometry_fields(receipt,step)
            actual=np.array([((1+chi)*rho/2).sum(),((1-chi)*rho/2).sum()])
            error=float(abs(actual/expected_masses-1).max())
            if error>3e-11:raise AssertionError('Prescribed slip changes a closed component mass')
            mass_checks.append(dict(name=receipt['name'],step=step,max_relative_component_mass_error=error))
    check_mass(reference);check_mass(parallel)
    comparisons=[compare(reference,parallel)]
    folder=Path(reference['output_files'][-1]).parent
    for ranks in (1,2):
        r=study.run(dict(cfg,name=f'front-resumed-{ranks}'),ranks=ranks,start=64,checkpoint=folder/'field.0000000064')
        check_mass(r)
        comparisons.append(compare(reference,r))
        if (Path(r['output_files'][-1]).parent/'field.0000000128.rst').read_bytes()!=(folder/'field.0000000128.rst').read_bytes():
            raise AssertionError('Slip restart changes outgoing populations')
    neutral=study.run(dict(cfg,name='neutral-fixed',steps=128,restart=128,
        indices=['0 "GV" "0" "1" "0"','100 "GV" "0" "1" "1"','255 "SD" "90" "N/A" "N/A"']),start=1,checkpoint=path)
    _,rho,chi,v=geometry_fields(neutral);color_error=float(abs((chi[0,1:-1]+1)/2-phi).max())
    if color_error>3e-11 or np.linalg.norm(v,axis=-1).max()>2e-14:raise AssertionError('Slip alters neutral rest equilibrium')
    pop=np.fromfile(str(path)+'.rst',dtype='=f8').reshape(2,9,18,64)
    swap=[int(np.where(np.all(VELOCITIES==e[::-1],axis=1))[0][0]) for e in VELOCITIES]
    rotated_path=study.root/'rotated-seed';rotated=pop[:,swap].transpose(0,1,3,2)
    Path(str(rotated_path)+'.rst').write_bytes(rotated.tobytes())
    Path(str(rotated_path)+'.rst.meta').write_text('DRLAB_OUTGOING_V1\n1 2 2 9\n64 18 \nNEXT_STREAM 2\nINITIAL_COLLISION_INCLUDED 2\nCONVERGENCE_HISTORY_RESET_ON_LOAD\n')
    rotated_cfg=dict(cfg,name='front-rotated',size=(18,64,1),
        label=lambda x,y,z:255 if x in (0,17) else 0 if 16<=y<48 else 100,
        indices=['0 "GV" "0.0000001" "2" "0"','100 "GV" "0.0000001" "2" "1"','255 "SD" "90" "N/A" "N/A"'])
    rr=study.run(rotated_cfg,ranks=2,start=1,checkpoint=rotated_path)
    left,right=fields(reference),fields(rr);rotation={}
    for name in left:
        data=right[name].transpose(0,2,1,3)
        if name=='v':data=data[...,[1,0]]
        rotation[name]=float(abs(left[name]-data).max())
    if max(rotation.values())>3e-12:raise AssertionError('Rotated slip/phase/halo transport changes the fields')
    guards=[]
    for name,change in [('negative',{'overrides':{'wallSlipLength':-1}}),
        ('nonneutral',{'indices':['0 "GV" "0" "1" "0"','100 "GV" "0" "1" "1"','255 "SD" "60" "N/A" "N/A"']}),
        ('unequal-tau',{'tau':(.8,1.)}),('srt',{'eq':'SRT'}),('cosine',{'recoloring_model':'cosine'}),
        ('bad-spectrum',{'mrt_spectrum':'liu2015'}),('unrecognized-number',{'overrides':{'wallSlipLength':'"bad"'}})]:
        r=study.run(dict(cfg,**change,name='invalid-'+name,expect_failure=True,steps=2,write=1,restart=2),start=1,checkpoint=path)
        guards.append(dict(name=name,exit=r['exit']))
    result=dict(hydraulic=hydraulic,comparisons=comparisons,component_mass_checks=mass_checks,
        neutral_color_error=color_error,rotation_errors=rotation,
        invalid_guards=guards,solver_sha256=study.solver_sha,
        scope='Neutral flat equal-tau1 D2Q9 only; both phase components and total populations retain impermeable reflection. Physical molecular slip is not calibrated.')
    (study.root/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(root=str(study.root),max_profile_error=max(r['relative_profile_error'] for r in hydraulic),color_error=color_error)))


if __name__=='__main__':main()
