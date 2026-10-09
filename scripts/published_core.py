#!/usr/bin/env python3
"""Compare specified published recoloring/MRT operators at matched transport."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path

import numpy as np

from benchmarks import Study, geometry_fields, linear_crossings, measure_capillary as bulk_measure, spec
from capillary_response import cases, measure
from phase_recurrence import coefficients, VELOCITIES, WEIGHTS
from porous_media import closed_history


def models():
    ratio=coefficients()['bulk_segregation']/coefficients()['bulk_diffusion']
    return [('cosine',.7,'cosine'),('projection',.7,'published-beta'),
            ('projection',.7*ratio,'matched-transport')]


def calibration(study, jobs):
    configs=[]
    for model,beta,label in models():
        for spectrum in ('legacy','liu2015'):
            configs.append(spec(f'core-planar-{label}-{spectrum}',nx=128,ny=4,nz=1,dimensions=2,
                tau=(1.,1.),capillary_model='csf',steps=8000,write=4000,recoloring_model=model,
                mrt_spectrum=spectrum,beta=beta,label=lambda x,y,z:0 if x<64 else 100,
                indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"']))
    with ThreadPoolExecutor(max_workers=jobs) as pool:receipts=list(pool.map(study.run,configs))
    records=[]
    for receipt in receipts:
        cfg=receipt['parameters'];_,rho,chi,v=geometry_fields(receipt)
        profile=chi[0,0];middle=linear_crossings(profile,0)[0]
        a=min(linear_crossings(profile,.9),key=lambda x:abs(x-middle))
        b=min(linear_crossings(profile,-.9),key=lambda x:abs(x-middle))
        if not 1<abs(a-b)<16:raise AssertionError('Planar width selected different interfaces')
        records.append(dict(name=receipt['name'],recoloring_model=cfg['recoloring_model'],beta=cfg['beta'],
            mrt_spectrum=cfg['mrt_spectrum'],width_09=abs(a-b),profile=profile.tolist(),
            max_speed=float(np.linalg.norm(v,axis=-1).max()),
            late_color_change=float(abs(chi-geometry_fields(receipt,4000)[2]).max()),
            density_range=[float(rho.min()),float(rho.max())],mass=closed_history(receipt)['max_relative_mass_drift']))
    errors=[]
    for spectrum in ('legacy','liu2015'):
        old=next(r for r in records if r['recoloring_model']=='cosine' and r['mrt_spectrum']==spectrum)
        matched=next(r for r in records if r['recoloring_model']=='projection' and r['beta']!=.7 and r['mrt_spectrum']==spectrum)
        error=float(abs(np.array(old['profile'])-matched['profile']).max())
        if error>2e-11:raise AssertionError('Analytically matched planar profiles differ')
        errors.append(dict(mrt_spectrum=spectrum,max_matched_profile_error=error))
    tensors={}
    for model in ('cosine','projection'):
        length=np.linalg.norm(VELOCITIES,axis=1);factor=np.ones(9)
        if model=='cosine':factor=np.divide(1.,length,out=np.zeros(9),where=length>0)
        weight=WEIGHTS*factor;second=np.einsum('q,qa,qb->ab',weight,VELOCITIES,VELOCITIES)
        fourth=np.einsum('q,qa,qb,qc,qd->abcd',weight,VELOCITIES,VELOCITIES,VELOCITIES,VELOCITIES)
        tensors[model]=dict(second=second.tolist(),fourth_xxxx=float(fourth[0,0,0,0]),
                            fourth_xxyy=float(fourth[0,0,1,1]),fourth_isotropy_defect=float(fourth[0,0,0,0]-3*fourth[0,0,1,1]))
    result=dict(protocol='Independent bulk planar calibration, exact one-dimensional beta match and recoloring moment tensors.',
                models=records,matched_profile_checks=errors,tensors=tensors,solver_sha256=study.solver_sha,
                scope='Specified published core operators; curvature/wall-normal continuation remains LBFlow-defined, not an exact full 2015 replication.')
    (study.root/'calibration.json').write_text(json.dumps(result,indent=2)+'\n')
    return receipts


def slits(study,args):
    data=json.loads(Path(args.calibration).read_text());configs=[]
    for model,beta,label in models():
        if label not in args.models:continue
        for spectrum in args.spectra:
            row=next(r for r in data['models'] if r['recoloring_model']==model and r['beta']==beta and r['mrt_spectrum']==spectrum)
            for cfg in cases(args.heights,args.angles,['ghost'],beta,capillary_models=['csf']):
                for boundary in args.boundaries:
                    configs.append(dict(cfg,name=cfg['name']+'-'+label+'-'+spectrum+'-'+boundary,
                        phase_boundary=boundary,recoloring_model=model,mrt_spectrum=spectrum,planar_width_09=row['width_09']))
    def run(cfg):
        receipt=study.run(cfg);result=measure(receipt)
        (study.root/'results'/f'{receipt["name"]}-measurement.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(name=receipt['name'],pressure_error=result['relative_pressure_error'],
                             max_speed=result['history'][-1]['max_speed'])),flush=True)
        return receipt
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:return list(pool.map(run,configs))


def channels(study,args):
    configs=[]
    for h in (8,16,32):
        for tau in (.6,1.,1.55):
            nu=(tau-.5)/3;steps=500*math.ceil(24*h*h/(nu*math.pi**2)/500)
            for spectrum in ('legacy','liu2015'):
                configs.append(spec(f'core-channel-H{h}-tau{tau:g}-{spectrum}',nx=4,ny=h+2,nz=1,
                    dimensions=2,tau=(.6,tau),steps=steps,write=steps//2,restart=steps,
                    mrt_spectrum=spectrum,recoloring_model='projection',capillary_model='csf',
                    wall_gradient='ghost',sigma=0,height=h,fluid_tau=tau,force=1e-6,
                    label=lambda x,y,z,h=h:255 if y in (0,h+1) else 75,
                    indices=['75 "GV" "0.000001" "1" "1"','255 "SD" "60" "N/A" "N/A"']))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(study.run,configs))
    records=[]
    for receipt in receipts:
        cfg=receipt['parameters'];h=cfg['height'];tau=cfg['fluid_tau'];nu=(tau-.5)/3
        sq=1.92 if cfg['mrt_spectrum']=='liu2015' else 8*(2-1/tau)/(8-1/tau)
        magic=(tau-.5)*(1/sq-.5);offset=4*magic/3-.25
        _,rho,chi,v=geometry_fields(receipt);profile=v[0,1:-1,:,0].mean(axis=1)
        y=np.arange(h)+.5;reference=cfg['force']/(2*nu)*(y*(h-y)+offset)
        unshifted=cfg['force']/(2*nu)*y*(h-y)
        error=float(abs(profile-reference).max()/abs(reference).max())
        if error>3e-6:raise AssertionError('Published MRT hydraulic wall offset does not match the independent reference')
        records.append(dict(name=receipt['name'],height=h,tau=tau,mrt_spectrum=cfg['mrt_spectrum'],
            magic_parameter=magic,constant_offset=offset,profile=profile.tolist(),reference=reference.tolist(),
            relative_profile_error=error,relative_unshifted_flux_error=float(profile.mean()/unshifted.mean()-1),
            predicted_unshifted_flux_error=float(12*offset/(2*h*h+1)),
            density_range=[float(rho[rho>0].min()),float(rho.max())],
            mass=closed_history(receipt)['max_relative_mass_drift']))
    (study.root/'channels.json').write_text(json.dumps(dict(protocol='Pure second-phase force-driven channel; prescribed MRT spectrum, physical half-force velocity and independent midpoint offset.',measurements=records,solver_sha256=study.solver_sha),indent=2)+'\n')
    return receipts


def oblique(study,args):
    configs=[]
    for model,beta,label in models():
        configs.append(spec('oblique-'+label,nx=128,ny=128,nz=1,dimensions=2,
            capillary_model='csf',recoloring_model=model,beta=beta,steps=8000,write=4000,
            label=lambda x,y,z:0 if (x+y)%128<64 else 100,
            indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"']))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(study.run,configs))
    records=[]
    for r in receipts:
        _,rho,chi,v=geometry_fields(r);profile=chi[0,0]
        mid=linear_crossings(profile,0)[0]
        left=min(linear_crossings(profile,.9),key=lambda x:abs(x-mid))
        right=min(linear_crossings(profile,-.9),key=lambda x:abs(x-mid))
        records.append(dict(name=r['name'],recoloring_model=r['parameters']['recoloring_model'],
            beta=r['parameters']['beta'],width_09_normal=abs(right-left)/math.sqrt(2),
            max_speed=float(np.linalg.norm(v,axis=-1).max()),profile=profile.tolist()))
    (study.root/'oblique.json').write_text(json.dumps(dict(normal_orientation_degrees=45,
        measurements=records,solver_sha256=study.solver_sha,
        sampling_caveat='Linear crossing samples are spaced by 1/sqrt(2) along the normal; axial spacing is 1. Width differences include interpolation bias.'),indent=2)+'\n')
    return receipts


def drops(study,args):
    configs=[]
    for radius in (12,24):
        n=3*radius+8;center=(n-1)/2
        for model,beta,label in models():
            if label not in ('cosine','matched-transport'):continue
            for spectrum in ('legacy','liu2015'):
                configs.append(spec(f'core-laplace-R{radius}-{label}-{spectrum}',nx=n,ny=n,nz=1,
                    dimensions=2,radius=radius,steps=16000,write=4000,
                    capillary_model='csf',recoloring_model=model,mrt_spectrum=spectrum,beta=beta,
                    label=lambda x,y,z,c=center,r=radius:100 if (x-c)**2+(y-c)**2<=r*r else 0,
                    indices=['0 "VD" "N/A" "N/A" "0"','100 "VD" "N/A" "N/A" "1"']))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:receipts=list(pool.map(study.run,configs))
    records=[]
    for r in receipts:
        row=bulk_measure(r);cfg=r['parameters']
        row.update(recoloring_model=cfg['recoloring_model'],mrt_spectrum=cfg['mrt_spectrum'],beta=cfg['beta'])
        last=geometry_fields(r)[2];previous=geometry_fields(r,12000)[2]
        row['late_color_change']=float(abs(last-previous).max())
        records.append(row)
    (study.root/'drops.json').write_text(json.dumps(dict(measurements=records,solver_sha256=study.solver_sha,
        protocol='Bulk droplet pressure/curvature calibration; input surface tension unchanged across independently matched transport controls.'),indent=2)+'\n')
    return receipts


def hold_drops(study,args):
    seeds=json.loads(Path(args.seeds).read_text())
    seeds=[r for r in seeds if r['parameters'].get('radius')==24]
    def run(seed):
        cfg=dict(seed['parameters']);folder=Path(seed['output_files'][-1]).parent;cut=cfg['steps']
        cfg.update(name=cfg['name']+'-held',steps=64000,write=16000,restart=64000,geometry_file=str(folder/'geometry.raw'))
        receipt=study.run(cfg,start=cut,checkpoint=folder/f'field.{cut:010d}')
        joined=dict(receipt,output_files=seed['output_files']+receipt['output_files'])
        row=bulk_measure(joined);row.update(seed_name=seed['name'],seed_checkpoint_sha256=receipt['checkpoint_sha256'],
            recoloring_model=cfg['recoloring_model'],mrt_spectrum=cfg['mrt_spectrum'],beta=cfg['beta'],
            late_color_change=float(abs(geometry_fields(receipt)[2]-geometry_fields(receipt,48000)[2]).max()))
        return receipt,row
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:runs=list(pool.map(run,seeds))
    (study.root/'drops.json').write_text(json.dumps(dict(measurements=[r[1] for r in runs],solver_sha256=study.solver_sha,
        protocol='Four-duration continuation of R24 droplets; initial mass from unchanged source snapshots.',seeds_sha256=__import__('benchmarks').sha(args.seeds)),indent=2)+'\n')
    return [r[0] for r in runs]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['calibrate','slits','channels','oblique','drops','hold-drops'])
    parser.add_argument('--solver',required=True);parser.add_argument('--root',required=True)
    parser.add_argument('--calibration');parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--seeds')
    parser.add_argument('--heights',nargs='+',type=int,default=[16])
    parser.add_argument('--angles',nargs='+',type=float,default=[60,90,120])
    parser.add_argument('--spectra',nargs='+',choices=['legacy','liu2015'],default=['legacy','liu2015'])
    parser.add_argument('--boundaries',nargs='+',choices=['bounce','flux'],default=['bounce','flux'])
    parser.add_argument('--models',nargs='+',choices=['cosine','published-beta','matched-transport'],default=['cosine','published-beta','matched-transport'])
    args=parser.parse_args()
    if args.stage=='slits' and not args.calibration:parser.error('Slit comparison requires independent calibration')
    if args.stage=='hold-drops' and not args.seeds:parser.error('Droplet holds require recorded seed receipts')
    study=Study(args.solver,args.root)
    action={'calibrate':lambda:calibration(study,args.jobs),'channels':lambda:channels(study,args),
            'slits':lambda:slits(study,args),'oblique':lambda:oblique(study,args),'drops':lambda:drops(study,args),
            'hold-drops':lambda:hold_drops(study,args)}
    receipts=action[args.stage]()
    (study.root/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')


if __name__=='__main__':main()
