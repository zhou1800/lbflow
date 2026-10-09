#!/usr/bin/env python3
"""Diagnose front resistance with audited forces and independent hydraulic fields.

The reciprocal work relation is a low-speed continuum diagnostic, not an
exact discrete LBM identity or a fitted moving-contact-line law.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil

import numpy as np

from benchmarks import Study, geometry_fields, sha
from front_mobility import donors
from parity_forcing import adjacent, force_audit


def auxiliary(root, height, length_factor):
    pair=[]
    for path,r in donors(root):
        cfg=r['parameters']
        if cfg['height']==height and cfg['length_factor']==length_factor:
            pair.append((path,r))
    if len(pair)!=2:raise ValueError('Reciprocity requires independent positive/negative pure controls')
    pair.sort(key=lambda item:item[1]['parameters']['pressure_difference'])
    velocities=[geometry_fields(r)[3][0] for path,r in pair]
    velocity=(velocities[1]-velocities[0])/2
    dp=(pair[1][1]['parameters']['pressure_difference']-pair[0][1]['parameters']['pressure_difference'])/2
    x=velocity.shape[1]//2
    flux=float(velocity[1:-1,x,0].sum())
    volume=float(velocity[1:-1,1:-1,0].mean())
    return velocity,flux,dict(height=height,length_factor=length_factor,
        pressure_difference=dp,flux=flux,volume_speed=volume,
        volume_conductance=volume/dp,flux_conductance=flux/dp,
        receipts=[dict(name=r['name'],sha256=sha(path)) for path,r in pair])


def diagnose(study,seed,path,driver,channels_root):
    cfg=seed['parameters'];h=cfg['height'];length=cfg['length_factor']
    aux,flux,aux_record=auxiliary(channels_root,h,length)
    replay=adjacent(study,seed,cfg['name']+'-reciprocal-adjacent',steps=2)
    points=[];force_fields=[];velocity_fields=[]
    for time in (cfg['steps']+1,cfg['steps']+2):
        audit=force_audit(study,replay,driver,time)
        folder=study.root/'audit'/replay['name']/str(time)
        data=json.loads(next(line[len('CAPILLARY_STATE '):]
            for line in (folder/'audit.log').read_text().splitlines()
            if line.startswith('CAPILLARY_STATE ')))
        shape=tuple(data['size'])
        force=np.array(data['force']).reshape((2,)+shape)[::-1].transpose(1,2,0)
        _,rho,chi,velocity=geometry_fields(replay,time)
        color=np.array(data['color']).reshape(shape)
        density=np.array(data['density']).reshape(shape)
        fluid=density>0
        state_error=max(float(abs(density-rho[0]).max()),float(abs(color[fluid]-chi[0][fluid]).max()))
        if state_error>3e-12:raise AssertionError('Outgoing audit differs from the physical force-generating scalar state')
        weighted=-float(np.sum(force*aux))/flux
        mean=float(velocity[0,1:-1,1:-1,0].mean())
        predicted=aux_record['volume_conductance']*(cfg['pressure_difference']-weighted)
        points.append(dict(audit,capillary_pressure_weighted=weighted,
            physical_mean_speed=mean,reciprocal_prediction=predicted,
            prediction_minus_raw=predicted-mean,max_force_state_error=state_error,
            max_speed=float(np.linalg.norm(velocity,axis=-1).max()),
            max_density_deviation=float(abs(rho[0][fluid]-1).max()),
            log_sha256=sha(folder/'audit.log')))
        force_fields.append(force);velocity_fields.append(velocity[0])
    force=np.mean(force_fields,axis=0);velocity=np.mean(velocity_fields,axis=0)
    work=-np.sum(force*aux,axis=2)/flux
    weighted=float(work.sum());mean=float(velocity[1:-1,1:-1,0].mean())
    predicted=aux_record['volume_conductance']*(cfg['pressure_difference']-weighted)
    return dict(name=seed['name'],height=h,length_factor=length,
        pressure_offset=cfg['pressure_offset'],pressure_difference=cfg['pressure_difference'],
        young_pressure=cfg['young_pressure'],seed_receipt_sha256=sha(path),
        replay_receipt_name=replay['name'],auxiliary=aux_record,points=points,
        capillary_pressure_weighted=weighted,physical_mean_speed=mean,
        reciprocal_prediction=predicted,prediction_minus_raw=predicted-mean,
        weighted_pressure_by_x=work.sum(axis=0).tolist(),
        weighted_pressure_by_y=work.sum(axis=1).tolist(),
        convention='Two adjacent physical states averaged only for analysis. Force is reconstructed from each outgoing checkpoint and its color/density checked against the physical VTI state. Auxiliary is the signed half-difference of independent pure-channel velocities. No solver filtering or force fitting.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',required=True);parser.add_argument('--audit-driver',required=True)
    parser.add_argument('--root',required=True);parser.add_argument('--front-root',required=True)
    parser.add_argument('--channels-root',required=True)
    parser.add_argument('--heights',type=int,nargs='+',default=[16,32,48])
    parser.add_argument('--jobs',type=int,default=3)
    args=parser.parse_args();study=Study(args.solver,args.root)
    source=Path(args.audit_driver).resolve();driver=study.root/'.bin'/sha(source)/source.name
    driver.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,driver);driver.chmod(source.stat().st_mode)
    seeds=[(r,path) for path,r in donors(args.front_root)
        if r['exit']==0 and r['finished'] and r['parameters']['height'] in args.heights]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows=list(pool.map(lambda item:diagnose(study,item[0],item[1],driver,args.channels_root),seeds))
    data=dict(measurements=rows,solver_sha256=study.solver_sha,audit_driver_sha256=sha(driver),
        protocol='Independent pure-channel auxiliary and actual capillary force work, with adjacent physical-state checks.',
        scope='Continuum creeping-flow reciprocity used as a diagnostic. Finite compressibility, inertia, lattice and pressure-boundary effects remain. It is not an exact discrete identity or a general dynamic contact-angle law.')
    (study.root/'reciprocity.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(cases=len(rows),root=str(study.root),
        maximum_absolute_prediction_residual=max(abs(r['prediction_minus_raw']) for r in rows))))


if __name__=='__main__':main()
