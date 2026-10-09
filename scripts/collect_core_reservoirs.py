#!/usr/bin/env python3
"""Reanalyze core/operator and guarded reservoir campaigns with raw hashes."""
import argparse
import json
from pathlib import Path

import numpy as np

from benchmarks import geometry_fields, linear_crossings, measure_capillary, sha
from capillary_response import measure as slit_measure
from collect_capillary_controls import verify
from pressure_reservoirs import channel_measure, front_measure
from staggered_diagnostics import measure as adjacent_measure
from phase_recurrence import VELOCITIES


SLITS=('core-pilot','core-h16-holds','core-refinement','core-fine-holds')
OPEN=('open-front-pilot','open-brackets','open-bracket-holds','open-cosine-brackets',
      'open-cosine-holds','open-narrow-brackets','open-narrow-holds','open-force-removed')
OTHER=('core-calibration','core-channels','core-drops','core-drop-holds','oblique-calibration',
       'parity-forcing','core-adjacent','open-pilot','open-channels','open-channels-tau1')


def receipts(root):
    for path in sorted((root/'results').glob('*.json')):
        data=json.loads(path.read_text())
        if not isinstance(data,dict) or 'parameters' not in data:continue
        yield path,data


def key(cfg):
    return cfg.get('recoloring_model','cosine'),cfg.get('mrt_spectrum','legacy'),cfg['height']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True);parser.add_argument('--previous-data',required=True)
    parser.add_argument('--core-checks',required=True);parser.add_argument('--pressure-checks',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args();root=Path(args.root);provenance=[];by_name={}
    for sub in SLITS+OPEN+OTHER:
        for path,r in receipts(root/sub):
            verify(r);folder=Path(r['output_files'][-1]).parent
            if r['start']:
                for item,suffix in [('checkpoint_sha256','.rst'),('checkpoint_metadata_sha256','.rst.meta')]:
                    if sha(folder/f'field.{r["start"]:010d}{suffix}')!=r[item]:raise AssertionError('Copied native checkpoint changed')
            provenance.append(dict(root=sub,case=r['name'],receipt_sha256=sha(path),
                solver_sha256=r['solver_sha256'],input_sha256=r['input_sha256'],geometry_sha256=r['geometry_sha256'],
                start=r['start'],checkpoint_sha256=r.get('checkpoint_sha256'),
                checkpoint_metadata_sha256=r.get('checkpoint_metadata_sha256'),
                fields=[dict(path=f'{Path(f).parent.parent.parent.name}/cases/{Path(f).parent.name}/{Path(f).name}',sha256=sha(f)) for f in r['output_files']]))
            by_name[r['name']]=r
    calibration=json.loads((root/'core-calibration/calibration.json').read_text())
    for path,r in receipts(root/'core-calibration'):
        _,rho,chi,v=geometry_fields(r);profile=chi[0,0];mid=linear_crossings(profile,0)[0]
        a=min(linear_crossings(profile,.9),key=lambda x:abs(x-mid));b=min(linear_crossings(profile,-.9),key=lambda x:abs(x-mid))
        row=next(m for m in calibration['models'] if m['name']==r['name'])
        if abs(row['width_09']-abs(a-b))>2e-14 or abs(np.array(row['profile'])-profile).max()>2e-14:
            raise AssertionError('Independent planar calibration no longer reproduces')
    closed=[]
    for sub in SLITS:
        for path,r in receipts(root/sub):
            cfg=r['parameters'];row=next(m for m in calibration['models'] if m['recoloring_model']==cfg['recoloring_model']
                and m['mrt_spectrum']==cfg['mrt_spectrum'] and m['beta']==cfg['beta'])
            # Width is diagnostic metadata, never a simulation input. Preserve
            # original receipts and attach the independent correction explicitly.
            copy=dict(r,parameters=dict(cfg,planar_width_09=row['width_09']))
            result=slit_measure(copy)
            result.update(campaign=sub,receipt_sha256=sha(path),
                width_diagnostic=dict(original=cfg.get('planar_width_09'),used=row['width_09'],
                    calibration_sha256=sha(root/'core-calibration/calibration.json'),
                    reason='Crossings must surround the same interface; the original 64-node width metadata in some receipts selected opposite interfaces. No simulation input, raw field or pressure measurement is changed.'))
            closed.append(result)
    prior=json.loads(Path(args.previous_data).read_text())
    baseline=[r for r in prior['long_holds'] if r['height'] in (16,32,48) and r['angle']==60
        and r['beta']==.7 and r['wall_gradient']=='ghost' and r['phase_boundary']=='flux']
    if len(baseline)!=3:raise AssertionError('Need the complete recorded cosine baseline')
    adjacent=[]
    for path,r in receipts(root/'core-adjacent'):
        seed_name=r['parameters']['name'].removesuffix('-adjacent')+'-1rank'
        seed=by_name[seed_name];adjacent.append(adjacent_measure(seed,r))
    drop=[]
    for sub in ('core-drops','core-drop-holds'):
        for path,r in receipts(root/sub):
            original=r
            if r['start']:
                seed_name=r['parameters']['name'].removesuffix('-held')+'-1rank'
                seed=by_name[seed_name];original=dict(r,output_files=seed['output_files']+r['output_files'])
            row=measure_capillary(original);cfg=r['parameters']
            row.update(campaign=sub,recoloring_model=cfg['recoloring_model'],mrt_spectrum=cfg['mrt_spectrum'],
                late_color_change=float(abs(geometry_fields(r)[2]-geometry_fields(r,cfg['steps']-cfg['write'])[2]).max()))
            drop.append(row)
    channels=json.loads((root/'core-channels/channels.json').read_text())
    pressure_channels=[dict(channel_measure(r),campaign=sub) for sub in ('open-channels','open-channels-tau1') for p,r in receipts(root/sub)]
    fronts=[dict(front_measure(r),campaign=sub) for sub in OPEN for p,r in receipts(root/sub)]
    pressure_fits=[]
    for group in [('cosine','legacy'),('projection','legacy'),('projection','liu2015')]:
        model,spectrum=group
        for height in (16,32):
            records=[]
            for campaign in (('open-cosine-brackets','open-cosine-holds') if model=='cosine' else ('open-brackets','open-bracket-holds')):
                rows=sorted([r for r in fronts if r['campaign']==campaign and r['height']==height and
                    r['recoloring_model']==model and r['mrt_spectrum']==spectrum],key=lambda r:r['pressure_ratio'])
                if len(rows)!=3:raise AssertionError('Missing broad pressure/duration sequence')
                ratios=np.array([r['pressure_ratio'] for r in rows]);speed=np.array([r['late_volume_front_speed'] for r in rows])
                slope,intercept=np.linalg.lstsq(np.column_stack([ratios,np.ones(3)]),speed,rcond=None)[0]
                records.append(dict(campaign=campaign,zero_drift_ratio=float(-intercept/slope),
                    pressure_mobility=float(slope/rows[0]['young_pressure']),
                    max_linear_fit_speed_residual=float(abs(speed-(slope*ratios+intercept)).max()),
                    max_front_speed=float(abs(speed).max())))
            narrow=[]
            for campaign in ('open-narrow-brackets','open-narrow-holds'):
                rows=sorted([r for r in fronts if r['campaign']==campaign and r['height']==height and
                    r['recoloring_model']==model and r['mrt_spectrum']==spectrum],key=lambda r:r['pressure_ratio'])
                if len(rows)!=2 or rows[0]['late_volume_front_speed']>=0 or rows[1]['late_volume_front_speed']<=0:
                    raise AssertionError('Prespecified narrow pressure pair does not bracket both directions')
                narrow.append(dict(campaign=campaign,ratios=[r['pressure_ratio'] for r in rows],
                    measured_speeds=[r['late_volume_front_speed'] for r in rows]))
            current=next(r for r in (baseline if model=='cosine' else closed) if r['height']==height and
                r.get('mrt_spectrum','legacy')==spectrum and r.get('phase_boundary')=='flux' and
                (model=='cosine' or r['beta']!=.7) and
                (model=='cosine' or r['campaign'] in ('core-h16-holds','core-fine-holds')) and r['angle']==60)
            pure=[r for r in pressure_channels if r['tau']==1. and r['height']==height and r['mrt_spectrum']==spectrum]
            conductance=float(np.mean([np.mean(r['profile'])/(r['nominal_pressure_gradient']*(4*height-1)) for r in pure]))
            pressure_fits.append(dict(recoloring_model=model,mrt_spectrum=spectrum,height=height,
                records=records,narrow_brackets=narrow,closed_pressure_ratio=1+current['relative_pressure_error'],
                zero_minus_closed=records[-1]['zero_drift_ratio']-(1+current['relative_pressure_error']),
                duration_zero_ratio_change=records[-1]['zero_drift_ratio']-records[0]['zero_drift_ratio'],
                measured_pure_channel_mobility=conductance,
                relative_front_mobility_to_pure=records[-1]['pressure_mobility']/conductance-1))
    # A force-free, phase-containing control tests total-fluid decoupling.
    removed=[]
    for row in [r for r in fronts if r['campaign']=='open-force-removed']:
        h=row['height'];spectrum=row['mrt_spectrum'];pure=[r for r in pressure_channels if r['tau']==1. and r['height']==h and r['mrt_spectrum']==spectrum]
        conductance=np.mean([np.mean(r['profile'])/(r['nominal_pressure_gradient']*(4*h-1)) for r in pure])
        removed.append(dict(name=row['name'],recoloring_model=row['recoloring_model'],mrt_spectrum=spectrum,height=h,
            predicted_pure_speed=float(conductance*row['pressure_difference']),measured_mean_speed=row['history'][-1]['mean_speed'],
            measured_front_speed=row['late_volume_front_speed'],
            relative_mean_error=float(row['history'][-1]['mean_speed']/(conductance*row['pressure_difference'])-1),
            relative_front_error=float(row['late_volume_front_speed']/(conductance*row['pressure_difference'])-1),
            scope='Changed sigma=0 physics, not a capillary result or velocity filter.'))
    parity=json.loads((root/'parity-forcing/parity-forcing.json').read_text())
    seed=by_name['held-force-adjacent-1rank'];force_checks=[];points=[]
    driver=root/'parity-forcing/.bin'/parity['audit_driver_sha256']/'NumericalAudit'
    if sha(driver)!=parity['audit_driver_sha256']:raise AssertionError('Frozen force audit executable changed')
    for point in parity['force_identity_points']:
        t=point['step'];folder=root/'parity-forcing/audit/held-force-adjacent-1rank'/str(t)
        log=folder/'audit.log';audit=json.loads(next(line[len('CAPILLARY_STATE '):] for line in log.read_text().splitlines() if line.startswith('CAPILLARY_STATE ')))
        shape=tuple(audit['size']);rho=np.array(audit['density']).reshape(shape);fluid=rho>0;y,x=np.indices(shape)
        force=np.array(audit['force']).reshape((2,)+shape)[::-1]
        g=np.array([np.mean((force[d]*(-1.)**coordinate)[fluid]) for d,coordinate in enumerate((x,y))])
        pop=np.fromfile(folder/f'field.{t:010d}.rst',dtype='=f8').reshape((2,9)+shape)
        momentum=np.einsum('qyx,qd->dyx',pop.sum(axis=0),VELOCITIES)
        out=np.array([np.mean((momentum[d]*(-1.)**coordinate)[fluid]) for d,coordinate in enumerate((x,y))])
        _,r,chi,v=geometry_fields(seed,t);physical=np.array([np.mean((r[0]*v[0,...,d]*(-1.)**coordinate)[fluid]) for d,coordinate in enumerate((x,y))])
        if max(abs(g-point['force_projection']).max(),abs(out-point['outgoing_projection']).max(),abs(physical-point['physical_projection']).max())>2e-18:
            raise AssertionError('Recorded force/alternating momentum projections do not reproduce')
        if sha(folder/'input.prm')!=point['input_sha256'] or sha(folder/f'field.{t:010d}.rst')!=point['checkpoint_sha256']:
            raise AssertionError('Force audit inputs changed')
        force_checks.append(dict(step=t,log_sha256=sha(log),input_sha256=point['input_sha256'],checkpoint_sha256=point['checkpoint_sha256']))
        points.append(dict(point,physical_projection=physical.tolist(),force_projection=g.tolist(),outgoing_projection=out.tolist()))
    parity['force_identity_points']=points;parity['audit_file_checks']=force_checks
    parity['identity_residuals']=[(np.array(b['physical_projection'])+a['physical_projection']-.5*(np.array(b['force_projection'])-a['force_projection'])).tolist() for a,b in zip(points,points[1:])]
    if abs(np.array(parity['identity_residuals'])).max()>3e-12:raise AssertionError('Alternating momentum identity failed reanalysis')
    parity['translation_sampling_caveat']='Fields translate exactly. The inherited fixed pressure windows do not translate, causing a 4.739e-6 relative sampled difference; this is not a physical translation defect.'
    data=dict(protocol='Published core controls with exact axial transport matching; parity-force identity; fixed-spectrum hydraulic offset; guarded pressure-collar fronts with prespecified pressure pairs and four-duration checks.',
        calibration=calibration,oblique=json.loads((root/'oblique-calibration/oblique.json').read_text()),
        closed=closed,baseline=baseline,adjacent=adjacent,drops=drop,force_channels=channels,
        pressure_channels=pressure_channels,fronts=fronts,pressure_fits=pressure_fits,force_removed_controls=removed,
        parity=parity,core_contracts=json.loads(Path(args.core_checks).read_text()),
        pressure_contracts=json.loads(Path(args.pressure_checks).read_text()),
        closed_compatibility=json.loads((root/'core-compatibility.json').read_text()),
        previous_data_sha256=sha(args.previous_data),provenance=provenance,physical_executions=len(provenance),
        scope='Equal-density D2Q9, flat nonsingular walls; long physical fronts at equal viscosities and theta60, independently matched beta. Joint grid/sharp-interface sequence. No full 2015 algorithm replica, dynamic contact-line law, curved-3D wetting, throat trapping or hysteresis claim.',
        verification='Raw input/geometry/executable/checkpoint hashes checked. Slit/pressure/front/drop diagnostics recomputed. Original diagnostic width metadata retained with independent correction.')
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(physical_executions=len(provenance),closed=len(closed),fronts=len(fronts),pressure_fits=len(pressure_fits)),indent=2))


if __name__=='__main__':main()
