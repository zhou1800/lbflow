#!/usr/bin/env python3
"""Recompute fixed-width front controls and reciprocal work from raw evidence."""
import argparse
import json
from pathlib import Path

import numpy as np

from benchmarks import geometry_fields, linear_crossings, sha
from capillary_response import measure as closed_measure
from collect_capillary_controls import verify
from front_mobility import donors, match_beta, zero_flux_profile, WIDTH_REFERENCE
from front_reciprocity import auxiliary
from phase_recurrence import bulk_profile, coefficients
from pressure_reservoirs import channel_measure, front_measure


CAMPAIGNS=('calibrate','closed','channels','fronts-coarse','fronts-fine',
    'holds-coarse','holds-fine','reciprocal-coarse','reciprocal-fine')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();root=Path(args.root);receipts={};provenance=[]
    for campaign in CAMPAIGNS:
        for path,r in donors(root/campaign):
            verify(r);folder=Path(r['output_files'][-1]).parent
            checkpoints=[]
            for p in sorted(folder.glob('*.rst*')):
                checkpoints.append(dict(path=f'{campaign}/cases/{folder.name}/{p.name}',sha256=sha(p)))
            if r['start']:
                for field,suffix in (('checkpoint_sha256','.rst'),('checkpoint_metadata_sha256','.rst.meta')):
                    if sha(folder/f'field.{r["start"]:010d}{suffix}')!=r[field]:
                        raise AssertionError('Copied outgoing checkpoint changed')
            provenance.append(dict(campaign=campaign,name=r['name'],receipt_sha256=sha(path),
                input_sha256=r['input_sha256'],geometry_sha256=r['geometry_sha256'],
                solver_sha256=r['solver_sha256'],start=r['start'],
                fields=[dict(path=f'{campaign}/cases/{folder.name}/{Path(p).name}',sha256=sha(p)) for p in r['output_files']],
                checkpoints=checkpoints,log_sha256=sha(folder/'run.log')))
            receipts[(campaign,r['name'])]=r
    calibration=[]
    for row in json.loads((root/'calibrate/calibration.json').read_text())['measurements']:
        r=receipts[('calibrate',row['name'])];cfg=r['parameters'];_,rho,chi,v=geometry_fields(r);profile=chi[0,0]
        middle=linear_crossings(profile,0)[0]
        a=min(linear_crossings(profile,.9),key=lambda x:abs(x-middle))
        b=min(linear_crossings(profile,-.9),key=lambda x:abs(x-middle));width=abs(a-b)
        beta=match_beta(WIDTH_REFERENCE*cfg['height']/16)
        ratio=coefficients()['bulk_segregation']/coefficients()['bulk_diffusion']
        marched=np.roll(bulk_profile(cfg['size'][0],beta/ratio),-cfg['size'][0]//4)
        scalar_error=float(abs((profile+1)/2-marched).max())
        if abs(width-row['width_09'])>2e-14 or beta!=cfg['beta'] or scalar_error!=row['max_scalar_profile_error'] or abs(profile-row['profile']).max()>2e-14:
            raise AssertionError('Independent width calibration does not reproduce')
        calibration.append(dict(row,independent_width=zero_flux_profile(beta)[0],
            relative_width_error=width/cfg['width_target']-1))
    closed=[closed_measure(r) for (campaign,name),r in receipts.items() if campaign=='closed']
    channels=[]
    for (campaign,name),r in receipts.items():
        if campaign!='channels':continue
        cfg=r['parameters'];row=channel_measure(r);_,rho,chi,v=geometry_fields(r)
        row.update(length_factor=cfg['length_factor'],pressure_difference=cfg['pressure_difference'],
            volume_conductance=float(v[0,1:-1,1:-1,0].mean()/cfg['pressure_difference']))
        channels.append(row)
    fronts=[]
    for (campaign,name),r in receipts.items():
        if not campaign.startswith(('fronts-','holds-')):continue
        row=front_measure(r);row.update(campaign=campaign,pressure_offset=r['parameters']['pressure_offset'])
        fronts.append(row)
    reciprocal=[]
    for campaign in ('reciprocal-coarse','reciprocal-fine'):
        data=json.loads((root/campaign/'reciprocity.json').read_text())
        driver=root/campaign/'.bin'/data['audit_driver_sha256']/'NumericalAudit'
        if sha(driver)!=data['audit_driver_sha256']:raise AssertionError('Frozen force audit executable changed')
        for row in data['measurements']:
            replay=receipts[(campaign,row['replay_receipt_name'])]
            seed_campaign=campaign.replace('reciprocal-','holds-')
            seed_path=root/seed_campaign/'results'/f'{row["name"]}.json'
            if sha(seed_path)!=row['seed_receipt_sha256']:
                raise AssertionError('Force audit no longer refers to the recorded held front')
            aux,flux,aux_record=auxiliary(root/'channels',row['height'],row['length_factor'])
            if abs(flux-row['auxiliary']['flux'])>2e-18 or abs(aux_record['volume_conductance']-row['auxiliary']['volume_conductance'])>2e-14:
                raise AssertionError('Independent hydraulic auxiliary does not reproduce')
            forces=[];speeds=[]
            for point in row['points']:
                t=point['step'];folder=root/campaign/'audit'/replay['name']/str(t)
                for field,filename in (('input_sha256','input.prm'),('checkpoint_sha256',f'field.{t:010d}.rst'),('log_sha256','audit.log')):
                    if sha(folder/filename)!=point[field]:raise AssertionError('Force audit evidence changed')
                audit=json.loads(next(line[len('CAPILLARY_STATE '):] for line in (folder/'audit.log').read_text().splitlines()
                    if line.startswith('CAPILLARY_STATE ')))
                force=np.array(audit['force']).reshape((2,)+tuple(audit['size']))[::-1].transpose(1,2,0)
                _,rho,chi,v=geometry_fields(replay,t);density=np.array(audit['density']).reshape(audit['size']);fluid=density>0
                color=np.array(audit['color']).reshape(audit['size'])
                if max(abs(density-rho[0]).max(),abs(color[fluid]-chi[0][fluid]).max())>3e-12:
                    raise AssertionError('Force generating state no longer matches physical fields')
                weighted=-float(np.sum(force*aux))/flux;speed=float(v[0,1:-1,1:-1,0].mean())
                if abs(weighted-point['capillary_pressure_weighted'])>2e-18 or abs(speed-point['physical_mean_speed'])>2e-18:
                    raise AssertionError('Reciprocal force work does not reproduce')
                forces.append(force);speeds.append(speed)
            work=-np.sum(np.mean(forces,axis=0)*aux,axis=2)/flux
            if max(abs(work.sum()-row['capillary_pressure_weighted']),
                abs(np.mean(speeds)-row['physical_mean_speed']),
                abs(work.sum(axis=0)-row['weighted_pressure_by_x']).max(),
                abs(work.sum(axis=1)-row['weighted_pressure_by_y']).max())>2e-18:
                raise AssertionError('Force-work spatial decomposition changed')
            reciprocal.append(dict(row,campaign=campaign,audit_driver_sha256=sha(driver)))
    fits=[]
    for h,l in sorted({(r['height'],r['length_factor']) for r in fronts}):
        g=float(np.mean([r['volume_conductance'] for r in channels if (r['height'],r['length_factor'])==(h,l)]))
        c=next(r for r in closed if r['height']==h);dx=16/h
        for held in (False,True):
            rows=sorted([r for r in fronts if (r['height'],r['length_factor'])==(h,l) and
                r['campaign'].startswith('holds-' if held else 'fronts-')],key=lambda r:r['pressure_offset'])
            for i in range(len(rows)//2):
                minus,plus=rows[i],rows[-1-i];dp=plus['pressure_difference']-minus['pressure_difference']
                u0,u1=minus['late_volume_front_speed'],plus['late_volume_front_speed'];m=(u1-u0)/dp
                if u0>=0 or u1<=0:raise AssertionError('Prespecified centered pair does not give both drift directions')
                history0,history1=minus['history'],plus['history']
                crossing_slopes={}
                for observable in ('centerline_front','first_fluid_row_front'):
                    rates=[]
                    for history in (history0,history1):
                        a,b=history[len(history)//2],history[-1]
                        rates.append((b[observable]-a[observable])/(b['step']-a['step']))
                    crossing_slopes[observable]=float((rates[1]-rates[0])/dp)
                record=dict(height=h,length_factor=l,held=held,offset=plus['pressure_offset'],
                    pure_mobility=g,front_mobility=m,relative_deficit=1-m/g,
                    extra_resistance=1/m-1/g,base_unit_extra_resistance=(1/m-1/g)/dx,
                    base_unit_front_mobility=m*dx,base_unit_pure_mobility=g*dx,
                    centerline_and_first_row_mobilities=crossing_slopes,
                    zero_pressure_ratio=(plus['pressure_ratio']*u0-minus['pressure_ratio']*u1)/(u0-u1),
                    closed_pressure_ratio=1+c['relative_pressure_error'],
                    names=[minus['name'],plus['name']])
                if held:
                    rr=sorted([r for r in reciprocal if (r['height'],r['length_factor'])==(h,l)
                        and abs(r['pressure_offset'])==plus['pressure_offset']],key=lambda r:r['pressure_offset'])
                    if len(rr)!=2:raise AssertionError('Each held pressure pair needs independent force audits')
                    a,b=rr;du=b['physical_mean_speed']-a['physical_mean_speed']
                    rmean=dp/du-1/g;rforce=(b['capillary_pressure_weighted']-a['capillary_pressure_weighted'])/du
                    record.update(mean_speed_extra_resistance=rmean,force_work_extra_resistance=rforce,
                        relative_force_resistance_residual=rforce/rmean-1,
                        reciprocal_incremental_speed_residual=g*(dp-(b['capillary_pressure_weighted']-a['capillary_pressure_weighted']))/du-1,
                        mean_speed_mobility=du/dp)
                fits.append(record)
    for row in fits:
        if row['held']:
            before=next(r for r in fits if not r['held'] and all(r[k]==row[k] for k in ('height','length_factor','offset')))
            row['relative_duration_mobility_change']=row['front_mobility']/before['front_mobility']-1
    data=dict(protocol='Fixed independently calibrated W09/H; unchanged projection/legacy MRT at tau1/theta60; prespecified rate, length and duration controls, with audited capillary force work.',
        physical_executions=len(provenance),calibration=calibration,closed=closed,channels=channels,
        fronts=fronts,reciprocal=reciprocal,response_fits=fits,provenance=provenance,
        solver_sha256=provenance[0]['solver_sha256'],
        units='Base length unit is the H16 grid. dx=16/H, dt=dx^2, sigma_lattice=.005*dx. U_base=U_lattice/dx, dp_base=dp_lattice/dx^2, M_base=M_lattice*dx, R_base=R_lattice/dx. Leading scalar diffusion and viscosity both remain 1/6 in base units. Matching width and leading diffusion does not verify a fixed effective contact-line slip length.',
        verification='Inputs, geometry, executable, every VTI and checkpoint hash checked. Closed/front/channel diagnostics recomputed; every audited force work recomputed and color/density matched to physical states. No velocity or force filtering in solver.',
        scope='Equal-density/equal-viscosity D2Q9, planar theta60 walls. Length localization tested at H16; rate interval .01 to .04 Young pressure. Finite-width grid sequence is not converged. Reciprocal work is a continuum diagnostic, not an exact discrete identity or a verified physical contact-line law.')
    if len({p['solver_sha256'] for p in provenance})!=1:raise AssertionError('Unexpected numerical kernel change')
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(physical_executions=len(provenance),fronts=len(fronts),response_fits=len(fits),
        maximum_force_resistance_residual=max(abs(r['relative_force_resistance_residual']) for r in fits if r['held']))))


if __name__=='__main__':main()
