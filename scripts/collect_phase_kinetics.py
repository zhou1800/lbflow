#!/usr/bin/env python3
"""Recompute phase spectra, neutral displacement and independent force controls."""
import argparse
import json
from pathlib import Path

import numpy as np

from benchmarks import compare, geometry_fields, sha
from collect_capillary_controls import verify
from force_stokes_audit import remeasure as force_measure
from periodic_contact import measure
from phase_wave_spectrum import matrix, spectrum, expansion_controls
from stokes_phase_reference import observation, projection_rest
from stokes_reference_checks import check as reference_check


CAMPAIGNS=('production-pilot','pure-pilot','waves-pilot','production-fine',
    'pure-fine','waves-fine','waves-mode2','production-slip1-forward-H16',
    'production-slip1-forward-fine','pure-slip1-forward-fine',
    'production-slip1-holds','compat-default-old','compat-default-new',
    'hydraulic-audit-canonical')
REFERENCES=('reference-slip1-H16','reference-slip1-H32','reference-slip1-H48',
    'reference-slip1-H32-halfstep-negative','reference-zero-sigma-H16',
    'reference-zero-sigma-H32','hybrid-slip1-H16','hybrid-slip1-H32',
    'hybrid-slip1-H48-canonical','hybrid-slip1-H32-halfstep',
    'reference-slip1-H32-projection-seed','hybrid-slip1-H32-logistic-seed')


def receipt_at_root(receipt, root):
    """Relocate paths in memory while keeping original receipt/input bytes."""
    folder=Path(receipt['output_files'][-1]).parent
    if folder.parent.name!='cases':raise ValueError('Unexpected recorded case topology')
    original=str(folder.parents[2])+'/'
    destination=str(root)+'/'
    def moved(value):
        if isinstance(value,dict):return {k:moved(v) for k,v in value.items()}
        if isinstance(value,list):return [moved(v) for v in value]
        if isinstance(value,str) and value.startswith(original):return destination+value[len(original):]
        return value
    return moved(receipt)


def reference_measure(folder):
    original=json.loads((folder/'result.json').read_text());p=original['parameters']
    frozen=folder/'source';source_hashes={f.name:sha(f) for f in frozen.glob('*.py')}
    if source_hashes['stokes_phase_reference.py']!=original['source_sha256']:
        raise AssertionError('Reference receipt does not identify its frozen source')
    if 'source_files_sha256' in original and source_hashes!=original['source_files_sha256']:
        raise AssertionError('Frozen reference dependency changed')
    history=[]
    for f in sorted(folder.glob('field.*.npz')):
        with np.load(f) as field:
            phi,u,v=field['phi'],field['u'],field['v'];h=float(field['h']);t=float(field['time'])
            row=observation(phi,u,v,h,t)
            row['max_divergence']=float(abs((np.roll(u,-1,axis=1)-u)/h+np.diff(v,axis=0)/h).max())
            history.append(row)
    if len(history)!=len(original['history']):raise AssertionError('Missing reference fields')
    positions=np.unwrap([r['centroid_modulo']*2*np.pi/64 for r in history])*64/(2*np.pi)
    mass0=history[0]['phase_mass']
    for row,x,old in zip(history,positions,original['history']):
        row['centroid_unwrapped']=float(x);row['relative_mass_change']=(row['phase_mass']-mass0)/mass0
        if max(abs(row[k]-old[k]) for k in row if isinstance(row[k],(int,float)))>3e-13:
            raise AssertionError('Saved reference fields do not reproduce recorded measurements')
        if row['phase_range']!=old['phase_range']:raise AssertionError('Reference range changed')
    first,last=history[len(history)//2],history[-1]
    speed=(last['centroid_unwrapped']-first['centroid_unwrapped'])/(last['time']-first['time'])
    if abs(speed-original['late_centroid_speed'])>2e-16:raise AssertionError('Reference centroid speed changed')
    scope=('Independent Stokes solve with retained projection phase stencil: kinetic ablation, not independent continuum validation.'
        if p.get('phase_stencil','fv')=='projection' else original['scope'])
    return dict(name=folder.name,parameters=p,history=history,late_centroid_speed=speed,
        late_mean_speed=last['mean_speed'],source_hashes=source_hashes,scope=scope,recorded_scope=original['scope'],
        relative_last_interval_mean_change=last['mean_speed']/history[-2]['mean_speed']-1)


def collect(root):
    receipts={};rows=[]
    for campaign in CAMPAIGNS:
        for path in sorted((root/campaign/'results').glob('*.json')):
            r=receipt_at_root(json.loads(path.read_text()),root);verify(r);folder=Path(r['output_files'][-1]).parent
            if r['start']:
                for key,suffix in (('checkpoint_sha256','.rst'),('checkpoint_metadata_sha256','.rst.meta')):
                    if sha(folder/f'field.{r["start"]:010d}{suffix}')!=r[key]:raise AssertionError('Outgoing donor changed')
            cfg=r['parameters']
            if cfg.get('seed_phase_sha256') and sha(cfg['seed_phase_file'])!=cfg['seed_phase_sha256']:
                raise AssertionError('Recorded initial phase profile changed')
            receipts[(campaign,r['name'])]=r
            if not campaign.startswith(('compat-','hydraulic-')):
                row=measure(r);row['campaign']=campaign
                mass0=row['history'][0]['phase_mass'];mass1=row['history'][-1]['phase_mass']
                if mass0==0 and mass1!=0:raise AssertionError('Initially absent phase acquired mass')
                row['relative_phase_mass_change']=mass1/mass0-1 if mass0 else 0.
                rows.append(row)
    waves=[];spectra={}
    for row in rows:
        if not row['campaign'].startswith('waves-'):continue
        p=row['parameters'];h,m=p['height'],p['transverse_mode']
        key=(h,m)
        if key not in spectra:spectra[key]=spectrum(h,m)
        late=row['history'][len(row['history'])//2:]
        times=np.array([a['time'] for a in late]);amplitudes=np.array([a['projected_displacement'] for a in late])
        fit=np.polyfit(times,np.log(abs(amplitudes)),1);rate=float(-fit[0])
        predicted=spectra[key]['physical_relaxation_rate']
        waves.append(dict(name=row['name'],height=h,mode=m,amplitude=p['physical_amplitude'],
            measured_rate=rate,predicted_rate=predicted,relative_rate_difference=rate/predicted-1,
            relative_phase_mass_change=row['relative_phase_mass_change'],history=row['history']))
    symbol_checks=[]
    for h in (16,32,48):
        s=spectra[(h,1)];p=np.array(s['profile']);v=np.array(s['eigenvector'])
        c=np.cos(np.pi*(np.arange(h)+.5)/h)[:,None];epsilon=1e-5
        difference=(projection_rest(p[None]+epsilon*c*v[None],s['beta'])-
            projection_rest(p[None]-epsilon*c*v[None],s['beta']))/(2*epsilon)
        exact=c*(matrix(p,s['beta'],np.pi/h)@v)[None]
        error=float(abs(difference-exact).max())
        if error>1e-9:raise AssertionError('Fourier matrix is not the full nonlinear stencil derivative')
        symbol_checks.append(dict(height=h,epsilon=epsilon,max_absolute_error=error))
    references=[reference_measure(root/name) for name in REFERENCES]
    initial_controls=[]
    for changed,baseline in (('reference-slip1-H32-projection-seed','reference-slip1-H32'),
        ('hybrid-slip1-H32-logistic-seed','hybrid-slip1-H32')):
        a=next(r for r in references if r['name']==changed);b=next(r for r in references if r['name']==baseline)
        initial_controls.append(dict(changed=changed,baseline=baseline,
            relative_mean_speed_change=a['late_mean_speed']/b['late_mean_speed']-1,
            relative_centroid_speed_change=a['late_centroid_speed']/b['late_centroid_speed']-1))
    dynamic=[]
    for h in (16,32,48):
        lbm=next(r for r in rows if r['campaign'].startswith('production-slip1-forward')
            and r['parameters']['height']==h and r['parameters']['physical_force']==1e-7)
        fv=next(r for r in references if r['name']==f'reference-slip1-H{h}')
        hybrid=next(r for r in references if r['name']==f'hybrid-slip1-H{h}'+('-canonical' if h==48 else ''))
        dynamic.append(dict(height=h,lbm_centroid_speed=lbm['late_centroid_speed'],lbm_mean_speed=lbm['late_mean_speed'],
            fv_centroid_speed=fv['late_centroid_speed'],fv_mean_speed=fv['late_mean_speed'],
            hybrid_centroid_speed=hybrid['late_centroid_speed'],hybrid_mean_speed=hybrid['late_mean_speed'],
            relative_lbm_fv_mean_difference=lbm['late_mean_speed']/fv['late_mean_speed']-1,
            relative_hybrid_lbm_mean_difference=hybrid['late_mean_speed']/lbm['late_mean_speed']-1,
            relative_hybrid_lbm_centroid_difference=hybrid['late_centroid_speed']/lbm['late_centroid_speed']-1))
    durations=[]
    for hold in (r for r in rows if r['campaign']=='production-slip1-holds'):
        p=hold['parameters']
        donor=next(r for r in rows if r['campaign'].startswith('production-slip1-forward')
            and all(r['parameters'][k]==p[k] for k in ('height','physical_force')))
        durations.append(dict(height=p['height'],force=p['physical_force'],
            relative_centroid_speed_change=hold['late_centroid_speed']/donor['late_centroid_speed']-1,
            relative_mean_speed_change=hold['late_mean_speed']/donor['late_mean_speed']-1))
    rate_control=[]
    for force in (1e-7,2e-7):
        pair=sorted([r for r in rows if r['campaign']=='production-slip1-forward-H16'
            and abs(r['parameters']['physical_force'])==force],key=lambda r:r['parameters']['physical_force'])
        rate_control.append(dict(force_magnitude=force,
            mean_mobility=(pair[1]['late_mean_speed']-pair[0]['late_mean_speed'])/(2*force),
            centroid_mobility=(pair[1]['late_centroid_speed']-pair[0]['late_centroid_speed'])/(2*force)))
    audit=json.loads((root/'hydraulic-audit-canonical/force-stokes.json').read_text())
    driver=root/'hydraulic-audit-canonical/.bin'/audit['audit_driver_sha256']/'NumericalAudit'
    if sha(driver)!=audit['audit_driver_sha256']:raise AssertionError('Frozen force-audit executable changed')
    force_rows=[]
    for old in audit['measurements']:
        receipt=receipts[('hydraulic-audit-canonical',old['replay_receipt_name'])]
        fresh=force_measure(receipt,root/'hydraulic-audit-canonical')
        for a,b in zip(fresh['points'],old['points']):
            if a!=b:raise AssertionError('Actual force/Stokes audit does not reproduce')
        force_rows.append(fresh)
    compatibility=[]
    for (campaign,name),left in receipts.items():
        if campaign!='compat-default-old':continue
        right=receipts[('compat-default-new',name)];result=compare(left,right)
        left_folder=Path(left['output_files'][-1]).parent;right_folder=Path(right['output_files'][-1]).parent
        step=left['parameters']['steps']
        result['bitwise_populations']=(left_folder/f'field.{step:010d}.rst').read_bytes()==(right_folder/f'field.{step:010d}.rst').read_bytes()
        if not result['bitwise_fields'] or not result['bitwise_populations']:raise AssertionError('Default closure changed')
        compatibility.append(result)
    checks=reference_check()
    zeros=[]
    for ref in (r for r in references if r['parameters']['sigma']==0):
        p=ref['parameters'];h=p['h'];force=p['force'];expected=3*force*(16**2/6+1/3*h*h+16*p['slip'])
        error=ref['late_mean_speed']/expected-1
        if abs(error)>2e-12:raise AssertionError('Two-phase zero-sigma hydraulic decoupling fails')
        zeros.append(dict(height=p['ny'],relative_mean_speed_error=error))
    # Hash every retained file, including failed/exploratory receipts. Only
    # CAMPAIGNS/REFERENCES above contribute to publication measurements.
    manifest=[dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=sha(p))
        for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts
        and p.name not in ('phase-kinetics.json','reanalyzed-phase-kinetics.json','raw-manifest.json')
        and not p.name.endswith('.tex.build.log') and p.suffix!='.zip']
    data=dict(protocol='Fixed physical width and D; exact zero-flow phase spectrum, neutral periodic displacement with prescribed Navier slip, independent Stokes/FV reference and shared-phase kinetic ablation.',
        units='Physical height 16 and length 64; h=16/H, dt_lbm=h^2, nu=D=1/6, sigma=.005, body force=1e-7, prescribed hydraulic slip=1. Width W09=4.765235021842287. Waves use zero sigma/force and default zero slip.',
        production_executions=len(receipts),waves=waves,spectra=list(spectra.values()),symbol_checks=symbol_checks,
        formal_expansion_checks=expansion_controls(),
        production=rows,references=references,dynamic_comparison=dynamic,duration_controls=durations,rate_controls=rate_control,
        initial_profile_controls=initial_controls,
        force_audit=force_rows,force_audit_driver_sha256=audit['audit_driver_sha256'],
        reference_checks=checks,zero_sigma_controls=zeros,default_compatibility=compatibility,
        raw_manifest=manifest,
        diagnostic_corrections=['All production clock diagnostics recomputed from (step-start)*h^2; initial pure-pilot stored measurements used step-1 despite start=0. Original receipts retained.',
            'Analytic symbol uses cutoff |grad chi|>1e-9, namely 2|grad phi|>1e-9; earlier exploratory spectra omitted the factor two.',
            'hybrid-slip1-H48 original end-time source hash differed from frozen launch source after a docstring/scope edit. Original retained and excluded; canonical H48 rerun with immutable launch hashes is used.'],
        exclusions='Primitive-reference failures, early same-cell slip pilots, no-slip exploratory FV references and superseded H48 hybrid are retained as raw evidence but are not counted in the publication comparisons.',
        scope='Equal-density/equal-tau1 D2Q9 neutral flat walls; explicit hydraulic slip is numerical and is not molecular calibration. Spectrum predictions verify discrete phase kinetics; the hybrid shares that operator. Three-grid FV/LBM mismatch remains and neither independent full dynamic validation nor a universal contact-line law is established.')
    # Publication summary paths are relative to the retained raw root.
    def relative(value):
        if isinstance(value,dict):return {k:relative(v) for k,v in value.items()}
        if isinstance(value,list):return [relative(v) for v in value]
        if isinstance(value,str) and value.startswith(str(root)+'/'):return value[len(str(root))+1:]
        return value
    return relative(data)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();result=collect(Path(args.root).resolve());path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(production_executions=result['production_executions'],waves=len(result['waves']),
        reference_and_ablation_runs=len(result['references']),raw_files=len(result['raw_manifest']),
        max_relative_spectrum_difference=max(abs(r['relative_rate_difference']) for r in result['waves']))))


if __name__=='__main__':main()
