#!/usr/bin/env python3
"""Reanalyze frozen major-revision fields and retain unfavorable controls."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True

import numpy as np

from benchmarks import geometry_fields,sha
from collect_capillary_controls import verify
from collect_phase_kinetics import receipt_at_root,reference_measure
from periodic_contact import measure
from phase_wave_spectrum import spectrum,matrix
from stokes_phase_reference import projection_rest
from full_force_fields import collect as collect_force_fields


def native_campaign(root,campaign):
    rows=[]
    for path in sorted((root/campaign/'results').glob('*.json')):
        r=receipt_at_root(json.loads(path.read_text()),root);verify(r)
        p=r['parameters'];folder=Path(r['output_files'][-1]).parent
        if p.get('seed_phase_sha256') and sha(p['seed_phase_file'])!=p['seed_phase_sha256']:
            raise AssertionError('Initial phase source changed')
        for key,suffix in (('checkpoint_sha256','.rst'),('checkpoint_metadata_sha256','.rst.meta')):
            if r.get(key) and sha(folder/f'field.{r["start"]:010d}{suffix}')!=r[key]:
                raise AssertionError('Outgoing donor or its time-stage metadata changed')
        row=measure(r);row['campaign']=campaign;row['receipt_sha256']=sha(path)
        _,_,chi,_=geometry_fields(r);phi=(1+chi[0,1:-1])/2;h=p['dx']
        shape=-h*phi[:,:phi.shape[1]//2].sum(axis=1);shape-=shape.mean()
        row['phase_volume_shape']=shape.tolist()
        row['max_relative_phase_mass_change']=max(abs(t['phase_mass']/row['history'][0]['phase_mass']-1) for t in row['history'])
        rows.append(row)
    if not rows:raise FileNotFoundError(campaign)
    return rows


def manufactured(root):
    summary=json.loads((root/'summary.json').read_text())
    for name,digest in summary['source_sha256'].items():
        if sha(root/'source'/name)!=digest:raise AssertionError('Manufactured frozen source changed')
    spec=importlib.util.spec_from_file_location('frozen_manufactured',root/'source/coupled_reference_checks.py')
    original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
    fresh=[]
    for folder in sorted(p for p in root.iterdir() if p.is_dir() and p.name!='source'):
        recorded=json.loads((folder/'result.json').read_text());p=recorded['parameters'];history=[];masses=[];ranges=[]
        for file in sorted(folder.glob('field.*.npz')):
            with np.load(file) as f:
                phi,u,v=f['phi'],f['u'],f['v'];h=float(f['h']);t=float(f['time'])
                if not all(np.isfinite(a).all() for a in (phi,u,v)):raise AssertionError('Nonfinite manufactured field')
                xc=(np.arange(p['nx'])+.5)[None]*h;xf=np.arange(p['nx'])[None]*h
                yc=(np.arange(p['ny'])+.5)[:,None]*h;yf=np.arange(p['ny']+1)[:,None]*h
                target=original.target(xc,yc,t,b=p['b']);uf=original.flow(xf,yc,p['slip'])['u'];vf=original.flow(xc,yf,p['slip'])['v']
                if not np.array_equal(f['exact_phi'],target['phi']) or not np.array_equal(f['exact_u'],uf) or not np.array_equal(f['exact_v'],vf):
                    raise AssertionError('Analytic targets do not reproduce saved fields')
                row=dict(time=t,relative_phase_error=float(np.linalg.norm(phi-target['phi'])/np.linalg.norm(target['phi'])),
                    interface_displacement_error=float(np.linalg.norm(phi-target['phi'])/np.linalg.norm(target['interface_derivative'])),
                    relative_velocity_error=float(np.sqrt(((u-uf)**2).sum()+((v-vf)**2).sum())/np.sqrt((uf*uf).sum()+(vf*vf).sum())),
                    max_divergence=float(abs((np.roll(u,-1,axis=1)-u)/h+np.diff(v,axis=0)/h).max()))
                history.append(row)
                masses.append(float(phi.sum()*h*h));ranges.append([float(phi.min()),float(phi.max())])
        if len(history)!=len(recorded['history']):raise AssertionError('Manufactured snapshots missing')
        for row,old in zip(history,recorded['history']):
            if max(abs(row[k]-old[k]) for k in row)>2e-15:raise AssertionError('Manufactured measurements changed')
        if ranges!=[r['phase_range'] for r in recorded['history']]:raise AssertionError('Manufactured phase ranges changed')
        # Reintegrate the independently prescribed continuous source. No
        # evolving numerical phase or discretized flux is used in this budget.
        cell_flow=original.flow(xc,yc,p['slip']);steps=round(p['end']/p['dt']);budget=0.;budgets={0:0.}
        def integral(clock):
            exact=original.target(xc,yc,clock,b=p['b'])
            source=exact['dt']+cell_flow['u']*exact['dx']+cell_flow['v']*exact['dy']-p['diffusion']*exact['phase_divergence']
            return float(source.sum()*h*h)
        previous=integral(0.)
        for step in range(1,steps+1):
            following=integral(step*p['dt']);budget+=p['dt']*(previous+following)/2;previous=following
            if step in (steps//2,steps):budgets[step]=budget
        residuals=[mass-masses[0]-budgets[round(row['time']/p['dt'])] for mass,row in zip(masses,history)]
        if max(abs(a-r['mass_budget_residual']) for a,r in zip(residuals,recorded['history']))>2e-12:
            raise AssertionError('Manufactured continuous-source mass budget changed')
        fresh.append(dict(name=folder.name,parameters=p,history=history,
            recomputed_mass_budget_residuals=residuals,phase_ranges=ranges))
    orders=[]
    for slip in (1.,4.):
        rows=sorted([r for r in fresh if r['parameters']['slip']==slip and r['parameters']['dt_factor']==.125],key=lambda r:r['parameters']['ny'])
        for key in ('relative_phase_error','interface_displacement_error','relative_velocity_error'):
            errors=[r['history'][-1][key] for r in rows]
            orders.append(dict(slip=slip,observable=key,errors=errors,
                orders=[float(np.log(a/b)/np.log(right['parameters']['ny']/left['parameters']['ny'])) for a,b,left,right in zip(errors,errors[1:],rows,rows[1:])]))
    if orders!=summary['spatial_orders']:raise AssertionError('Manufactured convergence orders changed')
    return dict(runs=fresh,spatial_orders=orders,source_sha256=summary['source_sha256'],
        scope='Continuous-source coupled consistency, with C2 exact pure plateaus; not unforced physical contact-line validation.')


def mechanism(root):
    recorded=json.loads((root/'summary.json').read_text())
    for name,digest in recorded['source_sha256'].items():
        if sha(root/'source'/name)!=digest:raise AssertionError('Reduced-model frozen source changed')
    env=dict(os.environ,PYTHONPATH=str(root/'source')+os.pathsep+os.environ.get('PYTHONPATH',''),OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    code='import json; from phase_mechanism_controls import calculate; print(json.dumps(calculate()))'
    result=json.loads(subprocess.check_output([sys.executable,'-c',code],cwd=root/'source',env=env,text=True))
    if result!={k:v for k,v in recorded.items() if k!='source_sha256'}:
        raise AssertionError('Independent reduced calculation differs from frozen receipt')
    result['source_sha256']=recorded['source_sha256']
    return result


def collect(root):
    before=root/'bound-inputs/phase-before-revision';old=json.loads((root/'bound-inputs/phase-kinetics-before-revision.json').read_text())
    old_native=[r for c in ('production-slip1-forward-H16','production-slip1-forward-fine','production-slip4','waves-pilot','waves-fine','waves-mode2') for r in native_campaign(before,c)]
    native=[r for c in ('lbm-sigma-controls','lbm-mode3','lbm-width6','lbm-low-sigma-hold') for r in native_campaign(root/'experiments',c)]
    references=[reference_measure(p) for p in sorted((root/'experiments').glob('fv-*')) if (p/'result.json').exists()]
    old_refs=[reference_measure(before/name) for name in ('reference-slip1-H16','reference-slip1-H32','reference-slip1-H48','reference-slip4-H16','reference-slip4-H32','reference-slip4-H48')]
    waves=[];spectra={}
    for row in old_native+native:
        if not row['parameters']['physical_amplitude']:continue
        p=row['parameters'];key=(p['height'],p['transverse_mode'],p.get('physical_width',4.765235021846117))
        if key not in spectra:spectra[key]=spectrum(key[0],key[1],width=key[2])
        checks=[]
        for fraction in (.4,.5,.6):
            late=row['history'][int(len(row['history'])*fraction):]
            fit=np.polyfit([t['time'] for t in late],np.log([abs(t['projected_displacement']) for t in late]),1)
            checks.append(dict(window_start_fraction=fraction,rate=float(-fit[0])))
        predicted=spectra[key]['physical_relaxation_rate'];rate=checks[1]['rate']
        waves.append(dict(name=row['name'],height=key[0],mode=key[1],width=key[2],amplitude=p['physical_amplitude'],
            measured_rate=rate,predicted_rate=predicted,relative_rate_difference=rate/predicted-1,
            fit_window_checks=checks,history=row['history']))
    cutoffs=[]
    for (height,mode,width),s in spectra.items():
        if mode!=1:continue
        profile=np.array(s['profile']);derivative=(np.roll(profile,-1)-np.roll(profile,1))/2
        for cutoff in (1e-10,1e-9,1e-8):
            residual=float(abs(projection_rest(np.broadcast_to(profile,(4,len(profile))),s['beta'],cutoff=cutoff)-profile[None]).max())
            row=dict(height=height,width=width,cutoff=cutoff,profile_residual=residual,equilibrium_eligible=residual<=3e-14)
            if row['equilibrium_eligible']:
                values,vectors=np.linalg.eig(matrix(profile,s['beta'],s['wavenumber_lattice'],cutoff=cutoff))
                overlap=abs(vectors.conj().T@derivative)/np.linalg.norm(vectors,axis=0)/np.linalg.norm(derivative)
                chosen=values[int(np.argmax(overlap))]
                if abs(chosen.imag)>1e-12 or not 0<chosen.real<1:raise AssertionError('Cutoff control is not a stable translation mode')
                row['rate']=float(-np.log(chosen.real)/(16/height)**2)
            cutoffs.append(row)
    response=[];shape_comparisons=[]
    model=mechanism(root/'experiments/mechanism-controls-canonical')
    for prediction in model['reduced_predictions']:
        H=prediction['height'];slip=prediction['slip'];sigma=prediction['sigma'];a=prediction['strength']
        candidates=[r for r in references+old_refs if r['parameters']['ny']==H and r['parameters']['slip']==slip
            and r['parameters']['sigma']==sigma and r['parameters']['force']==1e-7 and r['parameters'].get('transverse_relaxation',0.)==a
            and r['parameters'].get('transverse_localization','uniform')=='uniform'
            and 'fixed-b' not in r['name']]
        if candidates:
            ref=max(candidates,key=lambda r:r['parameters']['end'])
            folder=root/'experiments'/ref['name'] if ref in references else before/ref['name']
            with np.load(sorted(folder.glob('field.*.npz'))[-1]) as field:
                phi=field['phi'];h=float(field['h'])
                observed=-h*phi[:,:phi.shape[1]//2].sum(axis=1);observed-=observed.mean()
            from stokes_phase_reference import logistic_bridge,match_logistic_width
            profile=logistic_bridge(4*H,h,match_logistic_width(4.765235021846117,h))
            weight=h*(2*match_logistic_width(4.765235021846117,h)*profile*(1-profile))[:2*H].sum()
            expected=np.array(prediction['shape'])*weight
            response.append(dict(height=H,slip=slip,sigma=sigma,strength=a,reference_name=ref['name'],
                predicted_mean_speed=prediction['predicted_mean_speed'],measured_mean_speed=ref['late_mean_speed'],
                relative_prediction_error=prediction['predicted_mean_speed']/ref['late_mean_speed']-1,
                predicted_phase_speed=prediction['predicted_phase_speed'],measured_phase_speed=ref['late_centroid_speed'],
                relative_volume_shape_error=float(np.linalg.norm(expected-observed)/np.linalg.norm(observed)),
                reference_end=ref['parameters']['end']))
        lbm=[r for r in old_native+native if r['parameters']['height']==H and r['parameters']['physical_force']==1e-7
            and r['parameters']['physical_sigma']==sigma and r['parameters'].get('physical_slip',0)==slip]
        if lbm and a==1:
            row=max(lbm,key=lambda r:r['history'][-1]['time'])
            from stokes_phase_reference import logistic_bridge,match_logistic_width
            b=match_logistic_width(4.765235021846117,16/H);profile=logistic_bridge(4*H,16/H,b)
            weight=(16/H)*(2*b*profile*(1-profile))[:2*H].sum()
            expected=np.array(prediction['shape'])*weight;observed=np.array(row['phase_volume_shape'])
            shape_comparisons.append(dict(height=H,slip=slip,sigma=sigma,
                relative_reduced_lbm_mean_error=prediction['predicted_mean_speed']/row['late_mean_speed']-1,
                relative_volume_shape_error=float(np.linalg.norm(expected-observed)/np.linalg.norm(observed)),
                predicted_volume_shape=expected.tolist(),measured_volume_shape=observed.tolist()))
    force=collect_force_fields(before)
    return dict(protocol='Referee-driven major revision: source-backed spectral/causal controls, continuous-source coupled verification and a predictive low-force resistance reduction.',
        baseline_commit='1ed179a9092819dff8a3b17e934dbdac7e99fe91',baseline_tag='cpc-pre-major-revision-20261009',
        new_native_executions=len(native),new_reference_executions=len(references),full_solver_waves=len(waves),
        native=native,previous_native=old_native,references=references,previous_references=old_refs,waves=waves,
        cutoff_controls=cutoffs,manufactured=manufactured(root/'experiments/manufactured-canonical'),
        mechanism_controls=model,reduced_reference_comparison=response,reduced_lbm_comparison=shape_comparisons,
        full_force_fields=force,previous_dataset_sha256=sha(root/'bound-inputs/phase-kinetics-before-revision.json'),
        exclusions='Primitive failures retained in the baseline archive; three unavailable-width CLI pilots and the pre-C2 manufactured pilot remain separate from publication measurements.',
        scope='Neutral equal-viscosity D2Q9 and prescribed hydraulic slip; numerical regularization/resistance diagnosis, not a molecular or universal contact-line law or quantitative curved-rock wetting.')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',required=True);parser.add_argument('--output',required=True)
    parser.add_argument('--verify-manifest',help='Verify every retained raw file before reanalysis; paths are relative to --root.')
    args=parser.parse_args();root=Path(args.root).resolve()
    if args.verify_manifest:
        manifest=json.loads(Path(args.verify_manifest).read_text())
        for name,entry in manifest['files'].items():
            path=root/name
            if not path.is_file() or path.stat().st_size!=entry['bytes'] or sha(path)!=entry['sha256']:
                raise AssertionError('Raw evidence changed: '+name)
    result=collect(root)
    def relative(value):
        if isinstance(value,dict):return {k:relative(v) for k,v in value.items()}
        if isinstance(value,list):return [relative(v) for v in value]
        if isinstance(value,str) and value.startswith(str(root)+'/'):return value[len(str(root))+1:]
        return value
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(relative(result),indent=2)+'\n')
    print(json.dumps(dict(new_native_runs=result['new_native_executions'],new_reference_runs=result['new_reference_executions'],waves=result['full_solver_waves'],
        max_relative_wave_rate_difference=max(abs(w['relative_rate_difference']) for w in result['waves']))))


if __name__=='__main__':main()
