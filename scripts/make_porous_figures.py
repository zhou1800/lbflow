#!/usr/bin/env python3
"""CPC application figures/tables from recorded results or distributed summaries."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from benchmarks import geometry_fields, sha


def read(path):
    return json.loads(Path(path).read_text())


def case_receipt(folder, prefix):
    found = [read(path) for path in (Path(folder)/'results').glob('*.json')
             if read(path)['name'].startswith(prefix)]
    if len(found) != 1:
        raise ValueError(f'Expected one receipt matching {prefix}: {len(found)}')
    return found[0]


def section(result, at=None, refinement=1):
    arrays, rho, _, _ = geometry_fields(result, at)
    phi = np.divide(arrays['r0'][..., 0], rho, out=np.zeros_like(rho), where=rho > 0)
    # Average the adjacent fine-grid planes about the same original-voxel
    # center. The fine x/y resolution is retained for display.
    z = (rho.shape[0]//refinement)//2-1
    sl = slice(z*refinement, (z+1)*refinement)
    return dict(phi0=phi[sl].mean(axis=0).tolist(),
                pore=(rho[sl].mean(axis=0) > 0).tolist(),
                original_z_voxel_center=z+.5, refinement=refinement)


def assemble(args):
    root = Path(args.root)
    refinement_root = Path(args.refinement_root or root/'rock-refinement')
    if (refinement_root/'rock-continuation.json').exists():
        continued = read(refinement_root/'rock-continuation.json')
        refinement = dict(measurements=[read(root/'refined-single-measurement.json'),
                            continued['measurement']],
                          rank_change_checks=continued['rank_change_checks'],
                          protocol=continued['protocol'])
        fine_receipt = case_receipt(refinement_root, 'rock-continued-')
    else:
        refinement = read(refinement_root/'rock-refinement.json')
        fine_receipt = case_receipt(refinement_root, 'rock-refined-two-')
    reference = read(Path(args.reference_root or root/'lbpm-reference')/'lbpm-reference.json')
    if args.refined_reference_root:
        extra = read(Path(args.refined_reference_root)/'lbpm-refined.json')
        reference['measurements'].append(extra['measurement'])
        reference['refined_dataset'] = extra['dataset']
    coarse = case_receipt(root/'rock-two', 'rock-two-theta60-F2e-06-')
    data = dict(bundle=read(root/'bundle/bundle.json'),
        entry=read(root/'entry/entry.json'), entry3d=read(root/'entry-3d/entry.json'),
        bulk_interface=read(root/'rock-planar/bulk-interface.json'),
        pressure_report=read(root/'pressure-report/pressure-report.json'),
        rock_single=read(root/'rock-single/rock-single.json'),
        rock_two=read(root/'rock-two/rock-two.json'),
        rock_contracts=read(root/'rock-contracts/rock-contracts.json'),
        refinement=refinement, reference=reference,
        performance=read(root/'performance/performance.json'),
        phase_sections=dict(initial=section(coarse, 0), coarse=section(coarse, 8000),
                            refined=section(fine_receipt, 32000, 2)))
    data['provenance'] = dict(recorded_date='2026-10-08',
        source_revision_reference='See the source manifest and release metadata.',
        analysis_scripts_sha256={name:sha(Path(__file__).with_name(name)) for name in
            ('porous_media.py', 'rock_workflow.py', 'lbpm_reference.py', 'performance.py')},
        figures='Rendered only from included measured summaries and actual phase sections.',
        scope='Finite, periodically repeated cropped rock; two-phase flux remains transient.')
    return data


def scientific(value):
    if value == 0:
        return '0'
    mantissa, exponent = f'{value:.3g}'.split('e') if 'e' in f'{value:.3g}' else (None, None)
    return f'{value:.3g}' if exponent is None else mantissa+r'\times10^{'+str(int(exponent))+'}'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', help='Collect complete raw-study results; omit to rebuild from included data')
    parser.add_argument('--refinement-root')
    parser.add_argument('--reference-root')
    parser.add_argument('--refined-reference-root')
    parser.add_argument('--data', default=str(Path(__file__).resolve().parents[1]/
        'research/drlab-publication/paper/data/porous-media.json'))
    parser.add_argument('--paper', required=True)
    args = parser.parse_args()
    data = assemble(args) if args.root else read(args.data)
    paper = Path(args.paper)
    out, tables = paper/'figures', paper/'data'
    out.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    (tables/'porous-media.json').write_text(json.dumps(data, indent=2)+'\n')
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
                        'savefig.dpi':200, 'pdf.fonttype':42})

    def save(fig, name):
        fig.savefig(out/(name+'.pdf'), bbox_inches='tight', metadata={'CreationDate':None})
        fig.savefig(out/(name+'.png'), bbox_inches='tight')
        plt.close(fig)

    colors = {'legacy':'#a15b42', 'renormalized':'#187485'}
    bundles = data['bundle']['measurements']
    entries = data['entry']['measurements']+data['entry3d']['measurements']
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 6.3), constrained_layout=True)
    for ndim, marker in ((2, 'o'), (3, 's')):
        for taus, style in (([.6, 1.4], '-'), ([1.4, .6], '--')):
            cases = sorted([m for m in bundles if m['model']=='color' and
                m['dimensions']==ndim and m['tau']==taus], key=lambda m:m['scale'])
            ax[0, 0].semilogy([m['scale'] for m in cases],
                [max(m['relative_phase_flux_errors']) for m in cases], marker+style,
                label=f'D{ndim}, tau {taus[0]}/{taus[1]}')
    ax[0, 0].set(xlabel='Voxel refinement s', ylabel='Maximum relative phase-flux error')
    ax[0, 0].legend(fontsize=8)
    for method in colors:
        for ratio, style in ((.5, '--'), (1.5, '-')):
            m = next(m for m in entries if m.get('dimensions', 2)==2 and
                m['scale']==4 and m['wall_gradient']==method and m['pressure_ratio']==ratio)
            h = m['history']
            ax[0, 1].plot([p['step']/16 for p in h],
                [(p['phase0_equivalent_length']-h[0]['phase0_equivalent_length'])/4 for p in h],
                style, color=colors[method], label=f'{method}, dp/Pc={ratio}')
    ax[0, 1].axhline(0, color='.6', lw=.8)
    ax[0, 1].set(xlabel='Time in base diffusive units', ylabel='Equivalent front displacement / s')
    ax[0, 1].legend(fontsize=8)
    for col, ratio in enumerate((.5, 1.5)):
        for method in colors:
            ms = sorted([m for m in entries if m.get('dimensions', 2)==2 and
                m['wall_gradient']==method and m['pressure_ratio']==ratio], key=lambda m:m['height'])
            ax[1, col].plot([m['height'] for m in ms],
                [100*m['apparent_pressure_relative_error'] for m in ms], 'o-',
                color=colors[method], label=method+' slit')
            m = next(m for m in entries if m.get('dimensions', 2)==3 and
                m['wall_gradient']==method and m['pressure_ratio']==ratio)
            ax[1, col].plot(m['height'], 100*m['apparent_pressure_relative_error'],
                '^', mfc='none', ms=8, color=colors[method], label=method+' square tube')
        ax[1, col].axhline(0, color='.6', lw=.8)
        ax[1, col].set(xlabel='Fluid wall spacing H (nodes)',
                      ylabel='Hydraulic Pc deviation (%)', title=f'Pressure ratio {ratio}')
    ax[1, 0].legend(fontsize=8)
    save(fig, 'pore-controls')

    fig, ax = plt.subplots(2, 3, figsize=(8.4, 6.1), constrained_layout=True)
    for index, key in enumerate(('initial', 'coarse', 'refined')):
        s = data['phase_sections'][key]
        values = np.ma.array(s['phi0'], mask=~np.array(s['pore']))
        cmap = plt.get_cmap('coolwarm').copy();cmap.set_bad('#343a40')
        im = ax[0, index].imshow(values, origin='lower', cmap=cmap, vmin=0, vmax=1,
                                 extent=(0, 48, 0, 48), interpolation='nearest')
        ax[0, index].set(xlabel='x (original voxels)', ylabel='y (original voxels)',
            title=('Initial state', 'Coarse, base t=8000', 'Subdivided, base t=8000')[index])
    fig.colorbar(im, ax=list(ax[0]), label='Phase 0 fraction', shrink=.65)
    two = data['rock_two']['measurements']
    forced = next(m for m in two if m['force']==2e-6 and m['contact_angle']==60)
    zero = next(m for m in two if m['force']==0)
    fine = next(m for m in data['refinement']['measurements'] if m['two_phase'])
    for m, label, style in ((forced, 'coarse, F > 0', '-o'), (zero, 'coarse, F = 0', '--s'),
                            (fine, 'subdivided, F > 0', ':^')):
        scale = m['dataset']['refinement'];h=m['history']
        t = [p['step']/scale**2 for p in h]
        ax[1, 0].plot(t, [p['darcy_velocity'][0]*scale*1e6 for p in h], style, ms=3, label=label)
        ax[1, 2].plot(t, [p['diffuse_pore_fraction'] for p in h], style, ms=3, label=label)
    for m, style in ((forced, '-'), (fine, '--')):
        scale=m['dataset']['refinement'];h=m['history']
        for phase, color in ((0, '#b45b4d'), (1, '#3e78aa')):
            ax[1, 1].plot([p['step']/scale**2 for p in h],
                [p['phase_advective_darcy_flux'][phase][0]*scale*1e6 for p in h], style,
                color=color, label=f'phase {phase}, s={scale}')
    ax[1, 0].set(xlabel='Time in base diffusive units', ylabel='Mixture qx (base units × 10⁻⁶)')
    ax[1, 1].set(xlabel='Time in base diffusive units', ylabel='Advective phase qx (× 10⁻⁶)')
    ax[1, 2].set(xlabel='Time in base diffusive units', ylabel='Pore fraction with |chi| < 0.9')
    for a in ax[1]:a.legend(fontsize=8)
    save(fig, 'digital-rock')

    aggregates = data['performance']['aggregates']
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.4), constrained_layout=True)
    ranks=np.array([m['ranks'] for m in aggregates]);times=np.array([m['numerical_seconds_median'] for m in aggregates])
    errors=np.array([[m['numerical_seconds_median']-m['numerical_seconds_range'][0],
        m['numerical_seconds_range'][1]-m['numerical_seconds_median']] for m in aggregates]).T
    ax[0].errorbar(ranks, times, yerr=errors, fmt='o-', capsize=3, label='numerical work')
    ax[0].plot(ranks, [m['fields_diagnostics_seconds_median'] for m in aggregates], 's--', label='fields / diagnostics')
    ax[0].plot(ranks, [m['checkpoint_seconds_median'] for m in aggregates], '^--', label='checkpoint')
    ax[0].set(xlabel='MPI ranks', ylabel='Seconds', yscale='log', xticks=ranks);ax[0].legend(fontsize=8)
    ax[1].plot(ranks, times[0]/times, 'o-', label='measured numerical speedup')
    ax[1].plot(ranks, ranks/ranks[0], ':', color='.5', label='linear reference')
    ax[1].set(xlabel='MPI ranks', ylabel='Speedup', xticks=ranks);ax[1].legend(fontsize=8)
    save(fig, 'workstation-performance')

    single = data['rock_single']['measurements']
    baseline = next(m for m in single if m['model']=='scsp' and m['tau']==[1.] and m['force']==2e-6)
    comparisons = []
    for ref in data['reference']['measurements']:
        if ref['kind']=='rock':
            ours = next(m for m in single if m['model']=='scsp' and m['tau']==[ref['tau']] and m['force']==ref['force'])
        elif ref['kind']=='rock-refined':
            ours = next(m for m in data['refinement']['measurements'] if not m['two_phase'])
        else:continue
        comparisons.append(dict(name=ref['name'], refinement=ours['dataset']['refinement'],
            tau=ref['tau'], force=ref['force'], LBFlow_K_m2=ours['K_m2'], LBPM_K_m2=ref['K_m2'],
            relative_difference=ours['K_m2']/ref['K_m2']-1))
    (tables/'permeability-comparison.json').write_text(json.dumps(comparisons, indent=2)+'\n')
    ks = [m['K_m2'] for m in single if m['model']=='scsp' and m['tau']==[1.]]
    fine_k = next(m['K_m2'] for m in data['refinement']['measurements'] if not m['two_phase'])
    maximum_mass = max(m['max_relative_mass_drift']['phase_masses'] for m in two+[fine])
    values = dict(BundleFluxError=max(max(m['relative_phase_flux_errors']) for m in bundles if m['model']=='color'),
        RockPorosity=baseline['dataset']['porosity'], RockK=baseline['K_m2'],
        RockReferenceDifference=100*max(abs(m['relative_difference']) for m in comparisons),
        RockForceSpread=(max(ks)-min(ks))/baseline['K_m2'],
        RefinedKChange=100*abs(fine_k/baseline['K_m2']-1), RockMassDrift=maximum_mass,
        RockLateColor=forced['late_field_changes']['color_rms'],
        PerfFluidNodes=int(data['performance']['measurements'][0]['timing']['fluid_nodes']),
        PerfSteps=int(data['performance']['measurements'][0]['timing']['completed_steps']))
    macros=['% Generated from recorded measurements by make_porous_figures.py.']
    for key, value in values.items():
        rendered = str(value) if key in ('PerfFluidNodes', 'PerfSteps') else scientific(value)
        macros.append("\\newcommand{\\"+key+"}{\\ensuremath{"+rendered+"}}")
    (tables/'cpc-numbers.tex').write_text('\n'.join(macros)+'\n')
    lines=[r'\begin{center}\small',r'\begin{tabular}{rrrrrr}',r'\toprule',
        r'Ranks & Numerical (s) & Fluid MLUPS & Process (s) & Checkpoint (s) & Peak/rank (MiB)\\',r'\midrule']
    for m in aggregates:
        lines.append(f"{m['ranks']} & {m['numerical_seconds_median']:.3f} & {m['fluid_MLUPS_median']:.3f} & "
            f"{m['process_seconds_median']:.3f} & {m['checkpoint_seconds_median']:.3f} & {m['max_rank_rss_mib']:.1f}"+r'\\')
    lines += [r'\bottomrule',r'\end{tabular}',r'\end{center}']
    (tables/'performance-table.tex').write_text('\n'.join(lines)+'\n')
    print('Generated three CPC application figures, recorded summaries, macros and timing table.')


if __name__ == '__main__':
    main()
