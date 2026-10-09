#!/usr/bin/env python3
"""Generate the wall transport/location figure and measured TeX numbers."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def bound(value):
    if value <= 0: return '$0$'
    exponent = int(np.floor(np.log10(value)))
    mantissa = int(np.ceil(value/10.**exponent))
    return f'${mantissa}\\times10^{{{exponent}}}$'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True)
    parser.add_argument('--paper', required=True)
    args = parser.parse_args(); data = json.loads(Path(args.data).read_text()); paper = Path(args.paper)
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42, 'savefig.dpi': 200})
    fig, axes = plt.subplots(2, 2, figsize=(7.7, 6.3), constrained_layout=True)
    for item, color, label in zip(data['independent']['neutral_scalar_comparison'],
            ['#a26336', '#416592', '#248875'], ['Local bounce-back', 'Stationary balance', 'Relative-phase flux']):
        history = item['history']
        axes[0, 0].loglog([row['step'] for row in history],
                         [max(row['max_bulk_transient_error'], 1e-17) for row in history],
                         'o-', ms=3, color=color, label=label)
    axes[0, 0].set(xlabel='Step (phase-only calculation)', ylabel='Maximum bulk-recurrence discrepancy')
    axes[0, 0].set_title(r'$\beta=0.7$, $H=8$, $L=64$', fontsize=9)
    axes[0, 0].legend(fontsize=8)
    for mode, marker, label in [('bounce', 'o', 'Local bounce-back'), ('flux', 's', 'Relative-phase flux')]:
        rows = sorted([m for m in data['neutral'] if m['wall_gradient'] == 'ghost'
                       and m['phase_boundary'] == mode and m['sigma'] > 0], key=lambda m: m['beta'])
        axes[0, 1].semilogy([m['beta'] for m in rows], [m['max_speed'] for m in rows],
                           marker+'-', color='#a26336' if mode=='bounce' else '#248875', label=label)
    axes[0, 1].set(xlabel=r'Segregation strength, $\beta$', ylabel='Neutral raw maximum speed')
    axes[0, 1].legend(fontsize=8)
    axes[0, 1].set_title(r'Ghost color, $H=16$, $t=8000$', fontsize=9)
    colors = {'renormalized': '#a26336', 'halfway': '#416592', 'ghost': '#248875'}
    labels = {'renormalized': 'Center rotation', 'halfway': 'Halfway cosine', 'ghost': 'Ghost color'}
    chosen = {}
    for gradient in colors:
        for mode in ('bounce', 'flux'):
            rows = sorted([m for m in data['long_holds'] if m['wall_gradient'] == gradient
                           and m['phase_boundary'] == mode and m['angle'] == 60 and m['beta'] == .7],
                          key=lambda m: m['height'])
            if [m['height'] for m in rows] != [16, 32, 48]:
                raise ValueError('Pressure figure needs the complete matched three-width hold sequence')
            if any(m['steps'] != 32000*(m['height']/16)**2 for m in rows):
                raise AssertionError('Four-duration protocol differs')
            chosen[(gradient, mode)] = rows
            style = 'o-' if mode == 'bounce' else 's--'
            label = labels[gradient]+(' / BB' if mode == 'bounce' else ' / flux')
            axes[1, 0].plot([m['height'] for m in rows], [100*m['relative_pressure_error'] for m in rows],
                            style, color=colors[gradient], ms=4, label=label)
            axes[1, 1].plot([m['height'] for m in rows],
                            [m['history'][-1]['max_speed']/6/m['sigma'] for m in rows],
                            style, color=colors[gradient], ms=4, label=label)
    heights = np.linspace(16, 48, 100)
    axes[1, 0].plot(heights, 100/(heights-1), ':', color='.5', lw=1, label=r'$100/(H-1)$ estimate (c)')
    axes[1, 0].axhline(0, color='.5', lw=.6)
    axes[1, 0].set(xlabel=r'Slit width, $H$', ylabel='Signed pressure discrepancy (%)')
    handles, legend_labels = axes[1, 0].get_legend_handles_labels()
    axes[1, 1].legend(handles, legend_labels, fontsize=6.5, loc='upper right')
    axes[1, 1].set(xlabel=r'Slit width, $H$', ylabel=r'Nonneutral residual $\mu U_{max}/\sigma$')
    for axis, letter in zip(axes.flat, 'abcd'):
        axis.spines[['top', 'right']].set_visible(False)
        axis.text(-.15, 1.02, f'({letter})', transform=axis.transAxes)
    (paper/'figures').mkdir(parents=True, exist_ok=True)
    (paper/'data').mkdir(parents=True, exist_ok=True)
    for extension in ('png', 'pdf'): fig.savefig(paper/'figures'/f'wall-flux.{extension}')
    plt.close(fig)
    checks = data['contracts']
    neutral_speed = max([m['max_speed'] for m in data['neutral'] if m['phase_boundary'] == 'flux' and m['sigma'] > 0]
                        +[m['max_speed_every_step'] for m in checks['surface_transients']])
    values = dict(FluxTransientBound=bound(max(m['max_transient_error'] for m in checks['neutral_transients'])),
                  FluxNeutralSpeedBound=bound(neutral_speed))
    for item, suffix in zip(chosen[('ghost', 'flux')], ('Sixteen', 'ThirtyTwo', 'FortyEight')):
        values['FluxGhostPressure'+suffix] = f'{100*item["relative_pressure_error"]:.3f}'
    (paper/'data/wall-flux-numbers.tex').write_text('% Generated from wall-flux.json; measured bounds rounded upward.\n'+
        ''.join('\\newcommand{\\'+key+'}{'+value+'}\n' for key, value in values.items()))
    table = [r'\begin{table}[htbp]\centering\small',
             r'\caption{Ghost-color H16 nonneutral controls after four initial durations. Pressure uses the specified bulk sample; residual speed is unfiltered. Other sampling and contour-exclusion results remain in the data.}',
             r'\begin{tabular}{rrrrr}\toprule',
             r'$\beta$ & $\theta$ & Phase reflection & Pressure error (\%) & $\mu U_{max}/\sigma$ \\ \midrule']
    rows = sorted([m for m in data['long_holds'] if m['height'] == 16 and m['wall_gradient'] == 'ghost'],
                  key=lambda m: (m['beta'], m['angle'], m['phase_boundary']))
    for m in rows:
        mode = 'BB' if m['phase_boundary'] == 'bounce' else 'Flux'
        table.append(f'{m["beta"]:g} & {m["angle"]:g} & {mode} & {100*m["relative_pressure_error"]:+.3f} & '+
                     f'{m["history"][-1]["max_speed"]/6/m["sigma"]:.3g}'+r' \\')
    table += [r'\bottomrule\end{tabular}\end{table}']
    (paper/'data/wall-flux-table.tex').write_text('\n'.join(table)+'\n')
    print(json.dumps(values))


if __name__ == '__main__': main()
