#!/usr/bin/env python3
"""Plot measured transport biases and independent continuum predictions."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from interfacial_transport import diffuse_reference, sharp_reference, viscosity
from benchmarks import sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', default=str(Path(__file__).resolve().parents[1]/
                        'research/drlab-publication/paper/data/interfacial-transport.json'))
    parser.add_argument('--paper', required=True)
    parser.add_argument('--sensitivity', default=str(Path(__file__).resolve().parents[1]/
                        'research/drlab-publication/paper/data/reference-sensitivity.json'))
    args = parser.parse_args()
    data = json.loads(Path(args.data).read_text())
    sensitivity = json.loads(Path(args.sensitivity).read_text())
    if sensitivity['input_sha256'] != sha(args.data):
        raise ValueError('Sensitivity data belong to a different measurement summary')
    envelopes = {r['name']: r['profile_error_envelope'] for r in sensitivity['records']}
    folder = Path(args.paper)/'figures'
    folder.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42, 'savefig.dpi': 200})
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6), constrained_layout=True)
    colors = {'arithmetic': '#187485', 'harmonic': '#a46529'}
    markers = {.35: 's', .7: 'o', .9: '^'}
    for closure, color in colors.items():
        for beta, marker in markers.items():
            selected = sorted((m for m in data['measurements'] if m['closure'] == closure
                               and m['beta'] == beta), key=lambda m: m['width_over_height'])
            if not selected:
                continue
            x = [m['width_over_height'] for m in selected]
            axes[0].scatter(x, [100*m['measured_flux_shift'] for m in selected], color=color,
                            marker=marker, s=35, label=f'{closure}, beta={beta:g}')
            axes[1].loglog(x, [100*m['measured_diffuse_profile_relative_error'] for m in selected],
                          marker+'-', color=color, ms=5, label=f'{closure}, beta={beta:g}')
            values=np.array([100*m['measured_diffuse_profile_relative_error'] for m in selected])
            bounds=np.array([envelopes[m['name']] for m in selected])*100
            axes[1].errorbar(x,values,yerr=[values-bounds[:,0],bounds[:,1]-values],
                             fmt='none',color=color,capsize=3,alpha=.6)
        widths = np.linspace(.04, .62, 150)
        grid = np.linspace(0, 1, 32769)
        _, sharp_flux = sharp_reference(np.array([.5]), 1., (1/30, .3), 1.)
        predictions = []
        for width in widths:
            xi = width/(2*np.arctanh(.9))
            phi = (1-np.tanh((grid-.5)/xi))/2
            flux = diffuse_reference(np.array([.5]), grid, viscosity(phi, (1/30, .3),
                                                                     closure), 1.)['volume_flux']
            predictions.append(100*(flux/sharp_flux-1))
        axes[0].plot(widths, predictions, color=color, lw=1, alpha=.8)
    axes[0].axhline(0, color='.5', lw=.8, ls=':')
    axes[0].set(xlabel='Measured interface width / slit width, W/H',
                ylabel='Flux shift from sharp interface (%)')
    axes[0].legend(fontsize=7.3, loc='upper left')
    axes[1].set(xlabel='Measured interface width / slit width, W/H',
                ylabel='Discrepancy against reconstructed continuum (%)')
    axes[1].set_xticks([.06, .1, .2, .4, .6])
    axes[1].set_xticklabels(['0.06', '0.1', '0.2', '0.4', '0.6'])
    axes[1].xaxis.set_minor_formatter(NullFormatter())
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    for suffix in ('pdf', 'png'):
        fig.savefig(folder/f'interfacial-transport.{suffix}')
    plt.close(fig)


if __name__ == '__main__':
    main()
