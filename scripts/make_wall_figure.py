#!/usr/bin/env python3
"""Regenerate the wall-equilibrium figure and numerical bounds from data."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',required=True);p.add_argument('--paper',required=True)
    args=p.parse_args();data=json.loads(Path(args.data).read_text());paper=Path(args.paper)
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'savefig.dpi':200})
    fig,axes=plt.subplots(2,2,figsize=(7.6,6.2),constrained_layout=True)
    colors={'legacy':'#7b4e2d','tangential':'#187485'}
    for mode in colors:
        rows=[m for m in data['neutral'] if m['height']==16 and m['axis']=='x' and m['collision']=='MRT'
              and m['sigma']==0 and m['wall_recoloring']==mode]
        rows=sorted(rows,key=lambda m:m['beta'])
        axes[0,0].semilogy([m['beta'] for m in rows],[max(m['max_wall_color_difference'],1e-16) for m in rows],
                          'o-',color=colors[mode],label='Original recoloring' if mode=='legacy' else 'Balanced neutral wall')
    axes[0,0].set(xlabel=r'Segregation strength, $\beta$',ylabel='Zero-force wall color distortion')
    axes[0,0].legend(fontsize=8)
    rows=[m for m in data['staggered'] if m['parameters'].get('wall_recoloring')=='tangential'
          and m['parameters'].get('curvature_stencil','normal')=='normal'
          and m['parameters']['beta']==.7 and m['parameters']['size'][0]==4*m['parameters']['height']]
    rows=sorted(rows,key=lambda m:m['parameters']['height'])
    for key,label,marker in [('raw_speed','Raw corrected CSF','o'),
                             ('adjacent_pair_average_speed','Adjacent-step average (diagnostic)','s')]:
        axes[0,1].semilogy([m['parameters']['height'] for m in rows],[max(m[key],1e-17) for m in rows],
                          marker+'-',label=label)
    baseline=[m for m in data['staggered'] if m['parameters'].get('wall_recoloring','legacy')=='legacy'
              and m['parameters'].get('capillary_model')=='csf']
    baseline=sorted(baseline,key=lambda m:m['parameters']['height'])
    axes[0,1].semilogy([m['parameters']['height'] for m in baseline],[m['raw_speed'] for m in baseline],
                      '^--',color=colors['legacy'],label='Uncorrected CSF')
    fixed_max=max(m['max_speed'] for m in data['fixed_point_checks']['fixed_points'])
    axes[0,1].axhline(fixed_max,color='.5',ls=':',lw=.8,label='Manufactured equilibrium maximum')
    axes[0,1].set(xlabel=r'Slit width, $H$',ylabel='Maximum speed (lattice units)')
    axes[0,1].legend(fontsize=7)
    chosen=next(m for m in rows if m['parameters']['height']==16)
    velocity=np.array(chosen['tangent_velocity_rows'])
    axes[1,0].plot(np.arange(64),velocity[0]/chosen['raw_speed'],'o-',ms=2,lw=.7,label=r'Wall-adjacent $u_x/U_{max}$')
    axes[1,0].set(xlabel='Tangential cell index',ylabel='Alternating raw velocity',ylim=(-1.2,1.2))
    axes[1,0].legend(fontsize=8)
    for form,marker in [('normal','o'),('gradient','s')]:
        rows=[m for m in data['matched_bulk'] if m['parameters']['curvature_stencil']==form]
        rows=sorted(rows,key=lambda m:m['parameters']['radius'])
        axes[1,1].plot([m['parameters']['radius'] for m in rows],
                       [100*m['relative_laplace_error'] for m in rows],marker+'-',label=form.capitalize()+' derivatives')
    axes[1,1].set(xlabel='Bulk droplet radius (cells)',ylabel='Laplace discrepancy (%)')
    axes[1,1].legend(fontsize=8)
    for ax,label in zip(axes.flat,'abcd'):
        ax.spines[['top','right']].set_visible(False);ax.text(-.15,1.02,f'({label})',transform=ax.transAxes)
    (paper/'figures').mkdir(parents=True,exist_ok=True)
    for extension in ('png','pdf'):fig.savefig(paper/'figures'/f'wall-balance.{extension}')
    plt.close(fig)
    def bound(value):
        quantum=10.**np.floor(np.log10(value))
        return f'${np.ceil(value/quantum)*quantum:.0e}$'
    zero=max(m['max_wall_color_difference'] for m in data['neutral'] if m['sigma']==0 and m['wall_recoloring']=='tangential')
    checks=data['fixed_point_checks']['fixed_points']
    values=dict(WallZeroColorBound=bound(zero),
                WallFixedColorBound=bound(max(m['color_error'] for m in checks)),
                WallFixedSpeedBound=bound(max(m['max_speed'] for m in checks)),
                WallRawSpeed=f'${chosen["raw_speed"]:.3e}$',
                WallPairSpeedBound='$1\\times10^{-14}$',WallStaggeredEnergyPercent='99.9999999')
    if chosen['adjacent_pair_average_speed']>=1e-14:raise AssertionError('Adjacent diagnostic bound changed')
    if 100*chosen['staggered_component_energy_fractions'][0][0]<=99.9999999:raise AssertionError('Staggered-energy bound changed')
    # TeX scientific notation must be mathematical, not a literal programming exponent.
    for key,value in values.items():
        if 'e-' in value:
            mantissa,exponent=value.strip('$').split('e-');values[key]=f'${mantissa}\\times10^{{-{int(exponent)}}}$'
    (paper/'data/wall-numbers.tex').write_text('% Generated from wall-consistency.json; bounds checked and rounded upward.\n'+
        ''.join('\\newcommand{\\'+key+'}{'+value+'}\n' for key,value in values.items()))
    print(json.dumps(values))


if __name__=='__main__':main()
