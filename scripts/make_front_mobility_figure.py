#!/usr/bin/env python3
"""Rebuild the fixed-width/resistance scientific figure and numeric macros."""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from make_core_reservoir_figures import decorate, save, scientific


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',required=True);parser.add_argument('--paper',required=True)
    args=parser.parse_args();data=json.loads(Path(args.data).read_text());paper=Path(args.paper)
    for sub in ('data','figures'):(paper/sub).mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'savefig.dpi':200})
    fig,axes=plt.subplots(2,2,figsize=(7.7,6.3),constrained_layout=True)
    colors={16:'#248875',32:'#416592',48:'#a26336'}
    for row in data['calibration']:
        h=row['height'];profile=np.array(row['profile']);i=np.arange(len(profile))
        # The ascending periodic interface lies halfway between the last/first nodes.
        x=np.where(i<len(profile)/2,i+.5,i-len(profile)+.5)*16/h
        mask=abs(x)<=6;order=np.argsort(x[mask])
        axes[0,0].plot(x[mask][order],profile[mask][order],'o-',ms=3,color=colors[h],label=f'H={h}')
    axes[0,0].set(xlabel='Distance in H16 length units',ylabel='Planar color indicator')
    axes[0,0].legend(fontsize=8)
    fit=[r for r in data['response_fits'] if r['held'] and r['offset']==.02 and r['length_factor']==4]
    fit.sort(key=lambda r:r['height']);x=[r['height'] for r in fit]
    axes[0,1].plot(x,[100*r['relative_deficit'] for r in fit],'o-',color='#416592')
    axes[0,1].set(xlabel='Slit width, H',ylabel='Front mobility deficit (%)',title='Fixed measured W09/H')
    for length,marker,color in [(4,'o','#248875'),(8,'s','#a26336')]:
        rows=sorted([r for r in data['response_fits'] if r['held'] and r['height']==16 and r['length_factor']==length],key=lambda r:r['offset'])
        axes[1,0].plot([r['offset'] for r in rows],[r['extra_resistance'] for r in rows],marker+'-',color=color,label=f'L={length}H')
    axes[1,0].set(xlabel=r'Pressure halfwidth / $p_Y$',ylabel=r'Excess resistance, $R_e$ (lu)')
    axes[1,0].legend(fontsize=8);axes[1,0].ticklabel_format(axis='y',useOffset=False)
    axes[1,0].set_ylim(.06229,.06235)
    for h,length,marker in [(16,4,'o'),(16,8,'s'),(32,4,'^'),(48,4,'D')]:
        row=next(r for r in data['response_fits'] if r['held'] and r['height']==h and r['length_factor']==length and r['offset']==.02)
        axes[1,1].plot(row['mean_speed_extra_resistance'],row['force_work_extra_resistance'],marker,color=colors[h],ms=7,label=f'H={h}, L={length}H')
    limits=axes[1,1].get_xlim();axes[1,1].plot(limits,limits,'--',color='.5',lw=.8)
    axes[1,1].set(xlabel='Pressure / flow excess resistance (lu)',ylabel='Capillary work excess resistance (lu)')
    axes[1,1].legend(fontsize=7);decorate(axes);save(fig,paper,'front-resistance')
    longest=[r for r in data['response_fits'] if r['held']]
    four=next(r for r in longest if r['height']==16 and r['length_factor']==4 and r['offset']==.02)
    eight=next(r for r in longest if r['height']==16 and r['length_factor']==8 and r['offset']==.02)
    vals=dict(FrontWidthRelativeError=scientific(max(abs(r['relative_width_error']) for r in data['calibration'])),
        FrontLengthResistancePercent=f'{100*abs(eight["extra_resistance"]/four["extra_resistance"]-1):.4f}',
        FrontForceResistancePercent=f'{math.ceil(100000*max(abs(r["relative_force_resistance_residual"]) for r in longest))/1000:.3f}',
        FrontDurationMobilityPercent=f'{100*max(abs(r["relative_duration_mobility_change"]) for r in longest):.4f}',
        FrontMaxBudget=scientific(max(r['budget']['max_relative_budget_residual'] for r in data['fronts']+data['channels'])),
        FrontMaxForceResidual=scientific(max(abs(r['prediction_minus_raw']) for r in data['reciprocal'])),
        FrontFineDeficit=f'{100*fit[-1]["relative_deficit"]:.2f}')
    (paper/'data/front-mobility-numbers.tex').write_text('% Generated from verified front-mobility.json.\n'+
        ''.join('\\newcommand{\\'+name+'}{'+value+'}\n' for name,value in vals.items()))
    table=[r'\begin{table}[htbp]\centering\small',
        r'\caption{Fixed-width front response at the doubled duration and pressure halfwidth $0.02p_Y$. Mobilities and $R_e$ use native lattice units. The force-work residual compares resistance inferred from instantaneous physical mean flow with resistance from independently weighted force.}',
        r'\begin{tabular}{rrrrrr}\toprule',
        r'$H$ & $L/H$ & $G_0$ & $M_V$ & $R_e$ & Force residual (\%) \\ \midrule']
    for row in sorted([r for r in longest if r['offset']==.02],key=lambda r:(r['height'],r['length_factor'])):
        table.append(f'{row["height"]} & {row["length_factor"]} & {row["pure_mobility"]:.6f} & {row["front_mobility"]:.6f} & {row["extra_resistance"]:.6f} & {100*row["relative_force_resistance_residual"]:+.4f}'+r' \\')
    table.extend([r'\bottomrule\end{tabular}\end{table}'])
    (paper/'data/front-mobility-table.tex').write_text('\n'.join(table)+'\n');print(json.dumps(vals))


if __name__=='__main__':main()
