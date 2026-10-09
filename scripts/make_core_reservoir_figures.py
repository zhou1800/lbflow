#!/usr/bin/env python3
"""Scientific figures and numeric macros from the verified core/reservoir data."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',required=True);parser.add_argument('--paper',required=True)
    args=parser.parse_args();data=json.loads(Path(args.data).read_text());paper=Path(args.paper)
    for sub in ('data','figures'):(paper/sub).mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'savefig.dpi':200})
    matched=next(r['beta'] for r in data['calibration']['models'] if r['recoloring_model']=='projection' and r['beta']!=.7)
    groups=[('cosine','legacy','Cosine / legacy MRT','#a26336','o'),
            ('projection','legacy','Projection / legacy MRT','#248875','s'),
            ('projection','liu2015','Projection / published MRT','#416592','^')]
    fig,axes=plt.subplots(2,2,figsize=(7.7,6.3),constrained_layout=True)
    for model,beta,label,color,style in [('cosine',.7,'Cosine, beta=0.7','#a26336','-'),
            ('projection',.7,'Projection, beta=0.7','#416592',':'),
            ('projection',matched,'Projection, matched beta','#248875','--')]:
        row=next(r for r in data['calibration']['models'] if r['recoloring_model']==model and r['mrt_spectrum']=='legacy' and r['beta']==beta)
        profile=np.array(row['profile']);x=np.arange(len(profile))-63.5;mask=abs(x)<8
        axes[0,0].plot(x[mask],profile[mask],style,color=color,label=label)
    axes[0,0].set(xlabel='Distance from planar interface',ylabel='Color indicator');axes[0,0].legend(fontsize=7)
    for model,spectrum,label,color,marker in groups:
        rows=data['baseline'] if model=='cosine' else [r for r in data['closed'] if r['campaign'] in ('core-h16-holds','core-fine-holds') and r['beta']!=.7 and r['phase_boundary']=='flux' and r['mrt_spectrum']==spectrum]
        rows=sorted(rows,key=lambda r:r['height'])
        axes[0,1].plot([r['height'] for r in rows],[100*r['relative_pressure_error'] for r in rows],marker+'-',color=color,label=label)
        axes[1,0].plot([r['height'] for r in rows],[r['history'][-1]['max_speed']/6/r['sigma'] for r in rows],marker+'-',color=color,label=label)
    axes[0,1].axhline(0,color='.5',lw=.7);axes[0,1].set(xlabel='Slit width, H',ylabel='Signed pressure discrepancy (%)')
    axes[0,1].legend(fontsize=7)
    axes[1,0].set(xlabel='Slit width, H',ylabel='Nonneutral residual capillary number')
    for tau,marker,color in [(.6,'o','#416592'),(1.,'s','#248875'),(1.55,'^','#a26336')]:
        rows=sorted([r for r in data['force_channels']['measurements'] if r['mrt_spectrum']=='liu2015' and r['tau']==tau],key=lambda r:r['height'])
        axes[1,1].plot([r['height'] for r in rows],[100*r['relative_unshifted_flux_error'] for r in rows],marker,color=color,label=f'tau={tau:g}')
        hgrid=np.linspace(8,32,100)
        axes[1,1].plot(hgrid,1200*rows[0]['constant_offset']/(2*hgrid*hgrid+1),'-',color=color,lw=.8)
    axes[1,1].set(xlabel='Slit width, H',ylabel='Published MRT flux offset (%)');axes[1,1].legend(fontsize=8)
    decorate(axes)
    save(fig,paper,'published-core')
    fig,axes=plt.subplots(2,2,figsize=(7.7,6.3),constrained_layout=True)
    for h,axis in [(16,axes[0,0]),(32,axes[0,1])]:
        for model,spectrum,label,color,marker in groups:
            rows=sorted([r for r in data['fronts'] if r['height']==h and r['recoloring_model']==model and r['mrt_spectrum']==spectrum and
                r['campaign'] in ('open-bracket-holds','open-cosine-holds')],key=lambda r:r['pressure_ratio'])
            axis.plot([r['pressure_ratio'] for r in rows],[1e6*r['late_volume_front_speed'] for r in rows],marker+'-',color=color,label=label)
        axis.axhline(0,color='.5',lw=.7);axis.set(xlabel=r'Reservoir pressure, $\Delta p/p_Y$',ylabel=r'Front speed ($10^{-6}$ lu/ts)',title=f'H={h}, four durations')
    axes[0,0].legend(fontsize=6.5)
    for model,spectrum,label,color,marker in groups:
        fits=sorted([r for r in data['pressure_fits'] if r['recoloring_model']==model and r['mrt_spectrum']==spectrum],key=lambda r:r['height'])
        x=[r['height'] for r in fits]
        axes[1,0].plot(x,[100*(r['records'][-1]['zero_drift_ratio']-1) for r in fits],marker+'-',color=color,label=label)
        axes[1,0].plot(x,[100*(r['closed_pressure_ratio']-1) for r in fits],marker+'--',mfc='none',color=color)
        axes[1,1].plot(x,[100*r['relative_front_mobility_to_pure'] for r in fits],marker+'-',color=color)
    axes[1,0].axhline(0,color='.5',lw=.7);axes[1,0].set(xlabel='Slit width, H',ylabel='Pressure offset from Young (%)')
    axes[1,0].set_title('Filled: zero drift; open: closed pressure',fontsize=8)
    axes[1,1].set(xlabel='Slit width, H',ylabel='Front mobility shift (%)')
    decorate(axes);save(fig,paper,'pressure-fronts')
    fig,axes=plt.subplots(1,2,figsize=(7.7,2.7),constrained_layout=True)
    rows=data['parity']['force_identity_points'];step=[r['step']-rows[0]['step'] for r in rows]
    for axis,d,label in [(axes[0],0,'Tangent'),(axes[1],1,'Normal')]:
        p=np.array([r['physical_projection'][d] for r in rows]);g=np.array([r['force_projection'][d] for r in rows])
        pred=-p[:-1]+.5*np.diff(g)
        axis.plot(step,p,'o-',label='Physical momentum',ms=4)
        axis.plot(step[1:],pred,'x',label='Predicted next momentum',ms=6)
        axis.set(xlabel='Adjacent replay step',ylabel='Mean alternating momentum',title=label)
        axis.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    axes[0].legend(fontsize=7);decorate(axes);save(fig,paper,'parity-force')
    table=[r'\begin{table}[htbp]\centering\small',
        r'\caption{Flat-slit closed/dynamic pressure consistency and measured front mobility. The pressure pair has width $0.0005\,p_Y$ and is tested at both durations; the linear zero estimate uses the three broad pressure values. Mobility includes contact-line/interface effects and is not assigned to a universal numerical error.}',
        r'\begin{tabular}{llrrrr}\toprule',
        r'Recoloring & MRT & $H$ & Closed offset (\%) & Dynamic offset (\%) & Mobility shift (\%) \\ \midrule']
    for row in data['pressure_fits']:
        table.append(f'{row["recoloring_model"].title()} & '+('Published' if row['mrt_spectrum']=='liu2015' else 'Legacy')+
            f' & {row["height"]} & {100*(row["closed_pressure_ratio"]-1):+.3f} & {100*(row["records"][-1]["zero_drift_ratio"]-1):+.3f} & {100*row["relative_front_mobility_to_pure"]:+.2f}'+r' \\')
    table.extend([r'\bottomrule\end{tabular}\end{table}']);(paper/'data/core-reservoir-table.tex').write_text('\n'.join(table)+'\n')
    vals=dict(CoreMatchedBeta=f'{matched:.12f}',
        CoreMatchedProfileError=scientific(max(r['max_matched_profile_error'] for r in data['calibration']['matched_profile_checks'])),
        CoreMaxZeroClosedPercent=f'{100*max(abs(r["zero_minus_closed"]) for r in data["pressure_fits"]):.4f}',
        CoreMaxDurationRatioChange=scientific(max(abs(r['duration_zero_ratio_change']) for r in data['pressure_fits'])),
        CoreMaxPressureBudget=scientific(max(r['budget']['max_relative_budget_residual'] for r in data['fronts']+data['pressure_channels'])),
        CoreMaxHydraulicError=scientific(max(r['relative_profile_error'] for r in data['force_channels']['measurements'])))
    (paper/'data/core-reservoir-numbers.tex').write_text('% Generated from verified core-reservoirs.json.\n'+
        ''.join('\\newcommand{\\'+name+'}{'+value+'}\n' for name,value in vals.items()))
    print(json.dumps(vals))


def decorate(axes):
    for ax,letter in zip(np.asarray(axes).flat,'abcd'):
        ax.spines[['top','right']].set_visible(False);ax.text(-.16,1.02,f'({letter})',transform=ax.transAxes)


def scientific(value):
    coefficient,exponent=f'{value:.2e}'.split('e')
    return coefficient+r'\times10^{'+str(int(exponent))+'}'


def save(fig,paper,name):
    for extension in ('pdf','png'):fig.savefig(paper/'figures'/f'{name}.{extension}')
    plt.close(fig)


if __name__=='__main__':main()
