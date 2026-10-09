#!/usr/bin/env python3
"""Scientific figures and manuscript numbers from audited phase measurements."""
import argparse
import json
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
    colors={16:'#248875',32:'#416592',48:'#a26336'}
    fig,axes=plt.subplots(1,2,figsize=(7.7,3.1),constrained_layout=True)
    for height,mode in ((16,1),(32,1),(48,1),(16,2)):
        wave=next(r for r in data['waves'] if r['height']==height and r['mode']==mode and r['amplitude']==.01)
        late=wave['history'][len(wave['history'])//2:]
        time=np.array([r['time'] for r in late]);a=np.array([r['projected_displacement'] for r in late]);t=time-time[0]
        axes[0].plot(t,a/a[0],'o',ms=4,color=colors[height],markerfacecolor='none' if mode==2 else colors[height],
            label=f'H={height}, m={mode}')
        axes[0].plot(t,np.exp(-wave['predicted_rate']*t),'--' if mode==2 else '-',lw=1,color=colors[height])
    axes[0].set(xlabel='Time after start of late window (base units)',ylabel=r'$A(t)/A(t_0)$')
    axes[0].legend(fontsize=7)
    for mode,marker in ((1,'o'),(2,'s')):
        rows=sorted([r for r in data['spectra'] if r['mode']==mode],key=lambda r:r['height'])
        axes[1].loglog([(16/r['height'])**2 for r in rows],
            [r['physical_relaxation_rate']/(mode*np.pi/16)**4 for r in rows],marker+'-',label=f'm={mode}')
    x=np.array([1/9,1.]);axes[1].loglog(x,x/24,':',color='.5',label=r'$h^2/24$ guide')
    axes[1].set(xlabel=r'$h^2$, fixed physical width',ylabel=r'$Gamma/k_{phys}^4$ (base units)')
    axes[1].legend(fontsize=8);decorate(axes);save(fig,paper,'phase-spectrum')
    fig,axes=plt.subplots(1,2,figsize=(7.7,3.1),constrained_layout=True)
    rows=data['dynamic_comparison'];x=[r['height'] for r in rows]
    for key,marker,color,label in [('lbm_mean_speed','o','#248875','LBM'),
        ('fv_mean_speed','s','#a26336','Independent Stokes/FV'),
        ('hybrid_mean_speed','^','#416592','Shared-phase ablation')]:
        axes[0].plot(x,[r[key]*1e5 for r in rows],marker+'-',color=color,label=label)
    axes[0].set(xlabel='Slit width, H (lattice nodes)',ylabel=r'Mean speed $\times10^5$ (base units)')
    axes[0].legend(fontsize=8)
    axes[1].semilogy(x,[100*abs(r['relative_lbm_fv_mean_difference']) for r in rows],'s-',color='#a26336',label=r'$|U_{LB}/U_{FV}-1|$')
    axes[1].semilogy(x,[100*abs(r['relative_hybrid_lbm_mean_difference']) for r in rows],'^-',color='#416592',label=r'$|U_{hyb}/U_{LB}-1|$')
    audits=sorted(data['force_audit'],key=lambda r:r['height'])
    axes[1].semilogy(x,[100*abs(np.mean([p['relative_actual_force_speed_difference'] for p in r['points']])) for r in audits],
        'o-',color='#248875',label=r'$|U_{force}/U_{LB}-1|$')
    axes[1].set(xlabel='Slit width, H (lattice nodes)',ylabel='Relative speed difference (%)')
    axes[1].legend(fontsize=8);decorate(axes);save(fig,paper,'kinetic-resistance')
    half=next(r for r in data['references'] if r['name']=='reference-slip1-H32-halfstep-negative')
    full=next(r for r in data['references'] if r['name']=='reference-slip1-H32')
    halfhybrid=next(r for r in data['references'] if r['name']=='hybrid-slip1-H32-halfstep')
    fullhybrid=next(r for r in data['references'] if r['name']=='hybrid-slip1-H32')
    vals=dict(PhaseSpectrumPercent=f'{100*max(abs(r["relative_rate_difference"]) for r in data["waves"]):.5f}',
        PhaseFineFVDifference=f'{100*rows[-1]["relative_lbm_fv_mean_difference"]:.2f}',
        PhaseFineHybridDifference=f'{100*rows[-1]["relative_hybrid_lbm_mean_difference"]:.3f}',
        PhaseFineHydraulicDifference=f'{100*max(abs(p["relative_actual_force_speed_difference"]) for p in audits[-1]["points"]):.4f}',
        PhaseDurationPercent=f'{100*max(abs(r["relative_centroid_speed_change"]) for r in data["duration_controls"]):.5f}',
        PhaseFVTimeSignError=scientific(abs(half['late_mean_speed']/(-full['late_mean_speed'])-1)),
        PhaseHybridTimeError=scientific(abs(halfhybrid['late_mean_speed']/fullhybrid['late_mean_speed']-1)),
        PhaseSymbolError=scientific(max(r['max_absolute_error'] for r in data['symbol_checks'])))
    (paper/'data/phase-kinetics-numbers.tex').write_text('% Generated from audited phase-kinetics.json.\n'+
        ''.join('\\newcommand{\\'+name+'}{'+value+'}\n' for name,value in vals.items()))
    table=[r'\begin{table}[htbp]\centering\small',
        r'\caption{Physical mean speeds for the neutral periodic bridge at $F=10^{-7}$ and prescribed hydraulic slip $\lambda=1$. All speeds are multiplied by $10^5$. The shared-phase calculation is an ablation, not an independent continuum reference.}\label{tab:kineticspeeds}',
        r'\begin{tabular}{rrrrr}\toprule',r'$H$ & LBM & Stokes/FV & Shared phase & Actual-force Stokes \\ \midrule']
    for row,audit in zip(rows,audits):
        table.append(f'{row["height"]} & {1e5*row["lbm_mean_speed"]:.6f} & {1e5*row["fv_mean_speed"]:.6f} & {1e5*row["hybrid_mean_speed"]:.6f} & '+
            f'{1e5*np.mean([p["actual_force_stokes_speed"] for p in audit["points"]]):.6f}'+r' \\')
    table.extend([r'\bottomrule\end{tabular}\end{table}'])
    (paper/'data/phase-kinetics-table.tex').write_text('\n'.join(table)+'\n')
    print(json.dumps(vals))


if __name__=='__main__':main()
