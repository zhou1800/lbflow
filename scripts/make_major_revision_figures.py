#!/usr/bin/env python3
"""Generate scientific figures/tables from independently reanalyzed fields."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from make_core_reservoir_figures import decorate,save,scientific


def select(rows,H,slip=1,sigma=.005,a=None):
    chosen=[]
    for row in rows:
        p=row['parameters'];native='height' in p
        if (p.get('height',p.get('ny')),p.get('physical_slip',p.get('slip',0)),p.get('physical_sigma',p.get('sigma')))!=(H,slip,sigma):continue
        if p.get('physical_force',p.get('force'))!=1e-7 or p.get('physical_amplitude',0):continue
        if a is not None and p.get('transverse_relaxation',0)!=a:continue
        if p.get('transverse_localization','uniform')!='uniform' or 'fixed-b' in row['name']:continue
        chosen.append(row)
    if not chosen:raise LookupError((H,slip,sigma,a))
    return max(chosen,key=lambda row:row['history'][-1]['time'])


def generate(data,paper):
    for name in ('data','figures'):(paper/name).mkdir(parents=True,exist_ok=True)
    native=data['previous_native']+data['native'];refs=data['previous_references']+data['references']
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'savefig.dpi':200})
    colors={16:'#248875',32:'#416592',48:'#a26336'}
    fig,axes=plt.subplots(1,2,figsize=(7.7,3.2),constrained_layout=True)
    groups={(w['height'],w['mode'],w['width']):w for w in data['waves']}
    for mode,marker in ((1,'o'),(2,'s'),(3,'^')):
        rows=sorted([w for (H,m,W),w in groups.items() if m==mode and W<5],key=lambda w:w['height'])
        axes[0].loglog([(16/w['height'])**2 for w in rows],[w['predicted_rate']/(mode*np.pi/16)**4 for w in rows],marker+'-',label=f'W=4.765, m={mode}')
    rows=sorted([w for (H,m,W),w in groups.items() if W==6],key=lambda w:w['height'])
    axes[0].loglog([(16/w['height'])**2 for w in rows],[w['predicted_rate']/(np.pi/16)**4 for w in rows],'d--',label='W=6, m=1')
    x=np.array([1/9,1.]);axes[0].loglog(x,x/24,':',color='.5',label=r'$h^2/24$ formal guide')
    axes[0].set(xlabel=r'$h^2$ (base length units)',ylabel=r'$\Gamma/k^4$ (base units)');axes[0].legend(fontsize=7)
    for mode,width,marker in ((1,4.765235021846117,'o'),(2,4.765235021846117,'s'),(3,4.765235021846117,'^'),(1,6.,'d')):
        w=next(w for w in data['waves'] if w['height']==16 and w['mode']==mode and w['width']==width and w['amplitude']==.01)
        late=w['history'][len(w['history'])//2:];t=np.array([r['time'] for r in late]);t-=t[0]
        A=np.array([r['projected_displacement'] for r in late]);line=axes[1].semilogy(t,A/A[0],marker,ms=4,label=f'W={width:.3g}, m={mode}')[0]
        axes[1].semilogy(t,np.exp(-w['predicted_rate']*t),'-',lw=1,color=line.get_color())
    axes[1].set(xlabel='Time from late-window start (base units)',ylabel=r'$A(t)/A(t_0)$, H=16');axes[1].legend(fontsize=7)
    decorate(axes);save(fig,paper,'major-spectrum')

    fig,axes=plt.subplots(1,2,figsize=(7.7,3.2),constrained_layout=True)
    for ax,slip in zip(axes,(1,4)):
        Hs=[16,32,48];ax.plot(Hs,[1e5*select(native,H,slip)['late_mean_speed'] for H in Hs],'o-',color='#248875',label='Full LBM')
        for a,marker,color in ((0,'s','#a26336'),(1,'^','#416592')):
            ax.plot(Hs,[1e5*select(refs,H,slip,a=a)['late_mean_speed'] for H in Hs],marker+'-',color=color,label=f'FV, a={a}')
            r=[r for r in data['reduced_reference_comparison'] if r['slip']==slip and r['sigma']==.005 and r['strength']==a]
            ax.plot([r['height'] for r in r],[r['predicted_mean_speed']*1e5 for r in r],'x',ms=8,color=color,label=f'Reduced, a={a}')
        if slip==1:ax.plot([16,32],[1e5*select(refs,H,a=.5)['late_mean_speed'] for H in (16,32)],'d--',color='.5',label='FV, a=0.5')
        ax.set(title=rf'Prescribed $\lambda={slip}$',xlabel='Wall-normal lattice nodes, H',ylabel=r'Mean speed $\times10^5$ (base units)');ax.legend(fontsize=7)
    decorate(axes);save(fig,paper,'major-intervention')

    fig,axes=plt.subplots(1,2,figsize=(7.7,3.2),constrained_layout=True)
    sigmas=[.0025,.005,.01]
    for H,ax in zip((16,32),axes):
        ax.plot(np.array(sigmas)*1000,[1e5*select(native,H,sigma=s)['late_mean_speed'] for s in sigmas],'o-',color='#248875',label='Full LBM')
        for a,marker,color in ((0,'s','#a26336'),(1,'^','#416592')):
            ax.plot(np.array(sigmas)*1000,[1e5*select(refs,H,sigma=s,a=a)['late_mean_speed'] for s in sigmas],marker+'-',color=color,label=f'FV, a={a}')
            pred=[next(r for r in data['reduced_reference_comparison'] if (r['height'],r['slip'],r['sigma'],r['strength'])==(H,1,s,a)) for s in sigmas]
            ax.plot(np.array(sigmas)*1000,[1e5*r['predicted_mean_speed'] for r in pred],'x',ms=8,color=color,label=f'Reduced, a={a}')
        ax.set(title=f'H={H}',xlabel=r'Surface tension $\sigma\times10^3$ (base units)',ylabel=r'Mean speed $\times10^5$ (base units)');ax.legend(fontsize=7)
    decorate(axes);save(fig,paper,'major-sigma')

    fig,axes=plt.subplots(1,2,figsize=(7.7,3.2),constrained_layout=True)
    for i,H in enumerate((32,48)):
        vals=[select(refs,H,a=0)['late_mean_speed']]+[next(r for r in refs if r['name']==f'fv-H{H}-slip1-fourth1-{part}')['late_mean_speed'] for part in ('wall','bulk')]+[select(refs,H,a=1)['late_mean_speed']]
        axes[0].plot(range(4),np.array(vals)/select(native,H)['late_mean_speed'],'o-',color=colors[H],label=f'H={H}')
    axes[0].axhline(1,color='.5',ls=':');axes[0].set_xticks(range(4),['Off','Wall','Bulk','Uniform'])
    axes[0].set(xlabel='Specified diagnostic localization',ylabel=r'$\overline{u}/\overline{u}_{LB}$');axes[0].legend(fontsize=8)
    for H in (16,32,48):
        r=next(r for r in data['reduced_lbm_comparison'] if (r['height'],r['slip'],r['sigma'])==(H,1,.005));y=(np.arange(H)+.5)/H
        axes[1].plot(y,np.array(r['measured_volume_shape'])*1000,'o',ms=3,color=colors[H],label=f'LBM, H={H}')
        axes[1].plot(y,np.array(r['predicted_volume_shape'])*1000,'-',lw=1,color=colors[H],label=f'Reduced, H={H}')
    axes[1].set(xlabel=r'Wall-normal coordinate $y/16$',ylabel=r'Phase-volume displacement $\times10^3$');axes[1].legend(fontsize=6)
    decorate(axes);save(fig,paper,'major-localization')

    fig,axes=plt.subplots(1,2,figsize=(7.7,3.2),constrained_layout=True)
    for slip,marker,color in ((1,'o','#248875'),(4,'s','#416592')):
        rows=sorted([r for r in data['manufactured']['runs'] if r['parameters']['slip']==slip and r['parameters']['dt_factor']==.125],key=lambda r:r['parameters']['h'])
        for ax,key in zip(axes,('relative_phase_error','relative_velocity_error')):
            ax.loglog([r['parameters']['h'] for r in rows],[r['history'][-1][key] for r in rows],marker+'-',color=color,label=rf'$\lambda={slip}$')
            if slip==1:
                hs=np.array([.25,1]);ax.loglog(hs,rows[0]['history'][-1][key]*(hs/.25)**2/1.6,':',color='.5',label='Second-order guide')
    half=next(r for r in data['manufactured']['runs'] if r['parameters']['dt_factor']==.0625)
    for ax,key in zip(axes,('relative_phase_error','relative_velocity_error')):
        ax.loglog([half['parameters']['h']],[half['history'][-1][key]],'x',ms=8,color='#a26336',label=r'$\lambda=1$, half time step')
        ax.set(xlabel='Grid spacing h (base units)',ylabel='Relative '+('phase' if key=='relative_phase_error' else 'velocity')+' error');ax.legend(fontsize=8)
    decorate(axes);save(fig,paper,'major-manufactured')

    numbers=dict(MajorWaveCount=str(len(data['waves'])),MajorWavePercent=f'{100*max(abs(r["relative_rate_difference"]) for r in data["waves"]):.5f}',
        MajorWindowDifference=scientific(max(max(abs(w['rate']/r['measured_rate']-1) for w in r['fit_window_checks']) for r in data['waves'])),
        MajorReducedBaselinePercent=f'{np.ceil(1e7*max(abs(r["relative_prediction_error"]) for r in data["reduced_reference_comparison"] if r["sigma"]==.005))/1e5:.5f}',
        MajorReducedAllPercent=f'{np.ceil(1e6*max(abs(r["relative_prediction_error"]) for r in data["reduced_reference_comparison"]))/1e4:.4f}',
        MajorShapePercent=f'{np.ceil(1e5*max(r["relative_volume_shape_error"] for r in data["reduced_lbm_comparison"] if r["sigma"]==.005))/1e3:.3f}',
        MajorFVShapePercent=f'{100*max(r["relative_volume_shape_error"] for r in data["reduced_reference_comparison"]):.3f}',
        MajorReferenceCount=str(data['new_reference_executions']),MajorNewNativeCount=str(data['new_native_executions']))
    for H,label in ((16,'Coarse'),(32,'Middle'),(48,'Fine')):
        numbers['Major'+label+'Percent']=f'{100*(select(refs,H,a=1)["late_mean_speed"]/select(native,H)["late_mean_speed"]-1):.3f}'
    (paper/'data/major-revision-numbers.tex').write_text('% Generated from raw-field reanalysis.\n'+''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in numbers.items()))
    table=[r'\begin{table}[htbp]\centering\small',r'\caption{Prescribed-slip intervention. Mean speeds are multiplied by $10^5$; $\epsilon_1=100(\overline u_{FV,1}/\overline u_{LB}-1)$. The coefficient of the added operator is theoretical and is not fitted.}\label{tab:majorintervention}',r'\begin{tabular}{rrrrrr}\toprule',r'$\lambda$ & H & LBM & FV, a=0 & FV, a=1 & $\epsilon_1$ (\%)\\\midrule']
    for slip in (1,4):
        for H in (16,32,48):
            lb=select(native,H,slip)['late_mean_speed'];v0=select(refs,H,slip,a=0)['late_mean_speed'];v1=select(refs,H,slip,a=1)['late_mean_speed']
            table.append(f'{slip} & {H} & {1e5*lb:.6f} & {1e5*v0:.6f} & {1e5*v1:.6f} & {100*(v1/lb-1):.3f}'+r'\\')
    table.append(r'\bottomrule\end{tabular}\end{table}');(paper/'data/major-intervention-table.tex').write_text('\n'.join(table)+'\n')
    table=[r'\begin{table}[htbp]\centering\small',r'\caption{Full-field hydraulic audit on the same recorded phase: relative $L_2$ differences in percent. The first two columns use the LBM cell-velocity norm, with a specified two-face MAC average. Raw force uses the actual capillary face-force norm; the last column uses the actual-force MAC velocity norm. Both Stokes fields include the same body force.}\label{tab:majorfields}',r'\begin{tabular}{rrrrr}\toprule',r'H & Actual-force velocity & Reference-force velocity & Raw force & Combined-force velocity\\\midrule']
    for H in (16,32,48):
        r=next(r for r in data['full_force_fields']['measurements'] if r['height']==H)
        keys=('relative_actual_force_cell_velocity_error','relative_reference_force_cell_velocity_error','relative_capillary_face_force_difference','relative_force_driven_mac_velocity_difference')
        table.append(str(H)+' & '+' & '.join(f'{100*r[k]:.4f}' for k in keys)+r'\\')
    table.append(r'\bottomrule\end{tabular}\end{table}');(paper/'data/major-fields-table.tex').write_text('\n'.join(table)+'\n')
    table=[r'\begin{table}[htbp]\centering\small',r'\caption{Full-solver spectrum controls. The final column is the maximum signed/amplitude relative rate difference in each grid/mode/width group; no rate is fitted into the prediction.}\label{tab:majorwaves}',r'\begin{tabular}{rrrrr}\toprule',r'H & Mode & $W_{0.9}$ & Predicted rate & Maximum relative difference\\\midrule']
    for H,m,W in sorted(groups):
        rows=[r for r in data['waves'] if (r['height'],r['mode'],r['width'])==(H,m,W)]
        table.append(f'{H} & {m} & {W:.4f} & $'+scientific(rows[0]['predicted_rate'])+'$ & $'+scientific(max(abs(r['relative_rate_difference']) for r in rows))+'$'+r'\\')
    table.append(r'\bottomrule\end{tabular}\end{table}');(paper/'data/major-waves-table.tex').write_text('\n'.join(table)+'\n')
    print(json.dumps(numbers))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--data',required=True);parser.add_argument('--paper',required=True)
    args=parser.parse_args();generate(json.loads(Path(args.data).read_text()),Path(args.paper))


if __name__=='__main__':main()
