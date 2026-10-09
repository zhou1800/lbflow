#!/usr/bin/env python3
"""Regenerate capillary plots and manuscript numbers from recorded measurements."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


METHODS=[('legacy','perturbation','Clipped / perturbation','#7b4e2d'),
         ('renormalized','perturbation','Linear / perturbation','#187485'),
         ('renormalized','csf','Linear / CSF','#805b98')]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',default=str(Path(__file__).resolve().parents[1]/
                        'research/drlab-publication/paper/data/capillary-response.json'))
    parser.add_argument('--paper',required=True)
    args=parser.parse_args()
    d=json.loads(Path(args.data).read_text())
    initial={(m['height'],m['angle'],m['wall_gradient'],m['capillary_model']):m for m in d['measurements']}
    held={(m['height'],m['wall_gradient'],m['capillary_model']):m for m in d['long_holds']}
    folder=Path(args.paper)
    (folder/'figures').mkdir(parents=True,exist_ok=True)
    (folder/'data').mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'savefig.dpi':200})
    fig,axes=plt.subplots(1,2,figsize=(7.5,3.6),constrained_layout=True)
    for wall,model,label,color in METHODS:
        selected=[held[(h,wall,model)] for h in (48,32,16)]
        x=[d['planar_width_09']/m['height'] for m in selected]
        y=np.array([100*m['relative_pressure_error'] for m in selected])
        bounds=[]
        for m in selected:
            values=[100*(s['pressure_jump']/m['reference_pressure_jump']-1) for s in
                    m['pressure_sampling_sensitivity'] if 'pressure_jump' in s]
            values+=[100*(m['seed_pressure_coefficient']/m['reference_pressure_coefficient']-1),
                     100*m['relative_pressure_error']]
            bounds.append([min(values),max(values)])
        bounds=np.array(bounds)
        axes[0].errorbar(x,y,yerr=[y-bounds[:,0],bounds[:,1]-y],fmt='o-',color=color,capsize=3,label=label)
        for angle,marker,style in [(60,'o','-'),(90,'^','--')]:
            rows=selected if angle==60 else [initial[(h,angle,wall,model)] for h in (48,32,16)]
            speed=[d['viscosity']*m['history'][-1]['max_speed']/m['sigma'] for m in rows]
            axes[1].semilogy(x,speed,marker+style,color=color,mfc=color if angle==60 else 'white')
    axes[0].axhline(0,color='.5',lw=.8,ls=':')
    axes[0].set(xlabel=r'Bulk interface width / slit width, $W_{0.9}/H$',ylabel='Signed capillary pressure discrepancy (%)')
    axes[0].legend(fontsize=8)
    axes[1].set(xlabel=r'Bulk interface width / slit width, $W_{0.9}/H$',ylabel=r'Residual speed, $\mu U_{\max}/\sigma$')
    axes[1].plot([],[],'o-',color='.3',label='60 degrees, longer hold')
    axes[1].plot([],[],'^--',color='.3',mfc='white',label='90 degrees')
    axes[1].legend(fontsize=8)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    for suffix in ('png','pdf'):fig.savefig(folder/'figures'/f'closed-capillary.{suffix}')
    plt.close(fig)
    ren=held[(48,'renormalized','perturbation')]
    csf=held[(48,'renormalized','csf')]
    ratios=dict(CapillaryRenPressureErrorFine=100*ren['relative_pressure_error'],
                CapillaryCSFPressureErrorFine=100*csf['relative_pressure_error'],
                CapillaryNonneutralReductionFine=ren['history'][-1]['max_speed']/csf['history'][-1]['max_speed'],
                CapillaryNeutralReductionFine=initial[(48,90,'renormalized','perturbation')]['history'][-1]['max_speed']/initial[(48,90,'renormalized','csf')]['history'][-1]['max_speed'],
                CapillaryHoldPressureChange=100*max(abs(m['hold_pressure_coefficient_change']/m['reference_pressure_coefficient']) for m in d['long_holds']))
    formatted={key:format(float(np.ceil(value*1000)/1000) if key=='CapillaryHoldPressureChange'
                          else value,'.3f') for key,value in ratios.items()}
    (folder/'data/capillary-numbers.tex').write_text('% Generated from capillary-response.json; hold-change bound rounded upward.\n'+
        ''.join('\\newcommand{\\'+key+'}{'+value+'}\n' for key,value in formatted.items()))
    rows=[]
    for height in (16,32,48):
        points=[held[(height,w,m)] for w,m,_,_ in METHODS]
        cells=[str(height)]+[f"{100*m['relative_pressure_error']:+.3f}" for m in points]+[
            f"{m['history'][-1]['max_speed']:.3e}" for m in points]
        rows.append(' & '.join(cells)+' \\\\')
    table='''% Generated from the same held measurements as the main figure.
\\begin{table}[htbp]\\centering\\small
\\begin{tabular}{rrrrrrr}\\toprule
H & $E_{clip}$ (\\%) & $E_{lin}$ (\\%) & $E_{CSF}$ (\\%) & $U_{clip}$ & $U_{lin}$ & $U_{CSF}$ \\\\
\\midrule
'''+ '\n'.join(rows)+'''
\\bottomrule\\end{tabular}
\\caption{Direct 60-degree closed-slit results after the longer holds.
Errors are signed relative pressure discrepancies; velocities are in lattice
units. Clipped and linear modes use perturbation; CSF uses the linear gradient.
Input surface tension scales as $H^{-1}$. Sampling and late-field sensitivity
remain in the accompanying JSON rather than being called statistical error.}
\\end{table}
'''
    (folder/'data/capillary-table.tex').write_text(table)
    print(json.dumps(ratios))


if __name__=='__main__':main()
