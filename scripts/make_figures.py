#!/usr/bin/env python3
"""Publication figures from recorded measurements; no fitted or invented data."""
import argparse,json,math,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser()
p.add_argument('--root',help='Optional recorded raw-study root for rebuilding measurement summaries')
p.add_argument('--data',default=str(Path(__file__).resolve().parents[1]/'research/drlab-publication/paper/data/measurements.json'))
p.add_argument('--paper',required=True)
a=p.parse_args();root=Path(a.root) if a.root else None;paper=Path(a.paper)
recorded=json.loads(Path(a.data).read_text()) if root is None else None
out=paper/'figures';out.mkdir(parents=True,exist_ok=True)
dataout=paper/'data';dataout.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,
                    'figure.dpi':150,'savefig.dpi':200,'pdf.fonttype':42})
colors={'legacy':'#9b5040','renormalized':'#187485','arithmetic':'#9b5040','harmonic':'#187485'}
def read(path):return json.loads(path.read_text())
def save(fig,name):
    fig.savefig(out/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(out/(name+'.png'),bbox_inches='tight')
    plt.close(fig)
def measurements(folder,file):return read(root/folder/file)['measurements']

verification=recorded['verification'] if recorded is not None else read(root/('verification-release' if (root/'verification-release/verification.json').exists() else 'verification-final')/'verification.json')
channel=[m for m in verification['measurements'] if m['kind']=='channel' and 'scsp' in m['name']]
fig,axes=plt.subplots(1,2,figsize=(7.4,3.5),constrained_layout=True)
for tau in (.6,1.4):
    for eq,style in [('SRT','o-'),('MRT','s--')]:
        ms=sorted([m for m in channel if m['tau']==tau and eq in m['name']],key=lambda m:m['height'])
        axes[0].loglog([m['height'] for m in ms],[m['relative_l2_error'] for m in ms],style,
                       label=f'{eq}, tau={tau}')
ref=np.array([8,16,32]);axes[0].loglog(ref,.02*(8/ref)**2,':',color='.5',label='second-order slope')
axes[0].set(xlabel='Fluid wall spacing H (nodes)',ylabel='Relative L2 profile error')
axes[0].legend(fontsize=9,ncol=2)
chosen=next(m for m in channel if m['tau']==.6 and 'MRT' in m['name'] and m['height']==32)
y=(np.arange(32)+.5)/32
axes[1].plot(y,np.array(chosen['analytic_profile'])/max(chosen['analytic_profile']),color='.2',label='analytic')
axes[1].plot(y,np.array(chosen['measured_profile'])/max(chosen['analytic_profile']),'o',ms=3,mfc='none',label='MRT')
axes[1].set(xlabel='y / H',ylabel='u / analytic maximum');axes[1].legend()
save(fig,'channel-verification')

layered=recorded['layered'] if recorded is not None else measurements('layered-v3','layered-flow.json')+measurements('layered-harmonic-v4','layered-flow.json')
fig,axes=plt.subplots(1,2,figsize=(7.4,3.5),constrained_layout=True)
for rule in ('arithmetic','harmonic'):
    ms=sorted([m for m in layered if m['tau']==[.6,1.4] and m['interpolation']==rule],key=lambda m:m['height'])
    axes[0].loglog([m['height'] for m in ms],[m['relative_l2_error'] for m in ms],'o-',color=colors[rule],label=rule+' / sharp')
    axes[0].loglog([m['height'] for m in ms],[m['relative_diffuse_l2_error'] for m in ms],'s--',color=colors[rule],label=rule+' / diffuse')
    m=ms[-1];y=(np.arange(m['height'])+.5)/m['height'];scale=max(m['analytic_profile'])
    axes[1].plot(y,np.array(m['profile'])/scale,'-',color=colors[rule],label=rule)
    axes[1].plot(y,np.array(m['diffuse_continuum_profile'])/scale,'--',color=colors[rule],lw=.8)
axes[1].plot(y,np.array(m['analytic_profile'])/scale,':',color='.2',label='sharp reference')
axes[0].set(xlabel='H (nodes)',ylabel='Relative L2 profile error');axes[0].legend(fontsize=9)
axes[1].set(xlabel='y / H',ylabel='u / sharp maximum');axes[1].legend(fontsize=9)
save(fig,'layered-flow')

if recorded is not None:
    cap=recorded['capillary']
else:
    cap=measurements('capillary-relaxed-v2','capillary.json')
    cap+=measurements('capillary-v2','capillary.json')
    extra=measurements('capillary-extra-v2','capillary.json')
    cap=[m for m in cap if not (m['name']=='laplace-d3-R16-1rank')]+extra
fig,axes=plt.subplots(1,2,figsize=(7.4,3.5),constrained_layout=True)
for ndim,color in [(2,'#187485'),(3,'#a46529')]:
    ms=[m for m in cap if m.get('measured_radius') and m['dimensions']==ndim and
         ((ndim==2 and 'relaxed' in m['name'] and 'tau' not in m['name']) or ndim==3)]
    ms=sorted(ms,key=lambda m:m['measured_radius'])
    x=np.array([(ndim-1)/m['measured_radius'] for m in ms]);dp=np.array([m['pressure_jump'] for m in ms])
    axes[0].plot(x,dp,'o',color=color,label=f'D{ndim}Q{9 if ndim==2 else 19}')
    axes[0].plot(x,.005*x,'-',color=color,lw=.8)
    axes[1].loglog([m['measured_radius'] for m in ms],[m['relative_laplace_error'] for m in ms],'o-',color=color,label=f'D{ndim}')
axes[0].set(xlabel='Curvature (d - 1) / fitted R',ylabel='Bulk pressure jump');axes[0].legend()
axes[1].set(xlabel='Fitted radius (nodes)',ylabel='Relative Laplace error');axes[1].legend()
save(fig,'laplace-verification')

if recorded is not None:
    wet=recorded['wetting']
else:
    wet={}
    for m in measurements('wetting-v3','wetting.json'):
        wet[(m['method'],m['target_angle'],m['initial_radius'],m['beta'])]=m
    for file in sorted((root/'wetting-steady-v3/results').glob('*.measurement.json')):
        m=read(file);wet[(m['method'],m['target_angle'],m['initial_radius'],m['beta'])]=m
    wet=list(wet.values())
fig,axes=plt.subplots(2,3,figsize=(7.4,6.2),constrained_layout=True)
for col,theta in enumerate((60,90,120)):
    for rule in ('legacy','renormalized'):
        ms=sorted([m for m in wet if m['target_angle']==theta and m['method']==rule and m['beta']==.7],key=lambda m:m['initial_radius'])
        values=np.array([m['principal']['angle_degrees']-theta for m in ms])
        lower=np.array([m['principal']['angle_degrees']-m['sensitivity_min'] for m in ms])
        upper=np.array([m['sensitivity_max']-m['principal']['angle_degrees'] for m in ms])
        axes[0,col].errorbar([m['initial_radius'] for m in ms],values,
                            yerr=np.vstack((np.maximum(0,lower),np.maximum(0,upper))),
                            fmt='o-',capsize=3,color=colors[rule],label=rule)
        axes[1,col].plot([m['initial_radius'] for m in ms],[m['max_spurious_speed'] for m in ms],'o-',color=colors[rule],label=rule)
    axes[0,col].axhline(0,color='.5',ls=':',lw=.8);axes[0,col].set_title(f'Target {theta} degrees')
    axes[0,col].set(ylabel='Measured - target (degrees)',xlabel='Initial radius (nodes)')
    axes[1,col].set(xlabel='Initial radius (nodes)',ylabel='Maximum spurious speed')
    axes[1,col].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
axes[0,0].legend(fontsize=9)
save(fig,'wetting-sensitivity')

assembled=dict(verification=verification,layered=layered,capillary=cap,wetting=wet,
               provenance=recorded['provenance'] if recorded is not None else 'See case receipts and raw VTI fields under the research result root. Interrupted R48 outputs are excluded.')
(dataout/'measurements.json').write_text(json.dumps(assembled,indent=2))
with (dataout/'wetting.csv').open('w',newline='') as stream:
    writer=csv.writer(stream,lineterminator="\n");writer.writerow(['method','target','initial_radius','beta','W09','angle','sensitivity_min','sensitivity_max','late_change','max_spurious_speed'])
    for m in sorted(wet,key=lambda m:(m['beta'],m['initial_radius'],m['target_angle'],m['method'])):
        writer.writerow([m['method'],m['target_angle'],m['initial_radius'],m['beta'],m['width_09'],m['principal']['angle_degrees'],m['sensitivity_min'],m['sensitivity_max'],m['late_angle_change'],m['max_spurious_speed']])
print('Generated four PDF/PNG figures and assembled',len(wet),'wetting measurements.')
