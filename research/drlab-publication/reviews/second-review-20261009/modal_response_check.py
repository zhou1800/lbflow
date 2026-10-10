"""Read-only audit: resolve which transverse modes carry the reported response.
No evolved-front input, coefficient fit, or change to publication source.
"""
import argparse,hashlib,json,sys
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--sources',type=Path,default=Path(__file__).resolve().parents[4]/'scripts')
args=parser.parse_args()
if args.output.exists():raise FileExistsError(args.output)
sys.path.insert(0,str(args.sources.resolve()))
import numpy as np
from stokes_phase_reference import StokesSlit,capillary_faces,phase_rhs,logistic_bridge,match_logistic_width,transverse_relaxation_rhs
from front_mobility import WIDTH_REFERENCE
from phase_wave_spectrum import spectrum,matrix
rows=[]
for H,slip in ((16,1),(32,1),(48,1),(48,4)):
 h=16/H; nx=4*H; b=match_logistic_width(WIDTH_REFERENCE,h); solver=StokesSlit(H,slip=slip,fast=True)
 x=(np.arange(nx)+.5)*h; p=logistic_bridge(nx,h,b)
 derivative=-np.sign(x-32)*2*b*p*(1-p); weight=h*derivative[:nx//2].sum();phi=np.broadcast_to(p,(H,nx)).copy()
 def project(rate):return -h*rate[:,:nx//2].sum(axis=1)/weight
 def response(eta):
  delta=-eta[:,None]*derivative[None]; ep=1e-4
  plus=capillary_faces(phi+ep*delta,h,.005,b);minus=capillary_faces(phi-ep*delta,h,.005,b)
  u,v=solver.solve((plus[0]-minus[0])/(2*ep),(plus[1]-minus[1])/(2*ep))
  return project(phase_rhs(phi,u,v,h,b,diffusion=0.)),float(u.mean())
 col=[response(e) for e in np.eye(H)]; C=np.array([v[0] for v in col]).T; cmean=np.array([v[1] for v in col]);C-=C.sum(axis=1)[:,None]/H
 u0,v0=solver.solve(np.full((H,nx),1e-7),np.zeros((H+1,nx)));g=project(phase_rhs(phi,u0,v0,h,b,diffusion=0.));mean0=float(u0.mean())
 Q=np.cos(np.pi*(np.arange(H)[:,None]+.5)*np.arange(H)[None]/H); Q[:,0]/=np.sqrt(H);Q[:,1:]*=np.sqrt(2/H)
 rates=h*h/24*(4/h**2*np.sin(np.pi*np.arange(H)/(2*H))**2)**2
 def solve(op,Qsub):
  # Work in an orthonormal cosine subspace, remove the constant gauge.
  o=Qsub.T@op@Qsub;gd=Qsub.T@g;N=len(gd);mat=np.zeros((N+1,N+1));mat[:N,:N]=o;mat[:N,-1]=-(Qsub.T@np.ones(H));mat[-1,0]=1
  z=np.linalg.solve(mat,np.r_[-gd,0]);eta=Qsub@z[:-1]
  return {'phase_speed':float(z[-1]),'mean_speed':float(mean0+cmean@eta),'max_eta':float(abs(eta).max()),'eta':eta.tolist(),'cosine_amplitudes':(Q.T@eta).tolist()}
 base=solve(C-(Q*rates)@Q.T,Q)
 item={'height':H,'slip':slip,'full_fourth':base,'galerkin':[],'low_mode_intervention':[]}
 for cut in sorted(set((2,3,4,6,8,12,H//2,H-1))):
  if cut>=H:continue
  z=solve(C-(Q*rates)@Q.T,Q[:,:cut+1]);z.update(max_mode=cut,relative_mean_change=z['mean_speed']/base['mean_speed']-1); item['galerkin'].append(z)
  mask=np.arange(H)<=cut;z=solve(C-(Q*(rates*mask))@Q.T,Q);z.update(max_regularized_mode=cut,relative_mean_change=z['mean_speed']/base['mean_speed']-1);item['low_mode_intervention'].append(z)

 # Compare the formal all-mode regularization against the actual rest-map
 # eigenvalues and the instantaneous integral response; neither is an LBM run.
 sp=spectrum(H,1);pr=np.array(sp['profile']);dp=(np.roll(pr,-1)-np.roll(pr,1))/2
 mode_rows=[];spectral_rates=np.zeros(H);instantaneous_rates=np.zeros(H)
 for m in range(1,H):
  T=matrix(pr,sp['beta'],m*np.pi/H);values,vectors=np.linalg.eig(T)
  overlaps=abs(vectors.conj().T@dp)/np.linalg.norm(vectors,axis=0)/np.linalg.norm(dp)
  j=int(np.argmax(overlaps));val=values[j]
  if abs(val.imag)>1e-10 or not 0<val.real<1:raise ValueError((H,m,val))
  spectral_rates[m]=(1-val.real)/h**2
  instantaneous_rates[m]=float(-((T@dp-dp)[:nx//2].sum()/dp[:nx//2].sum())/h**2)
  mode_rows.append(dict(mode=m,eigenvalue=float(val.real),translation_overlap=float(overlaps[j]),exact_semidiscrete_rate=float(spectral_rates[m]),exact_log_rate=float(-np.log(val.real)/h**2),instantaneous_projection_rate=float(instantaneous_rates[m]),formal_discrete_rate=float(rates[m]),relative_formal_vs_spectral=float(rates[m]/spectral_rates[m]-1)))
 item['all_mode_spectrum']=mode_rows
 for name,rr in [('exact_eigenvalue_semidiscrete',spectral_rates),('instantaneous_projection',instantaneous_rates)]:
  z=solve(C-(Q*rr)@Q.T,Q);z.update(relative_mean_change=z['mean_speed']/base['mean_speed']-1)
  item[name]=z
 item['untreated_galerkin']=[]
 for cut in sorted(set((2,4,6,8,12,H//2,H-1))):
  if cut>=H:continue
  z=solve(C,Q[:,:cut+1]);z['max_mode']=cut;item['untreated_galerkin'].append(z)
 rows.append(item)
root=args.output.parent
root.mkdir(parents=True,exist_ok=True)
out={'scope':'Additional reduced-matrix sensitivity, not new nonlinear reference or LBM trajectories','cases':rows,'source_sha256':{Path(m.__file__).name:hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest() for m in (sys.modules['stokes_phase_reference'],sys.modules['front_mobility'])}}
args.output.write_text(json.dumps(out,indent=2)+'\n')
for z in rows:
 print('H',z['height'],'slip',z['slip'],'full',z['full_fourth']['mean_speed'])
 print('Galerkin',[(v['max_mode'],round(v['relative_mean_change']*100,5)) for v in z['galerkin']])
 print('Spectral replacements',[(k,z[k]['relative_mean_change']*100) for k in ('exact_eigenvalue_semidiscrete','instantaneous_projection')])
 print('Relevant spectra',[(v['mode'],v['translation_overlap'],v['relative_formal_vs_spectral']*100) for v in z['all_mode_spectrum'] if v['mode'] in (2,4,6,8,12)])
 print('Untreated Galerkin',[(v['max_mode'],v['mean_speed']) for v in z['untreated_galerkin']])
 print('Low-mode intervention',[(v['max_regularized_mode'],round(v['relative_mean_change']*100,5)) for v in z['low_mode_intervention']])
