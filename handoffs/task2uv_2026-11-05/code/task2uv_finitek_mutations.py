"""Targeted 8x8 interface-row mutations at k>0.
Uses prior independent finite-k bulk solver; NOT a production assembly test.
Mutant replaces one physical interface row by a constraint on an outgoing coefficient,
solves the changed system, then measures its physical violation of the ORIGINAL law.
"""
import json,argparse
from pathlib import Path
import mpmath as mp
from independent_scattering_limit import scatter,scaled_solve

BASE=Path(__file__).resolve().parents[1]
def run(dps=40):
 mp.mp.dps=dps
 mat=BASE/'source'/'material_constants.csv'
 from independent_finite_k_state_space import load, scale
 from independent_scattering_limit import eig_modes,group,incident
 w=2*mp.pi*mp.mpf('1e6');p1=load(mat,'PZT-5H');p2=load(mat,'BaTiO3')
 kr,_=scale(p1,w);k=kr*mp.mpf('1e-4')
 m1,_=eig_modes(p1,k,w,dps);m2,_=eig_modes(p2,k,w,dps)
 group(m1,k,kr);group(m2,k,kr)
 inc=incident(m1,'P');out1=[q for q in m1 if q['left_out']];out2=[q for q in m2 if q['right_out']]
 assert len(out1)==len(out2)==4
 waves=[('L',q) for q in out1]+[('R',q) for q in out2]
 B=mp.matrix(8);y=mp.matrix(8,1)
 for i in range(8):
  y[i]=-inc['Y'][i]
  for j,(side,m) in enumerate(waves):B[i,j]=(1 if side=='L' else -1)*m['Y'][i]
 x,res=scaled_solve(B,y)
 assert res<mp.mpf('1e-20'),mp.nstr(res,10)
 out={'digits':dps,'pair':'PZT-5H|BaTiO3','incident':'P','k_over_kref':'1e-4',
      'pristine_max_equilibrated_residual':mp.nstr(res,12),'mutants':[]}
 # Note row i support and typical scale are defined by its original entries and incident datum.
 for i,name in [(2,'remove_phi_contact'),(6,'remove_Dz_contact'),(0,'remove_ux_contact'),(4,'remove_shear_traction_contact')]:
  # select a coefficient of nonzero pristine magnitude; ensure row replacement is independent.
  options=sorted(range(8),key=lambda j:abs(x[j]),reverse=True)
  candidate=None
  for jfix in options:
   Bt=mp.matrix(B);yt=mp.matrix(y)
   for j in range(8):Bt[i,j]=1 if j==jfix else 0
   yt[i]=0
   try:
    xt,rest=scaled_solve(Bt,yt)
    candidate=(jfix,xt,rest)
    break
   except Exception: continue
  if candidate is None:
   out['mutants'].append({'name':name,'detected':False,'reason':'all constrained systems singular'});continue
  jfix,xt,rest=candidate
  # Original physical interface row discontinuity, normalized by its uncancelled term magnitudes
  err=abs(sum(B[i,j]*xt[j] for j in range(8))-y[i])
  scale_i=max(sum(abs(B[i,j]*xt[j]) for j in range(8))+abs(y[i]),mp.mpf('1e-100'))
  rel=err/scale_i
  # Also report ROW AGGREGATE error for other seven original rows so the injection is traceable.
  other=max(abs(sum(B[r,j]*xt[j] for j in range(8))-y[r])/
     max(sum(abs(B[r,j]*xt[j]) for j in range(8))+abs(y[r]),mp.mpf('1e-100'))
     for r in range(8) if r!=i)
  delta=max(abs(xt[j]-x[j]) for j in range(8))/max(max(abs(z) for z in x),mp.mpf('1e-100'))
  detected=(rel>mp.mpf('1e-7') and other<mp.mpf('1e-18') and delta>mp.mpf('1e-7'))
  out['mutants'].append({'name':name,'detected':bool(detected),'fixed_column':jfix,
        'original_removed_row_backward_error':mp.nstr(rel,15),'other_rows_max_backward_error':mp.nstr(other,12),
        'amplitude_relative_change':mp.nstr(delta,12),'modified_solve_residual':mp.nstr(rest,8)})
 return out

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--digits',type=int,nargs='+',default=[40]);a.add_argument('--out',default=str(BASE/'evidence'/'finitek_mutation_results.json'));z=a.parse_args()
 res=[run(d) for d in z.digits]
 Path(z.out).write_text(json.dumps({'runs':res},indent=2)+'\n')
 for r in res: print('digits',r['digits'],[(m['name'],m['detected'],m.get('original_removed_row_backward_error')) for m in r['mutants']],flush=True)
 if not all(m['detected'] for r in res for m in r['mutants']): raise SystemExit(1)
