"""Finite-frequency *sampled* root-separation certificate for exact-axis 1D subspace.
Not a proof of global root simplicity or radiation correctness.
"""
import json
from pathlib import Path
from copy import deepcopy
import mpmath as mp
from task2uv_analytic_oracle import load_material,reduced,analytic_quartic
BASE=Path(__file__).resolve().parents[1]

def run():
 mp.mp.dps=60
 cases=[]
 for name in ('Al6061','PZT-5H','BaTiO3'):
  p=load_material(name)
  for f in ('1e5','1e6','1e7'):
   w=2*mp.pi*mp.mpf(f)
   for case in ('source_tau','synthetic_wtau_0.1','synthetic_wtau_1'):
    pp=deepcopy(p)
    if case!='source_tau':pp['tauq']=mp.mpf(case.split('_')[-1])/w
    r,_=analytic_quartic(reduced(pp),w)
    kvals=[w*mp.sqrt(pp['rho']/pp['C44']), *[abs(z) for z in r]]
    sep=min(abs(a-b) for i,a in enumerate(r) for b in r[i+1:])
    relative=sep/max(kvals)
    # independent 1D SV roots: ±iωsqrt(rho/C44) nondegenerate and distinct from P/thermal.
    shearlam=mp.j*w*mp.sqrt(pp['rho']/pp['C44'])
    mix=min(abs(z-shearlam) for z in r)/max(kvals)
    entry={'material':name,'frequency_hz':f,'tau_case':case,'omega_tau':mp.nstr(w*pp['tauq'],10),
           'minimum_1D_P_thermal_root_gap_normalised':mp.nstr(relative,12),
           'minimum_shear_to_P_thermal_gap_normalised':mp.nstr(mix,12),
           'discriminant_nonzero_sampled':bool(relative>mp.mpf('1e-10')),
           'shear_separate_sampled':bool(mix>mp.mpf('1e-10'))}
    cases.append(entry)
 Path(BASE/'evidence'/'root_separation_certificate.json').write_text(json.dumps({'scope':'sampled exact-axis 1D spectra only, not 8-state k>0 theorem','precision':60,'cases':cases},indent=2)+'\n')
 return cases
if __name__=='__main__':
 z=run()
 print('sampled frequency/tau cases',len(z),'pass',sum(c['discriminant_nonzero_sampled'] and c['shear_separate_sampled'] for c in z))
 print('minimum P-thermal gap',min(float(c['minimum_1D_P_thermal_root_gap_normalised']) for c in z))
 print('minimum SV-to-others gap',min(float(c['minimum_shear_to_P_thermal_gap_normalised']) for c in z))
 if not all(c['discriminant_nonzero_sampled'] and c['shear_separate_sampled'] for c in z):raise SystemExit(1)
