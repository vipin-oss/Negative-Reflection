"""Independent 8-state k>0 and reduced k=0 scattering, built from raw CSV alone.
Does NOT import production Side, polynomial roots, eigenvectors, scattering matrix, or state.py.
Experimental validation; no production promotion, energy/gauge assumptions remain documented.
"""
import argparse,csv,json
from pathlib import Path
import mpmath as mp
from independent_finite_k_state_space import load,operator,balanced,scale,residual,STATE

PAIRS=(('Al6061','Al6061'),('PZT-5H','PZT-5H'),('Al6061','PZT-5H'),('PZT-5H','Al6061'))
ROWS=(0,1,3,4,5,7)
PHI_SCALE=mp.mpf('1e9')

def eig_modes(p,k,w,dps):
    A=operator(p,k,w);B,kr,s=balanced(A,p,w)
    eigvals,V=mp.eig(B)
    out=[]
    for j,a in enumerate(eigvals):
        Y=mp.matrix([s[i]*V[i,j] for i in range(8)])
        la=a*kr
        er=residual(B,a,mp.matrix([V[i,j] for i in range(8)]))
        # Mechanical spectral identity, original eigenvector (not the normalised field).
        electric=abs(a)<=mp.mpf(5)*abs(k)/kr if k else abs(a)<mp.mpf('1e-15')
        # for k=0 exact 2 gauges, only one actual eigenvector; ignore zero eigenpairs
        if electric:
            reference=Y[2]; desired=PHI_SCALE
        elif abs(a)>mp.mpf(100):
            reference=Y[3];desired=mp.mpf(1e3)
        else:
            jdisp=0 if abs(Y[0])>=abs(Y[1]) else 1
            reference=Y[jdisp];desired=mp.mpf(1)
        if reference==0: raise ValueError('normalization reference zero')
        Y=[t*desired/reference for t in Y]
        out.append({'lam':la,'a':a,'Y':Y,'electric':bool(electric),'err':er})
    return out,kr

def group(modes,k,kr):
    """Propagating and decaying mode selection for the chosen time convention.
    Propagating: |Re(lambda)/Im(lambda)|<1e-4, distinguish group via phase sign
    after checking no anomalous flux branch (separate limitation).
    Nonpropagating: decay sign of Re(lambda).
    """
    for q in modes:
        la=q['lam']; rel=abs(mp.re(la))/max(abs(mp.im(la)),mp.mpf('1e-100'))
        isprop=(rel<mp.mpf('1e-4') and abs(mp.im(la))>kr*mp.mpf('1e-3'))
        q['prop']=isprop
        q['left_out']=(mp.im(la)<0) if isprop else (mp.re(la)>0)
        q['right_out']=(mp.im(la)>0) if isprop else (mp.re(la)<0)
        q['incident']=(mp.im(la)>0 and isprop)
    return modes

def incident(modes,kind):
    cands=[x for x in modes if x['incident']]
    # incoming +z modes, longitudinal acoustic has smaller |lambda| than SV
    cands=sorted(cands,key=lambda x:abs(x['lam']))
    if len(cands)!=2: raise RuntimeError('unexpected incident branches '+str([(mp.nstr(x['lam'],5),x['prop']) for x in modes]))
    return cands[0] if kind=='P' else cands[1]

def scaled_solve(B,y):
    n,m=B.rows,B.cols
    if n!=m: raise ValueError(f'non-square: {n}x{m}')
    rr=[max([abs(B[i,j]) for j in range(m)]+[abs(y[i]),mp.mpf('1e-150')]) for i in range(n)]
    A=mp.matrix(n,m);f=mp.matrix(n,1)
    for i in range(n):
        f[i]=y[i]/rr[i]
        for j in range(m):A[i,j]=B[i,j]/rr[i]
    cc=[max([abs(A[i,j]) for i in range(n)]+[mp.mpf('1e-150')]) for j in range(m)]
    for i in range(n):
        for j in range(m):A[i,j]/=cc[j]
    z=mp.lu_solve(A,f)
    x=[z[j]/cc[j] for j in range(m)]
    norm=mp.mpf(0)
    for i in range(n):
        norm=max(norm,abs(sum(B[i,j]*x[j] for j in range(m))-y[i])/rr[i])
    return x,norm

def scatter(matpath,left,right,kind,k_over,dps,exact=False):
    mp.mp.dps=dps
    w=2*mp.pi*mp.mpf('1e6')
    p1,p2=load(matpath,left),load(matpath,right)
    kr,_=scale(p1,w)
    k=kr*k_over if not exact else mp.mpf(0)
    M1,_=eig_modes(p1,k,w,dps);M2,_=eig_modes(p2,k,w,dps)
    group(M1,k,kr);group(M2,k,kr)
    inc=incident(M1,kind)
    side1=[m for m in M1 if m['left_out'] and (not exact or abs(m['a'])>mp.mpf('1e-15'))]
    side2=[m for m in M2 if m['right_out'] and (not exact or abs(m['a'])>mp.mpf('1e-15'))]
    rows=ROWS if exact else tuple(range(8))
    if len(side1)+len(side2)!=len(rows):
        raise RuntimeError(f'count {left} {right} {kind} {k_over} exact={exact}: {len(side1)}+{len(side2)} rows={len(rows)}; left={[mp.nstr(m["a"],5) for m in M1]}')
    modes=[('L',m) for m in side1]+[('R',m) for m in side2]
    B=mp.matrix(len(rows),len(modes));y=mp.matrix(len(rows),1)
    for i,r in enumerate(rows):
        y[i]=-inc['Y'][r]
        for j,(side,m) in enumerate(modes):B[i,j]=(1 if side=='L' else -1)*m['Y'][r]
    x,res=scaled_solve(B,y)
    fL=inc['Y'][:]
    fR=[mp.mpc(0)]*8
    for j,(side,m) in enumerate(modes):
        for i in range(8):
            if side=='L':fL[i]+=x[j]*m['Y'][i]
            else:fR[i]+=x[j]*m['Y'][i]
    C=0
    if exact:
        # under no-static-field sector, phi constants do not alter stresses or E_z.
        C=fL[2]-fR[2]
        fR[2]+=C
    jump=[fL[i]-fR[i] for i in range(8)]
    # electric derivatives: Ez=-phi', where φ' is from full physical operator
    AL=operator(p1,k,w);AR=operator(p2,k,w)
    EzL=-(AL*mp.matrix(fL))[2]
    EzR=-(AR*mp.matrix(fR))[2]
    efields={'ExL':-mp.j*k*fL[2],'ExR':-mp.j*k*fR[2],'EzL':EzL,'EzR':EzR,'DzL':fL[6],'DzR':fR[6]}
    def scale_phys(i):
        return max(abs(fL[i]),abs(fR[i]),mp.mpf('1e-80'))
    maxcont=max(abs(jump[i])/scale_phys(i) for i in rows)
    return {'left':left,'right':right,'kind':kind,'digits':dps,'k_over_kref':mp.nstr(k_over,12),
      'k':mp.nstr(k,14),'exact':exact,'nL':len(side1),'nR':len(side2),
      'solve_residual':mp.nstr(res,9),'max_contact_residual':mp.nstr(maxcont,9),
      'electric_amplitudes':[mp.nstr(x[j],12) for j,(side,m) in enumerate(modes) if m['electric']],
      'electrostatic_mode_phi':[mp.nstr(m['Y'][2],12) for side,m in modes if m['electric']],
      'incident_norm':mp.nstr(max(abs(inc['Y'][i]) for i in (0,1)),7),
      'phi_left':mp.nstr(fL[2],15),'phi_right':mp.nstr(fR[2],15),
      'phi_const':mp.nstr(C,15),
      'fields':{STATE[i]:mp.nstr(fL[i],20) for i in range(8)},
      'efields':{key:mp.nstr(val,20) for key,val in efields.items()},
      'raw':{'fL':fL,'fR':fR,'x':x,'modes':modes,'kr':kr,'E':efields},
      'status':'PASS_BC' if res<mp.mpf('1e-20') and maxcont<mp.mpf('1e-17') else 'WARN_BC'}

def run(matpath,out,dps_list=(40,70),ratios=('1e-2','1e-4','1e-6')):
    results=[];exceptions=[]
    for dps in dps_list:
      for left,right in PAIRS:
       for kind in ('P','SV'):
        try:
          zero=scatter(matpath,left,right,kind,mp.mpf('0'),dps,exact=True)
          for ratio in ratios:
            f=scatter(matpath,left,right,kind,mp.mpf(ratio),dps)
            ref=zero['raw']
            # Physical field convergence, fixed scale at k=0
            compare=[0,1,3,4,5,7]
            errors=[]
            for i in compare:
                a=f['raw']['fL'][i];b=ref['fL'][i]
                scalev=max(abs(b),abs(ref['fR'][i]), {0:mp.mpf(1),1:mp.mpf(1),3:mp.mpf('1e6'),4:mp.mpf('1e14'),5:mp.mpf('1e14'),7:mp.mpf('1e9')}[i])
                errors.append(abs(a-b)/scalev)
            # k-dependent quantity normalized to a nonzero *fixed* electric scale
            kr=zero['raw']['kr']; phi_scale=mp.mpf('1e9')
            Ez_delta=abs(f['raw']['E']['EzL']-ref['E']['EzL'])/max(abs(ref['E']['EzL']),kr*phi_scale,mp.mpf('1e-20'))
            dz_scale=kr*mp.mpf('1e9')*max(load(matpath,left)['eps33'],load(matpath,right)['eps33'])
            dz_norm=max(abs(f['raw']['E']['DzL']),abs(f['raw']['E']['DzR']))/dz_scale
            phi_jump=abs(f['raw']['fL'][2]-f['raw']['fR'][2])/phi_scale
            row={key:val for key,val in f.items() if key!='raw' and key not in ('fields','efields')}
            row.update({'max_mech_thermal_diff_fixed_scale':mp.nstr(max(errors),14),
                        'Ez_diff_fixed_scale':mp.nstr(Ez_delta,14),
                        'Dz_over_fixed_scale':mp.nstr(dz_norm,14),
                        'phi_jump_over_fixed_scale':mp.nstr(phi_jump,14),
                        'exact_6x6_residual':zero['solve_residual'],
                        'exact_relative_phi_C':zero['phi_const']})
            results.append(row)
            print(left,right,kind,dps,ratio,'M',mp.nstr(max(errors),4),'Dz',mp.nstr(dz_norm,4),'A',f['electric_amplitudes'],flush=True)
        except Exception as e:
          exceptions.append({'pair':f'{left}|{right}','kind':kind,'digits':dps,'error':repr(e)})
          print('ERROR',left,right,kind,dps,str(e),flush=True)
    full={'scope':'independent bulk 8-state eigenbasis and contact solve, not rigorously validated causal radiation and not production',
          'results':results,'exceptions':exceptions,'limitations':['no material conductors modeled','eps0 Al is a dielectric surrogate','radiation selection uses eigenvalue/phase rule rather than independent coupled flux','some electric mode amplitudes depend on phi=1e9 normalization']}
    Path(out).write_text(json.dumps(full,indent=2)+'\n')
    return full

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--materials',required=True);a.add_argument('--out',required=True);a.add_argument('--digits',type=int,nargs='+',default=[40]);a.add_argument('--ratios',nargs='+',default=['1e-2','1e-4']);o=a.parse_args()
    run(o.materials,o.out,o.digits,o.ratios)
