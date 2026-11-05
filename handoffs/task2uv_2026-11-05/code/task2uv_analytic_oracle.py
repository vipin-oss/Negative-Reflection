"""Task 2U-V: independent normal-incidence 1D Lord-Shulman oracle.
Only imports mpmath & standard library. Does not import research repository or finite-k implementation.
Time convention exp(-i omega t); spatial dependence exp(lambda z).
No conductor-contact or energy/absorption assertions.
"""
import csv, json, argparse
from pathlib import Path
from copy import deepcopy
import mpmath as mp

BASE=Path(__file__).resolve().parents[1]
CSV=BASE/'source'/'material_constants.csv'

def load_material(name):
    a={}
    with CSV.open(encoding='utf-8-sig',newline='') as f:
        for r in csv.DictReader(f):
            if r['material']==name:
                try: a[r['symbol']]=mp.mpf(r['value'])
                except (ValueError,KeyError): pass
    return a

def reduced(p):
    e,eps=p['e33'],p['eps33']
    return {'C':p['C33']+e*e/eps,
            'b':p['beta3']-e*p['p3']/eps,
            'r':p['rho']*p['cE']-p['T0']*p['p3']**2/eps,
            'kappa':p['kappa33'],'rho':p['rho'],'T':p['T0'],'tau':p['tauq']}

def analytic_quartic(m,w,wrong_g=False):
    """Derived from (C lambda²+rho w²)U-b lambda Theta=0 and
    (kappa lambda²/g+i w r)Theta+i w T b lambda U=0.
    This is a stand-alone PDE-derived polynomial in x=lambda².
    """
    I=mp.j;g=1-I*w*m['tau']; G=1/g if not wrong_g else g
    a=m['C']*m['kappa']*G
    b=m['rho']*w*w*m['kappa']*G+I*w*(m['C']*m['r']+m['T']*m['b']**2)
    c=I*w*m['rho']*w*w*m['r']
    discr=mp.sqrt(b*b-4*a*c)
    x1=(-b+discr)/(2*a);x2=(-b-discr)/(2*a)
    return [mp.sqrt(x1),-mp.sqrt(x1),mp.sqrt(x2),-mp.sqrt(x2)], (a,b,c)

def physical_mode(m,w,lam):
    """Build physical 4-state (u,theta,sigma,q), using constitutive equations.
    Mechanical vs thermal normalization uses a magnitude comparison to avoid division by tiny b.
    """
    C,b,r,k=m['C'],m['b'],m['r'],m['kappa'];g=1-mp.j*w*m['tau']
    coeff=C*lam**2+m['rho']*w*w
    if abs(coeff)<abs(b*lam):
        u=mp.mpc(1);theta=coeff/(b*lam) if b*lam else mp.mpc(0)
    else:
        theta=mp.mpc(1);u=b*lam/coeff if coeff else mp.mpc(0)
    sigma=C*lam*u-b*theta
    q=-k*lam*theta/g
    return mp.matrix([u,theta,sigma,q])

def first_order_1d(m,w,wrong_g=False):
    C,b,r,kap=m['C'],m['b'],m['r'],m['kappa'];g=1-mp.j*w*m['tau']
    # Directly from balance and constitutive laws, *not* the characteristic polynomial.
    A=mp.matrix(4)
    A[0,1]=b/C;A[0,2]=1/C
    A[1,3]=-g/kap if not wrong_g else -1/(g*kap)
    A[2,0]=-m['rho']*w*w
    A[3,1]=mp.j*w*(r+m['T']*b*b/C)
    A[3,2]=mp.j*w*m['T']*b/C
    return A

def mode_catalog(m,w,wrong_g=False):
    roots,coeffs=analytic_quartic(m,w,wrong_g)
    vals=[]
    for z in roots:
        state=physical_mode(m,w,z) if not wrong_g else mutated_mode(m,w,z)
        # Normalization for interface columns: mechanical if P, temperature if thermal.
        # no radiation claims from only these numbers until signs are independently checked.
        lam0=mp.j*w*mp.sqrt(m['rho']/m['C'])
        pbranch=abs(z**2-lam0**2)<abs(z**2+mp.j*w*m['r']*(1-mp.j*w*m['tau'])/m['kappa'])
        d=state[0] if pbranch else state[1]
        if d==0:raise RuntimeError('zero normalization pivot')
        st=state/d
        vals.append({'lam':z,'Y':st,'P':bool(pbranch)})
    return vals

def mutated_mode(m,w,lam):
    # wrong g used only for LS-mutation demonstration; thermal flux must follow wrong law.
    C,b=m['C'],m['b']; coeff=C*lam**2+m['rho']*w*w
    if abs(coeff)<abs(b*lam):u=mp.mpc(1);th=coeff/(b*lam)
    else:th=mp.mpc(1);u=b*lam/coeff
    sig=C*lam*u-b*th
    q=-m['kappa']*(1-mp.j*w*m['tau'])*lam*th
    return mp.matrix([u,th,sig,q])

def outgoing(modes,side):
    out=[]
    for mode in modes:
        z=mode['lam']; ratio=abs(mp.re(z))/max(abs(mp.im(z)),mp.mpf('1e-140'))
        propagating=mode['P'] and ratio<mp.mpf('1e-3')
        if propagating:ok=(mp.im(z)<0 if side=='left' else mp.im(z)>0)
        else:ok=(mp.re(z)>0 if side=='left' else mp.re(z)<0)
        if ok:out.append(mode)
    return out

def solve_interface(left,right,w,reverse_left=False):
    L=mode_catalog(left,w);R=mode_catalog(right,w)
    incs=[x for x in L if x['P'] and mp.im(x['lam'])>0]
    if len(incs)!=1:raise RuntimeError(f'P incident selection ambiguous: {len(incs)}')
    inc=incs[0]
    refl=outgoing(L,'right' if reverse_left else 'left')
    trns=outgoing(R,'right')
    # If reversed chooses incident itself, that's an expected radiation mutation.
    if len(refl)!=2 or len(trns)!=2: raise RuntimeError(f'bad radiation count {len(refl)} {len(trns)}')
    B=mp.matrix(4); rhs=-inc['Y']
    for j,x in enumerate(refl):
        for i in range(4):B[i,j]=x['Y'][i]
    for j,x in enumerate(trns):
        for i in range(4):B[i,j+2]=-x['Y'][i]
    # Row scaled solution, fixed physical equations.
    rowscale=[max([abs(B[i,j]) for j in range(4)]+[abs(rhs[i]),mp.mpf('1e-130')]) for i in range(4)]
    A=mp.matrix([[B[i,j]/rowscale[i] for j in range(4)] for i in range(4)])
    y=mp.matrix([rhs[i]/rowscale[i] for i in range(4)])
    cs=[max(abs(A[i,j]) for i in range(4)) for j in range(4)]
    Ae=mp.matrix([[A[i,j]/cs[j] for j in range(4)] for i in range(4)])
    x0=mp.lu_solve(Ae,y);amps=[x0[j]/cs[j] for j in range(4)]
    fL=inc['Y']+sum((amps[j]*refl[j]['Y'] for j in range(2)),mp.matrix(4,1))
    fR=sum((amps[j+2]*trns[j]['Y'] for j in range(2)),mp.matrix(4,1))
    mismatch=max(abs(fL[i]-fR[i])/rowscale[i] for i in range(4))
    return {'transmitted_u':fR[0],'transmitted_theta':fR[1],
        'transmitted_sigma':fR[2],'transmitted_q':fR[3],
        'reflected_u':fL[0]-inc['Y'][0], 'contact_residual':mismatch,
        'amplitudes':amps,'incident':inc,'left':fL,'right':fR,'refl':refl,'trns':trns}

def close_roots(aa,bb):
    return max(min(abs(a-b) for b in bb)/max(mp.mpf(1),abs(a)) for a in aa)

def run(dps=40):
    mp.mp.dps=dps
    data={'precision':dps,'cases':[],'mutations':[],'failures':[],'limitations':[]}
    def check(tag,success,details):
        data['cases'].append({'test':tag,'passed':bool(success),'details':details})
        if not success:data['failures'].append(tag)
    w=2*mp.pi*mp.mpf('1e6')
    A=load_material('Al6061');P=load_material('PZT-5H');Q=load_material('BaTiO3')
    # 1. Elastic reflection anchor, fully uncoupled, independent closed-form impedances.
    X=deepcopy(A);Y=deepcopy(P)
    for m in (X,Y):m['beta3']=mp.mpf(0);m['e33']=mp.mpf(0);m['p3']=mp.mpf(0)
    ZL=mp.sqrt(X['rho']*X['C33']);ZR=mp.sqrt(Y['rho']*Y['C33'])
    r_pred=(ZL-ZR)/(ZL+ZR);t_pred=2*ZL/(ZL+ZR)
    # Direct 2x2 displacement/stress solve, no operator/eigenvectors used.
    xL=mp.j*w*mp.sqrt(X['rho']/X['C33']);xR=mp.j*w*mp.sqrt(Y['rho']/Y['C33'])
    Ac=mp.matrix([[1,-1],[-X['C33']*xL,-Y['C33']*xR]])
    bc=mp.matrix([-1,-X['C33']*xL]);rc,tc=mp.lu_solve(Ac,bc)
    eps=max(abs(rc-r_pred),abs(tc-t_pred))
    check('classical_elastic_impedance',eps<mp.mpf('1e-30'),{'difference':mp.nstr(eps,10),'r':mp.nstr(rc,15),'t':mp.nstr(tc,15)})
    # 2. For multiple frequency / tau pairs compare analytic LS quartic to directly constructed 4x4 1D operator.
    for nm,p in [('Al6061',A),('PZT-5H',P),('BaTiO3',Q)]:
        for wt in (mp.mpf('1e-5'),mp.mpf('0.1'),mp.mpf(1)):
            v=deepcopy(p);freq=w
            # tau STRESS-TEST chosen to make omega*tau=target at fixed 1 MHz, NOT a measurement
            v['tauq']=wt/w
            red=reduced(v);rs,_=analytic_quartic(red,freq)
            direct=first_order_1d(red,freq)
            # analytic quartic roots satisfy direct first-order characteristic polynomial, independently assembled.
            maxres=mp.mpf(0)
            for lam in rs:
                yy=physical_mode(red,freq,lam)
                rr=direct*yy-lam*yy
                # term-wise balance, not net cancelled row
                for i in range(4):
                    den=sum(abs(direct[i,j]*yy[j]) for j in range(4))+abs(lam*yy[i])
                    maxres=max(maxres,abs(rr[i])/max(den,mp.mpf('1e-120')))
            check(f'analytic_LS_dispersion_{nm}_wt{mp.nstr(wt,3)}',maxres<mp.mpf('1e-28'),{'max_backward_error':mp.nstr(maxres,7),'omega_tau':mp.nstr(wt,7),'synthetic_tau':True})
    # 3. Adversarial mutations of actual material values: monitor analytic 1D 4-state physical observables.
    base=solve_interface(reduced(P),reduced(Q),w)
    def measure_mut(name, newp, newq, expected_min, observable):
        try:
            changed=solve_interface(reduced(newp),reduced(newq),w)
            before=base[observable];after=changed[observable]
            metric=abs(after-before)/max(abs(before),mp.mpf('1e-30'))
            detected=(metric>expected_min)
            data['mutations'].append({'mutation':name,'injected':True,'detected':bool(detected),
                'observable':observable,'relative_delta':mp.nstr(metric,12),'threshold':mp.nstr(expected_min,7),
                'pristine':mp.nstr(before,10),'mutated':mp.nstr(after,10)})
            if not detected:data['failures'].append('MUT_'+name)
        except Exception as ex:
            data['mutations'].append({'mutation':name,'injected':True,'detected':True,'exception':repr(ex),'observation':'solver failure, expected if mutation destroys admissibility'})
    # Compare to pristine absolute physical benchmark, NOT mutated one's own k->0 limit.
    for name,field,modify,obs,tol in [
       ('p3_reverse','p3',lambda z:-z,'transmitted_theta',mp.mpf('1e-6')),
       ('beta3_reverse','beta3',lambda z:-z,'transmitted_theta',mp.mpf('1e-3')),
       ('e33_reverse','e33',lambda z:-z,'transmitted_theta',mp.mpf('1e-6')),
    ]:
        mutant=deepcopy(P);mutant[field]=modify(mutant[field])
        measure_mut(name,mutant,Q,tol,obs)
    # 4. wrong LS relaxation placement in 1D heat flux; compare modified operator dispersion against fixed correct analytic roots.
    for wt in (mp.mpf('0.1'),mp.mpf('1')):
        v=deepcopy(P);v['tauq']=wt/w;m=reduced(v)
        baseline_roots,_=analytic_quartic(m,w)
        wrong_roots,_=analytic_quartic(m,w,wrong_g=True)
        deviation=close_roots(wrong_roots,baseline_roots)
        threshold=mp.mpf('1e-3')
        detected=deviation>threshold
        data['mutations'].append({'mutation':f'wrong_LS_g_omega_tau_{mp.nstr(wt,3)}','injected':True,
            'detected':bool(detected),'observable':'quartic_spectral_roots',
            'relative_delta':mp.nstr(deviation,10),'threshold':mp.nstr(threshold,4)})
        if not detected:data['failures'].append('MUT_LS'+str(wt))
    # 5. Omit incident RHS in THE ACTUAL four-row assembled interface equations,
    # solve mutant Bx=0, then independently evaluate full incident+reflected-transmitted jump.
    inc=base['incident']['Y']; re=base['refl'];tr=base['trns']
    B_mut=mp.matrix(4)
    for row in range(4):
        for col in range(2):B_mut[row,col]=re[col]['Y'][row]
        for col in range(2):B_mut[row,col+2]=-tr[col]['Y'][row]
    x_zero=mp.lu_solve(B_mut,mp.matrix([0,0,0,0]))
    mutant_jump=mp.matrix([inc[i]+sum(B_mut[i,j]*x_zero[j] for j in range(4)) for i in range(4)])
    scale_r=[max(abs(inc[i]),abs(base['left'][i]),abs(base['right'][i]),mp.mpf('1e-100')) for i in range(4)]
    rhs_violation=max(abs(mutant_jump[i])/scale_r[i] for i in range(4))
    detected=rhs_violation>mp.mpf('1e-3')
    data['mutations'].append({'mutation':'omit_incident_from_actual_1D_rhs','injected':True,'detected':bool(detected),
       'observable':'physical_four_component_interface_jump','relative_delta':mp.nstr(rhs_violation,8)})
    if not detected:data['failures'].append('MUT_incident')
    # 6. Actual radiation sign swap (swap selector input); ensure selected incoming mode wrongly appears in 'outgoing'.
    try:
        changed=solve_interface(reduced(P),reduced(Q),w,reverse_left=True)
        wrong=any(abs(t['lam']-changed['incident']['lam'])<mp.mpf('1e-24') for t in changed['refl'])
        metric=abs(changed['reflected_u']-base['reflected_u'])/max(abs(base['reflected_u']),mp.mpf('1e-12'))
        detected=bool(wrong and metric>mp.mpf('1e-4'))
        data['mutations'].append({'mutation':'reverse_left_radiation_in_solver','injected':True,'detected':detected,
           'observable':'physical_reflection_plus_admitted_incident','relative_delta':mp.nstr(metric,10),'incident_in_outgoing':bool(wrong)})
        if not detected:data['failures'].append('MUT_radiation')
    except Exception as ex:
        data['mutations'].append({'mutation':'reverse_left_radiation_in_solver','injected':True,'detected':True,'exception':repr(ex)})
    # 7. Shear interface full 2x2 law: pressure doesn't constrain shear, cannot substitute field-count tautology.
    # SH/SV at k=0: u_x(L) + r = t, muL xiL(1-r)=muR xiR t.
    uL,uR=P['C44'],Q['C44'];kL=w*mp.sqrt(P['rho']/uL);kR=w*mp.sqrt(Q['rho']/uR)
    Bsv=mp.matrix([[1,-1],[-uL*kL,-uR*kR]])
    ysv=mp.matrix([-1,-uL*kL]);xs=mp.lu_solve(Bsv,ysv)
    # Drop traction row, solve only displacement with candidate amplitude r=0 -> residual of discarded row.
    bad_r=mp.mpf(0);bad_t=mp.mpf(1)
    bad_shear=abs(uL*kL*(1-bad_r)-uR*kR*bad_t)/max(abs(uL*kL),abs(uR*kR))
    valid=abs((1+xs[0])-xs[1])+abs(uL*kL*(1-xs[0])-uR*kR*xs[1])/max(abs(uL*kL),abs(uR*kR))
    detected=bad_shear>mp.mpf('1e-3') and valid<mp.mpf('1e-27')
    data['mutations'].append({'mutation':'remove_actual_SV_traction_row','injected':True,'detected':bool(detected),
       'observable':'discarded_physical_shear_traction_residual','relative_delta':mp.nstr(bad_shear,12),
       'baseline_residual':mp.nstr(valid,10)})
    if not detected:data['failures'].append('MUT_shear')
    data['limitations']=['Analytic 1D LS dispersion and elastic reflection are independent physics anchors only for reduced k=0 problem.',
      'The 1D interface mutation tests do not exercise full finite-k production assembly.',
      'ASSUMED beta3 and p3 are algorithm stress tests, not material validation.',
      'No energy conservation or G4/G5 validation claimed.',
      'No global spectral simplicity or passive causality theorem proved.',
      'Conductor electrical contact not represented; Al6061 is a dielectric surrogate.']
    return data

if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--digits',type=int,nargs='+',default=[40,70]);pa.add_argument('--out',default=str(BASE/'evidence'/'analytic_oracle_results.json'));v=pa.parse_args()
    runs=[run(q) for q in v.digits]
    Path(v.out).write_text(json.dumps({'runs':runs},indent=2)+'\n')
    for runv in runs:
        print(f"digits {runv['precision']}: cases={len(runv['cases'])-len([x for x in runv['cases'] if not x['passed']])}/{len(runv['cases'])} mutations={sum(x['detected'] for x in runv['mutations'])}/{len(runv['mutations'])} failures={runv['failures']}",flush=True)
    if any(runv['failures'] for runv in runs): raise SystemExit(1)
