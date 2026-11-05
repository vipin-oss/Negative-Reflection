"""Independent full-8-state Fourier-z first-order LS piezothermoelastic operator.

NO imports of repository source / scattering code, no use of cached roots.
Conventions: exp(ik*x-i*omega*t), y'=A*y, y=(ux,uz,phi,theta,txz,tzz,Dz,qz).
Material parameters read as literal decimal strings from the source CSV.
Scientific scope: constitutive/operator verification, k->0 electrostatic-sector diagnostics.
This is not a validated scattering solver unless radiation and conditioning independently pass.
"""
import csv, argparse, json
from pathlib import Path
import mpmath as mp

KEYS=('C11','C13','C33','C44','e31','e33','e15','eps11','eps33',
      'beta1','beta3','p3','rho','cE','kappa11','kappa33','tauq','T0')
STATE=('ux','uz','phi','theta','txz','tzz','Dz','qz')

def load(path, name):
    v={}
    with open(path,encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f):
            if row['material']==name and row['symbol'] in KEYS:
                v[row['symbol']]=mp.mpf(row['value'])
    if set(v)!=set(KEYS): raise ValueError(f'{name} missing {set(KEYS)-set(v)}')
    return v

def operator(p,k,w):
    """Evaluate derivative y' by linearity on 8 basis vectors, not characteristic roots.

    From txz=C44(ux'+ik uz)+i k e15 phi,
      tzz=i k C13 ux+C33 uz'+e33 phi'-beta3 theta,
      Dz=i k e31 ux+e33 uz'-eps33 phi'+p theta,
    plus momentum, Gauss and Fourier Lord-Shulman entropy balance.
    """
    I=mp.j;g=1-I*w*p['tauq']; C=p['C33']; e=p['e33'];eps=p['eps33']
    det=C*eps+e**2
    def deriv(y):
        ux,uz,ph,th,txz,tzz,Dz,qz=y
        # solve sigma_zz and D_z equations independently of repository M
        sig=tzz-I*k*p['C13']*ux+p['beta3']*th
        dis=Dz-I*k*p['e31']*ux-p['p3']*th
        uzp=(eps*sig+e*dis)/det
        php=(e*sig-C*dis)/det
        uxp=txz/p['C44']-I*k*uz-I*k*p['e15']*ph/p['C44']
        thp=-g*qz/p['kappa33']
        sxx=I*k*p['C11']*ux+p['C13']*uzp+p['e31']*php-p['beta1']*th
        Dx=p['e15']*(uxp+I*k*uz)-I*k*p['eps11']*ph
        txzp=-I*k*sxx-p['rho']*w*w*ux
        tzzp=-I*k*txz-p['rho']*w*w*uz
        Dzp=-I*k*Dx
        # Heat: -(kappa11 k² +kappa33 xi²) th +iωg[rho*c th +T(β1 ik ux+β3 uz' -p phi')]=0
        qzp=-p['kappa11']*k*k*th/g +I*w*(p['rho']*p['cE']*th+p['T0']*(I*k*p['beta1']*ux+p['beta3']*uzp-p['p3']*php))
        return mp.matrix([uxp,uzp,php,thp,txzp,tzzp,Dzp,qzp])
    A=mp.matrix(8)
    for j in range(8):
        v=[mp.mpc(0)]*8;v[j]=1
        r=deriv(v)
        for i in range(8):A[i,j]=r[i]
    return A

def scale(p,w):
    kr=w/mp.sqrt(p['C33']/p['rho']);
    phi=max(abs(p['e33'])/p['eps33'],mp.mpf(1))
    theta=max(abs(p['beta3'])*kr/(p['rho']*p['cE']),mp.mpf(1))
    scales=[mp.mpf(1),mp.mpf(1),phi,theta,p['C44']*kr,p['C33']*kr,
            max(p['eps33']*phi*kr,mp.mpf('1e-9')),p['kappa33']*theta*kr]
    return kr,scales

def balanced(A,p,w):
    kr,s=scale(p,w)
    B=mp.matrix(8)
    for i in range(8):
        for j in range(8):B[i,j]=A[i,j]*s[j]/(s[i]*kr)
    return B,kr,s

def residual(A,lam,y):
    r=A*y-lam*y
    num=max(abs(t) for t in r)
    normy=max(abs(t) for t in y)
    normA=max(sum(abs(A[i,j]) for j in range(A.cols)) for i in range(A.rows))
    return num/max(normA*normy+abs(lam)*normy,mp.mpf('1e-150'))

def run(paths,out):
    summary=[]
    for dps in (40,70):
      mp.mp.dps=dps
      w=2*mp.pi*mp.mpf('1e6')
      for name in ('Al6061','PZT-5H'):
       p=load(paths,name)
       B0,kr,s=balanced(operator(p,mp.mpf(0),w),p,w)
       # Verify the exactly-zero constant potential eigenvector and generalized Jordan rank qualitatively
       phi=mp.matrix([0,0,1,0,0,0,0,0])
       zres=max(abs(t) for t in B0*phi)
       for ratio in (mp.mpf('1e-2'),mp.mpf('1e-4')):
        k=kr*ratio;A=operator(p,k,w);B,_,_=balanced(A,p,w)
        # eigenpairs of an independent 8x8 operator; spectral eigenvalues measured in k_ref
        vals,V=mp.eig(B)
        errs=[]
        for j,la in enumerate(vals):
            y=mp.matrix([V[i,j] for i in range(8)])
            errs.append(residual(B,la,y))
        # inspect two roots nearest zero: gauge pair approx k
        order=sorted(range(8),key=lambda j:abs(vals[j]))
        low=[vals[j] for j in order[:2]]
        entry={'material':name,'digits':dps,'k_over_kref':mp.nstr(ratio,10),
         'max_eigen_residual':mp.nstr(max(errs),12),
         'small_roots_over_kref':[mp.nstr(v,18) for v in low],
         'small_roots_divided_by_k':[mp.nstr(v/ratio,18) for v in low],
         'zero_constant_phi_residual':mp.nstr(zres,12)}
        summary.append(entry);print(name,dps,mp.nstr(ratio,3),'eigerr',mp.nstr(max(errs),3),'roots',[mp.nstr(v/ratio,5) for v in low],flush=True)
    Path(out).write_text(json.dumps({'scope':'independent 8-state bulk spectral diagnostics, NOT scattering/LAP','entries':summary},indent=2)+'\n')
    return summary

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--materials',required=True);a.add_argument('--out',required=True);o=a.parse_args()
    run(o.materials,o.out)
