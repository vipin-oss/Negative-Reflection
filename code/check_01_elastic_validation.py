"""
check_01_elastic_validation.py  (v2)
====================================
Validation of the spectral plate solver against exact isotropic theory.

V1  Cutoff frequencies at k=0 vs analytic n*cL/(4h), n*cT/(4h).
V2  Dispersive frequencies vs the EXACT entire Rayleigh-Lamb function
    (checked both ways: spectral frequencies must be zeros of the exact
    function, and exact roots must be reproduced by the spectral solver).
V3  Mesh convergence of the S1 ZGV point; comparison with the published
    Duralumin scaling f_ZGV * 2h = 2.86 MHz*mm (Gerardin et al., Sci. Rep. 2019).
V4  Independent cross-check of the ZGV frequency with an exact-theory root
    tracker that never touches the spectral solver.

These are *solver* validations. None of them demonstrates negative reflection.

Run:  python3 code/check_01_elastic_validation.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_spectral import (Plate, isotropic_tensors, rayleigh_lamb_entire,
                           exact_cutoffs, rl_roots)   # noqa: E402

E, nu, rho = 71.7e9, 0.34, 2790.0
mu = E / (2 * (1 + nu))
lam = E * nu / ((1 + nu) * (1 - 2 * nu))
cL, cT = np.sqrt((lam + 2 * mu) / rho), np.sqrt(mu / rho)
h = 1.5e-3 / 2
C, e, eps, rho_m = isotropic_tensors(E, nu, rho)

print("=" * 78)
print("VALIDATION OF SPECTRAL PLATE SOLVER -- isotropic Duralumin-like material")
print(f"  E={E:.4g} Pa  nu={nu}  rho={rho} kg/m^3   cL={cL:.1f} m/s  cT={cT:.1f} m/s")
print(f"  plate total thickness 2h = {2*h*1e3:.2f} mm")
print("=" * 78)
res = {"material": {"E": E, "nu": nu, "rho": rho, "cL": cL, "cT": cT, "thickness_mm": 2 * h * 1e3}}

# ---------------------------------------------------------------- V1 cutoffs
pl = Plate(h, C, e, eps, rho, N=28)
f_spec = pl.spectrum(0.0)
ex = exact_cutoffs(cL, cT, h, nmax=7)
print("\n[V1] cutoff frequencies at k = 0")
print(f"{'exact f (kHz)':>15} {'spectral (kHz)':>16} {'rel.err':>11}")
errs = []
for fe in ex:
    fk = fe / 1e3
    j = int(np.argmin(np.abs(f_spec / 1e3 - fk)))
    fs_ = f_spec[j] / 1e3
    er = abs(fs_ - fk) / fk
    errs.append(er)
    print(f"{fk:15.4f} {fs_:16.4f} {er:11.1e}")
v1 = max(errs)
print(f"  max rel. error = {v1:.3e}  ->  {'PASS' if v1 < 1e-4 else 'FAIL'}")
res["V1_max_rel_err"] = v1

# ------------------------------------------------------- V2 vs exact function
print("\n[V2] dispersive frequencies vs exact Rayleigh-Lamb theory")
print("     (a) each spectral frequency must annihilate the exact dispersion")
print("         function: normalised residual |F| / (|Q*A| + |P*B|)   [0 = root]")
print("     (b) each exact root must be reproduced by the spectral solver")
worst_a, worst_b = 0.0, 0.0
FK = np.linspace(600.0, 2600.0, 11)
for k in FK:
    fsp = pl.spectrum(k)
    fsp = fsp[(fsp > 1.2e6) & (fsp < 8.0e6)]      # p, q both real here
    W = 2 * np.pi * fsp
    p_ = np.sqrt((W / cL) ** 2 - k ** 2 + 0j)
    q_ = np.sqrt((W / cT) ** 2 - k ** 2 + 0j)
    P = np.sin(p_ * h) * np.cos(q_ * h)
    Q = np.sin(q_ * h) * np.cos(p_ * h)
    A = (q_ ** 2 - k ** 2) ** 2
    B = 4 * k ** 2 * p_ * q_
    for sym in (True, False):
        F = (Q * A + P * B) if sym else (P * A + Q * B)
        norm = np.abs(Q * A) + np.abs(P * B) if sym else np.abs(P * A) + np.abs(Q * B)
        r = np.abs(F) / (norm + 1e-300)
        # each spectral mode is either symmetric or antisymmetric; take the better
        if sym:
            rbest = r.copy()
        else:
            rbest = np.minimum(rbest, r)
    worst_a = max(worst_a, float(np.max(rbest)))
    # (b) exact roots vs spectral list
    S, Aa = rl_roots(k, cL, cT, h, 1.2e6, 8.0e6)
    exact = np.sort(np.concatenate([S, Aa]))
    if len(exact) and len(fsp):
        d = np.min(np.abs(exact[:, None] - fsp[None, :]), axis=1)
        worst_b = max(worst_b, float(np.max(d / exact)))
print(f"  (a) worst normalised residual at spectral frequencies: {worst_a:.2e}")
print(f"  (b) worst relative root mismatch (exact vs spectral):  {worst_b:.2e}")
print(f"  -> {'PASS' if (worst_a < 1e-6 and worst_b < 1e-6) else 'FAIL'}")
res["V2a_norm_residual"] = worst_a
res["V2b_root_mismatch"] = worst_b

# ------------------------------------------------------ V3 ZGV convergence
print("\n[V3] S1 ZGV point: mesh convergence")
def zgv_spectral(N, k_lo=700.0, k_hi=2400.0, nk=300):
    p = Plate(h, C, e, eps, rho, N=N)
    ks = np.linspace(k_lo, k_hi, nk)
    spec = [p.spectrum(k) for k in ks]
    ref, traj = None, []
    for f in spec:
        f = f[(f > 1.4e6) & (f < 4.0e6)]
        if not len(f):
            traj.append(np.nan); continue
        if ref is None:
            ref = f[0]
        j = int(np.argmin(np.abs(f - ref)))
        ref = f[j]
        traj.append(ref)
    traj = np.array(traj)
    ok = np.isfinite(traj)
    ks_o, tr = ks[ok], traj[ok]
    d = np.diff(tr)
    cand = [i for i in range(1, len(d)) if d[i - 1] < 0 <= d[i]]
    if not cand:
        return np.nan, np.nan
    i = cand[0]
    y0, y1, y2 = tr[i - 1], tr[i], tr[i + 1]
    dd = (y0 - 2 * y1 + y2)
    dk = 0.5 * (y0 - y2) / dd * (ks_o[1] - ks_o[0]) if dd != 0 else 0.0
    return ks_o[i] + np.clip(dk, -1, 1) * (ks_o[1] - ks_o[0]), y1

conv = []
for N in (14, 18, 22, 26, 30, 34):
    kz, fz = zgv_spectral(N)
    conv.append((N, kz, fz))
    if np.isfinite(fz):
        print(f"  N={N:3d}   f_ZGV = {fz/1e6:8.5f} MHz   k_ZGV = {kz:9.2f} 1/m")
good = [c for c in conv if np.isfinite(c[2])]
if len(good) >= 2:
    rel = abs(good[-1][2] - good[-2][2]) / good[-1][2]
    print(f"  mesh-to-mesh relative change = {rel:.2e}")
    res["V3_zgv_f_MHz"] = good[-1][2] / 1e6
    res["V3_zgv_k"] = good[-1][1]
    res["V3_mesh_rel_change"] = float(rel)

# --------------------------------------- V4 independent exact-theory ZGV
print("\n[V4] independent ZGV from exact dispersion theory (no spectral solver)")
# Track the 2nd symmetric branch (S1) in k and locate its minimum (d f/dk = 0).
kks = np.linspace(700.0, 1600.0, 240)
traj = []
for k in kks:
    S, _ = rl_roots(k, cL, cT, h, 1.5e6, 2.6e6, n=6000)
    S = S[(S > 1.6e6) & (S < 2.5e6)]
    traj.append(S)
ref, series = None, []
for S in traj:
    if not len(S):
        series.append(np.nan); continue
    if ref is None:
        ref = S[0]
    j = int(np.argmin(np.abs(S - ref)))
    ref = S[j]
    series.append(ref)
series = np.array(series)
good = np.isfinite(series)
ks_o, tr = kks[good], series[good]
d = np.diff(tr)
mi = [i for i in range(1, len(d)) if d[i - 1] < 0 <= d[i]]
if mi:
    i = mi[0]
    y0, y1, y2 = tr[i - 1], tr[i], tr[i + 1]
    den = (y0 - 2 * y1 + y2)
    dk = 0.5 * (y0 - y2) / den * (ks_o[1] - ks_o[0]) if den else 0.0
    k_zgv = ks_o[i] + np.clip(dk, -1, 1) * (ks_o[1] - ks_o[0])
    print(f"  exact-theory S1 ZGV:  f = {y1/1e6:8.5f} MHz   k = {k_zgv:9.2f} 1/m")
    res["V4_zgv_f_MHz"] = y1 / 1e6
    res["V4_zgv_k"] = k_zgv
    if "V3_zgv_f_MHz" in res:
        agree = abs(res["V4_zgv_f_MHz"] - res["V3_zgv_f_MHz"]) / res["V4_zgv_f_MHz"]
        print(f"  spectral:            f = {res['V3_zgv_f_MHz']:8.5f} MHz   k = {res['V3_zgv_k']:9.2f} 1/m")
        print(f"  relative difference = {agree:.2e}  -> {'PASS' if agree < 5e-3 else 'FAIL'}")
        res["V4_vs_V3_rel_diff"] = float(agree)

# ------------------------------------------------ literature scaling check
print("\n[V5] comparison with published Duralumin scaling")
print("  Gerardin et al. (2019) report f_ZGV * 2h = 2.86 MHz.mm for Duralumin")
ref = 2.86e6 * 1e-3 / (2 * h)
print(f"  -> implies f_ZGV = {ref/1e6:.4f} MHz at 2h = {2*h*1e3:.2f} mm")
if "V3_zgv_f_MHz" in res:
    print(f"  computed        f_ZGV = {res['V3_zgv_f_MHz']:.4f} MHz")
    print(f"  deviation = {100*(res['V3_zgv_f_MHz']-ref/1e6)/(ref/1e6):+.1f} %")
    print("  NOTE: a residual difference is expected because the published value")
    print("        is for the plate they used and because the ZGV frequency is")
    print("        strongly Poisson-ratio dependent; the check is order-of-magnitude.")
    res["V5_lit_ref_MHz"] = ref / 1e6
    res["V5_rel_dev"] = float((res["V3_zgv_f_MHz"] - ref / 1e6) / (ref / 1e6))

with open("data/check_01_elastic_validation.json", "w") as f:
    json.dump(res, f, indent=2)
print("\nwrote data/check_01_elastic_validation.json")
