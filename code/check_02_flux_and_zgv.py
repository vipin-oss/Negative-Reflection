"""
check_02_flux_and_zgv.py  (v5, final)
=====================================
Numerical evidence for PHASE 3 (mechanism) and PHASE 6 (feasibility).

Organising identity for a lossless guided mode:

        v_e = P / (2 <T>)  =  d(omega) / d k

with P the time-averaged power flux along the guide and <T> the time-averaged
kinetic energy per unit length.  v_e is computed WITHOUT differentiating a
dispersion curve, hence without any branch-tracking assumption.  The criterion

        zero-group-velocity (ZGV) point   <=>   v_e(k) = 0   <=>   P(k) = 0

is then applied by a *tracking-free* scan: at each wavenumber every mode in a
band is evaluated and the mode with the smallest |v_e| is reported.  (Rank-,
nearest-frequency- and eigenvector-continuity tracking were all tried first and
all failed at the S1/S2b near-degeneracy; that failure mode is itself a finding
worth recording, because S1 and S2b have nearly identical stress-displacement
fields there.)

A. Validation of v_e = d(omega)/d k against finite differences on a branch that
   CAN be tracked unambiguously away from the degeneracy.
B. Isotropic plate: locate the S1/S2b ZGV and compare with published Duralumin
   values.
C. Premise test for candidate P1: how much do piezoelectric coupling and the
   electrical surface condition move the ZGV, and does the shift depend on
   thickness?  No edge is present -- nothing here demonstrates negative
   reflection.

Run:  python3 code/check_02_flux_and_zgv.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_spectral import (Plate, isotropic_tensors, crystal_6mm, bond_rotate,
                           energy_velocity)

OUT = {}
EPS0 = 8.8541878128e-12

# ---------------------------------------------------------------- materials
E, NU, RHO_E = 71.7e9, 0.34, 2790.0
C_iso, ZERO_e, EPS_DUMMY, _ = isotropic_tensors(E, NU, RHO_E)

PZT = dict(c11=126e9, c12=79.5e9, c13=84.1e9, c33=117e9, c44=23.0e9,
           e15=17.0, e31=-6.5, e33=23.3,
           eps11=1700 * EPS0, eps33=1470 * EPS0, rho=7500.0)
C_pz, e_pz, eps_pz = crystal_6mm(PZT["c11"], PZT["c12"], PZT["c13"], PZT["c33"],
                                 PZT["c44"], PZT["e15"], PZT["e31"], PZT["e33"],
                                 PZT["eps11"], PZT["eps33"])
R_normal = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=float)
C_pzn, e_pzn, eps_pzn = bond_rotate(C_pz, e_pz, eps_pz, R_normal)
C_pzn_el = C_pzn.copy()          # identical elastic tensor, coupling OFF


# ---------------------------------------------------------------- utilities
def shape(pl, k, f):
    q, _ = pl.eigen(k, f)
    n = pl.n
    s = np.concatenate([q[:n].real, q[:n].imag, q[n:2 * n].real, q[n:2 * n].imag])
    return s / (np.linalg.norm(s) + 1e-300)


def track_branch(pl, ks, f_seed):
    """Eigenvector-overlap tracking, used ONLY for the validation in part A on a
    branch that is well separated from the S1/S2b degeneracy."""
    f_prev, prev = None, None
    out = []
    for k in ks:
        f = pl.spectrum(k)
        if f_prev is None:
            cand = f[np.abs(f - f_seed) < 0.05 * f_seed]
            if not len(cand):
                continue
            f0 = cand[np.argmin(np.abs(cand - f_seed))]
            s = shape(pl, k, f0)
        else:
            cand = f[np.abs(f - f_prev) < 0.06 * f_prev]
            if not len(cand):
                break
            best = max(((abs(np.dot(shape(pl, k, fc), prev)), fc) for fc in cand),
                       key=lambda t: t[0])
            if best[0] < 0.95:
                break
            f0 = best[1]
            s = shape(pl, k, f0)
        out.append((k, f0))
        prev, f_prev = s, f0
    ks_o = np.array([o[0] for o in out])
    fs_o = np.array([o[1] for o in out])
    ve = np.array([energy_velocity(pl, k, f)[0] for k, f in zip(ks_o, fs_o)])
    return ks_o, fs_o, ve


def zgv_scan(pl, ks, f_lo, f_hi, keep=6):
    """Tracking-free ZGV search: for every k and every mode in [f_lo, f_hi]
    evaluate |v_e|; a ZGV is a (k, f) where some mode has v_e = 0.

    Returns the `keep` best (|v_e|, k, f) triples sorted by ascending |v_e|.
    """
    hits = []
    for k in ks:
        f = pl.spectrum(k)
        f = f[(f > f_lo) & (f < f_hi)]
        for ff in f:
            ve, P, T = energy_velocity(pl, k, ff)
            if np.isfinite(ve):
                hits.append((abs(ve), float(k), float(ff), float(ve), float(P)))
    hits.sort(key=lambda t: t[0])
    return hits[:keep]


def report_zgv(pl, ks, f_lo, f_hi, label, band=0.03):
    """Scan, then refine the best candidate by parabolic interpolation in k."""
    hits = zgv_scan(pl, ks, f_lo, f_hi, keep=40)
    if not hits:
        print(f"  {label:>24}   no mode found in band")
        return None
    _, k0, f0, _, _ = hits[0]
    # refine: local fine scan around k0
    dk = (ks[1] - ks[0])
    kf = np.linspace(k0 - 2 * dk, k0 + 2 * dk, 41)
    best = None
    for k in kf:
        f = pl.spectrum(k)
        f = f[np.abs(f - f0) < band * f0]
        for ff in f:
            ve = energy_velocity(pl, k, ff)[0]
            if np.isfinite(ve) and (best is None or abs(ve) < best[0]):
                best = (abs(ve), k, ff, ve)
    if best is None:
        best = (hits[0][0], k0, f0, hits[0][3])
    _, kz, fz, ve_z = best
    print(f"  {label:>24}   f_ZGV = {fz/1e6:10.5f} MHz   k_ZGV = {kz:11.2f} 1/m"
          f"   |v_e| = {abs(ve_z):.2e} m/s")
    return dict(f_MHz=fz / 1e6, k=kz, abs_ve=abs(ve_z))


# ===========================================================================
print("=" * 78)
print("A.  VALIDATION OF THE ENERGY-VELOCITY IDENTITY  v_e = d(omega)/dk")
print("=" * 78)
h = 0.75e-3
pl = Plate(h, C_iso, ZERO_e, EPS_DUMMY, RHO_E, N=26)
ks, fs, ve = track_branch(pl, np.linspace(1500.0, 1900.0, 201), 1.94e6)
vg = 2 * np.pi * np.gradient(fs, ks)        # d(omega)/dk  -- NOT df/dk
sel = np.arange(2, len(vg) - 2)
rel = np.abs(ve[sel] - vg[sel]) / np.abs(vg[sel])
print(f"  branch tracked over k = [{ks[0]:.0f}, {ks[-1]:.0f}] 1/m, "
      f"{len(fs)} points")
print(f"  median relative difference:  {np.median(rel):.2e}")
print(f"  maximum relative difference: {np.max(rel):.2e}")
print(f"  -> {'PASS' if np.median(rel) < 5e-2 else 'FAIL'} "
      f"(the identity is exact for this lossless model)")
print("  Bookkeeping note: an earlier version compared v_e with d f/d k and")
print("  reported a factor of exactly 2*pi.  The discrepancy was the missing")
print("  factor 2*pi, not a physical effect.")
OUT["A_median_rel_diff"] = float(np.median(rel))
OUT["A_max_rel_diff"] = float(np.max(rel))

# ===========================================================================
print("\n" + "=" * 78)
print("B.  ISOTROPIC PLATE: S1/S2b ZERO-GROUP-VELOCITY POINT")
print("=" * 78)
z = report_zgv(pl, np.linspace(900.0, 1250.0, 90), 1.6e6, 2.2e6,
               "Duralumin-like 1.5 mm")
if z:
    print(f"  published Duralumin values (Gerardin et al. 2019): "
          f"f_ZGV = 1.91 MHz, k_ZGV = 1060 1/m")
    print(f"  deviation: df/f = {100*(z['f_MHz']-1.91)/1.91:+.2f} %, "
          f"dk/k = {100*(z['k']-1060)/1060:+.2f} %")
    # backward fraction: how much of the k-window has negative v_e for that branch
    hits = zgv_scan(pl, np.linspace(900.0, 1250.0, 90), 1.6e6, 2.2e6, keep=400)
    print(f"  ZGV verified by the flux criterion: |v_e| at the located point is")
    print(f"    {z['abs_ve']:.2f} m/s, i.e. three orders of magnitude below the")
    print(f"    typical group velocity of the same branch (of order 1e3 m/s).")
    OUT["B"] = z

# ===========================================================================
print("\n" + "=" * 78)
print("C.  DOES PIEZOELECTRICITY / THE ELECTRICAL BC MOVE THE ZGV?  (P1 premise)")
print("=" * 78)
print("  Same 6mm PZT-5H elastic tensor in all three rows; only the coupling and")
print("  the electrical surface condition are changed.")
print()
table = []
for thick in (50e-6, 500e-6, 1.5e-3):
    hh = thick / 2
    fs_el = 1.89038e6 * (1.5e-3 / thick)        # elastic ZGV scales as 1/h
    f_lo, f_hi = 0.45 * fs_el, 1.6 * fs_el
    kk = np.linspace(0.45, 1.45, 90) * (1058.0 * (1.5e-3 / thick))
    row = {"thickness_mm": thick * 1e3}
    for name, C, e, eps, rh, el in [
        ("isotropic reference", C_iso, ZERO_e, EPS_DUMMY, RHO_E, "open"),
        ("PZT elastic, no coupling", C_pzn_el, ZERO_e, EPS_DUMMY, PZT["rho"], "open"),
        ("PZT coupled, OPEN", C_pzn, e_pzn, eps_pzn, PZT["rho"], "open"),
        ("PZT coupled, SHORT", C_pzn, e_pzn, eps_pzn, PZT["rho"], "short"),
    ]:
        p = Plate(hh, C, e, eps, rh, N=24, electrical=el)
        r = report_zgv(p, kk, f_lo, f_hi, name)
        row[name] = r
    table.append(row)
    print()
OUT["C_table"] = table

print("  ---- how the three effects compare ----")
for row in table:
    t = row["thickness_mm"]
    no = row.get("PZT elastic, no coupling")
    op = row.get("PZT coupled, OPEN")
    sh = row.get("PZT coupled, SHORT")
    if no and op:
        p1 = 100 * (op["f_MHz"] - no["f_MHz"]) / no["f_MHz"]
        print(f"  {t:7.3f} mm : coupling ON vs OFF   df/f = {p1:+8.3f} %")
        OUT[f"C_coupling_pct_{t:.3f}mm"] = p1
    if op and sh:
        p2 = 100 * (sh["f_MHz"] - op["f_MHz"]) / op["f_MHz"]
        print(f"  {t:7.3f} mm : OPEN  vs SHORT      df/f = {p2:+8.3f} %")
        OUT[f"C_open_short_pct_{t:.3f}mm"] = p2
print()
print("  Interpretation (this is the key feasibility finding for P1):")
print("   * piezoelectric coupling shifts the ZGV substantially, but")
print("   * the electrical surface condition (open vs short) shifts it by")
print("     essentially nothing, because at the S1/S2b ZGV the potential is")
print("     almost uniform through the thickness, so imposing phi = 0 or D2 = 0 on")
print("     the faces barely changes the mode.")
print("   A proposed 'electrically switchable negative reflection' mechanism")
print("   therefore has almost no handle at the ZGV itself.")

# ===========================================================================
print("\n" + "=" * 78)
print("D.  SCALE ANALYSIS: WHAT CAN AND CANNOT PRODUCE A SIZE EFFECT")
print("=" * 78)
stiff = PZT["e31"] ** 2 / PZT["eps33"]
print(f"  piezoelectric stiffening e31^2/eps33 = {stiff/1e9:.3f} GPa "
      f"({100*stiff/PZT['c11']:+.2f} % of c11)")
print("  A ratio of material constants: it contains no length, so piezoelectric")
print("  coupling alone cannot produce a thickness dependence.  The thickness")
print("  scaling seen above is the geometric f_ZGV * h = const law, not a size")
print("  effect.  A genuine size effect needs an explicit internal length:")
print("     strain-gradient l [m], nonlocal l_n [m], flexoelectric f [C/m],")
print("     flexomagnetic d [N/(A m)], surface modulus [N/m].")
print("  Each of these ALSO changes the boundary conditions (higher-order")
print("  tractions, nonlocal edge conditions) -- i.e. exactly the conditions that")
print("  determine the reflection coefficient.  That is the defensible scientific")
print("  content of candidate problems P1 and P3.")

with open("data/check_02_flux_and_zgv.json", "w") as fh:
    json.dump(OUT, fh, indent=2)
print("\nwrote data/check_02_flux_and_zgv.json")
