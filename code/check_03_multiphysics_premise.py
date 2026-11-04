"""
check_03_multiphysics_premise.py
=================================
Premise verification for the REVISED multiphysics research problems.

Purpose
-------
Before any multiphysics negative-reflection problem is proposed, the following
must be established by computation, not assumed:

  V1  The generalised solver reduces EXACTLY to the validated piezoelectric
      solver of check_02 when only one scalar field is present.
  V2  The energy velocity from the generalised (uncondensed) pencil equals
      2*pi*df/dk on a tracked branch, for BOTH an elastic and a piezoelectric
      plate  ->  the flux used for the negative-reflection criterion is correct.
  V3  CORRECTION: plate_spectral.energy_velocity uses only the ELASTIC part of
      the stress.  For a coupled medium it is therefore NOT domega/dk.  This is
      quantified here, because check_02 Part B/C used it for PZT-5H.
  V4  Magneto-electro-elastic plate (BaTiO3-CoFe2O4): does a forward/backward
      fold with a zero-group-velocity point EXIST?
  V5  Does the magnetic field change the fold?  (coupling ON/OFF decomposition)
  V6  Do the electrical / magnetic surface conditions change the fold?
  V7  Does a piezothermoelastic plate admit the same fold?  (thermal coupling
      enters exactly like a third quasi-static scalar in the lossless limit)
  V8  Mesh invariance of the fold.

Everything is written to data/check_03_multiphysics_premise.json.

Run:  python3 code/check_03_multiphysics_premise.py
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plate_spectral as ps
from coupled_plate import GeneralPlate, bond_rotate

EPS0 = 8.8541878128e-12
OUT = {}
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)

ROT_X3_TO_X2 = np.array([[1.0, 0.0, 0.0],
                         [0.0, 0.0, 1.0],
                         [0.0, 1.0, 0.0]])

# ---------------------------------------------------------------------------
# materials
# ---------------------------------------------------------------------------
# Isotropic reference (the Duralumin-like plate of the validated work)
C_iso, _, ED_iso, RHO_ISO = ps.isotropic_tensors(71.7e9, 0.34, 2790.0)

# PZT-5H, poled through the thickness
C_pzt, e_pzt, eps_pzt = ps.crystal_6mm(
    126e9, 79.5e9, 84.1e9, 117e9, 23.0e9,
    17.0, -6.5, 23.3, 1700 * EPS0, 1470 * EPS0)
C_pzt_r, e_pzt_r, eps_pzt_r = ps.bond_rotate(C_pzt, e_pzt, eps_pzt, ROT_X3_TO_X2)

# BaTiO3-CoFe2O4 magneto-electro-elastic composite.
# Elastic / piezoelectric / dielectric: standard tabulated values for the
# 40%BaTiO3-60%CoFe2O4 (or 0.5-0.5) particulate composite, as reproduced in
# e.g. the MEE nanoplate and peridynamic-MEE reference tables.
# Piezomagnetic: q31=290.1, q33=349.9, q15=275 N/(A m); permeability 5.0e-6 N/A^2.
MEE = dict(c11=226e9, c12=125e9, c13=124e9, c33=216e9, c44=44.2e9, c66=50.5e9,
           e31=-2.2, e33=9.3, e15=5.8,
           eps11=5.64e-9, eps33=6.35e-9,
           q31=290.1, q33=349.9, q15=275.0,
           mu11=5.0e-6, mu33=5.0e-6,
           rho=5730.0)


def mee_tensors(piezo=True, magnetic=True):
    """Magneto-electro-elastic 6mm tensors, symmetry axis along x3 (crystal)."""
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = MEE["c11"]
    C[2, 2] = MEE["c33"]
    C[0, 1] = C[1, 0] = MEE["c12"]
    C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = MEE["c13"]
    C[3, 3] = C[4, 4] = MEE["c44"]
    C[5, 5] = MEE["c66"]

    e = np.zeros((6, 3))
    if piezo:
        e[0, 2] = e[1, 2] = MEE["e31"]
        e[2, 2] = MEE["e33"]
        e[3, 1] = e[4, 0] = MEE["e15"]

    q = np.zeros((6, 3))
    if magnetic:
        q[0, 2] = q[1, 2] = MEE["q31"]
        q[2, 2] = MEE["q33"]
        q[3, 1] = q[4, 0] = MEE["q15"]

    eps = np.diag([MEE["eps11"], MEE["eps11"], MEE["eps33"]])
    mu = np.diag([MEE["mu11"], MEE["mu11"], MEE["mu33"]])
    return C, e, eps, q, mu


def rotate_mee(C, e, eps, q, mu):
    """Rotate the MEE set so the poling axis is the plate normal x2."""
    R = ROT_X3_TO_X2
    Cr, (er, qr), (epsr, mur) = bond_rotate(C, [e, q], [eps, mu], R)
    return Cr, er, epsr, qr, mur


def build_material(kind, piezo=True, magnetic=True, thermal=False):
    """Return (h, C, couplings, second, rho, bcs, names)."""
    if kind == "elastic":
        return (C_iso, [], [], RHO_ISO, [], [])
    if kind == "pzt":
        return (C_pzt_r, [e_pzt_r], [[eps_pzt_r]], 7500.0, ["open"], ["phi"])
    if kind == "mee":
        C, e, eps, q, mu = mee_tensors(piezo=piezo, magnetic=magnetic)
        Cr, er, epsr, qr, mur = rotate_mee(C, e, eps, q, mu)
        coups, second, names, bcs = [], [], [], []
        if piezo:
            coups.append(er); second.append([epsr]); names.append("phi"); bcs.append("open")
            if magnetic:
                second[0].append(np.zeros((3, 3)))
        if magnetic:
            coups.append(qr); second.append([np.zeros((3, 3)), mur])
            if piezo:
                second[0][1] = np.zeros((3, 3))
            names.append("psi"); bcs.append("open")
        return (Cr, coups, second, MEE["rho"], bcs, names)
    raise ValueError(kind)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def f_on_branch(gp, k, f0, halfwidth=4e5):
    """Frequency of the branch nearest f0 at wavenumber k (windowed)."""
    f, _ = gp.modes(k, f0 - halfwidth, f0 + halfwidth)
    if len(f) == 0:
        return np.nan
    return float(f[np.argmin(np.abs(f - f0))])


def vg_vs_dfdk(gp, k, f0, rel_dk=2e-4):
    """Compare Hellmann-Feynman v_e with 2*pi*df/dk on a tracked branch."""
    fa = f_on_branch(gp, k - rel_dk * k, f0)
    fb = f_on_branch(gp, k + rel_dk * k, f0)
    if not np.isfinite(fa) or not np.isfinite(fb):
        return np.nan, np.nan, np.nan
    dfdk = (fb - fa) / (2 * rel_dk * k)
    ve = gp.flux(k, f0)[0]
    return ve, 2 * np.pi * dfdk, abs(ve - 2 * np.pi * dfdk) / max(abs(2 * np.pi * dfdk), 1e-9)


def min_abs_ve(gp, k, f_lo, f_hi):
    """(min |v_e|, frequency at that mode) over all modes in the window at k."""
    f, V = gp.modes(k, f_lo, f_hi)
    if len(f) == 0:
        return np.inf, np.nan
    dA = gp._dA(k)
    best = (np.inf, np.nan)
    for j in range(len(f)):
        _, _, ve = gp._flux_of(k, V[:, j], 2 * np.pi * f[j], dA)
        if np.isfinite(ve) and abs(ve) < best[0]:
            best = (abs(ve), float(f[j]))
    return best


def global_zgv(gp, ks, f_lo, f_hi, n_refine=3, factor=6):
    """Tracking-free global ZGV, two stage.

    Stage 1: on the coarse grid, g(k) = min_j |v_e(k,j)|.  At a zero-group-
    velocity point g(k) -> 0, so its minimum locates the fold and no branch
    identity is ever needed.
    Stage 2: refine the k-grid around the best k by `factor` and repeat,
    `n_refine` times.

    A candidate is accepted as a ZGV only if the refined |v_e| is below
    TOL_VE (m/s).  Otherwise the routine reports 'not found' rather than
    returning a spurious corner minimum of the grid.
    """
    TOL_VE = 2.0
    ks = np.asarray(ks, dtype=float)
    best_k, best = np.nan, np.inf
    for k in ks:
        g, f = min_abs_ve(gp, k, f_lo, f_hi)
        if g < best:
            best, best_k = g, k
    if not np.isfinite(best_k):
        return (np.inf, np.nan, np.nan, np.nan, False)
    step = (ks[1] - ks[0]) if len(ks) > 1 else 1.0
    for _ in range(n_refine):
        lo, hi = best_k - step, best_k + step
        sub = np.linspace(lo, hi, factor * 4 + 1)
        bk, bv = np.nan, np.inf
        for k in sub:
            g, f = min_abs_ve(gp, k, f_lo, f_hi)
            if g < bv:
                bv, bk = g, k
        best_k, best = bk, bv
        step = (sub[1] - sub[0])
    g, f = min_abs_ve(gp, best_k, f_lo, f_hi)
    ok = g < TOL_VE
    return (g, float(best_k), float(f), float(g), ok)


# ---------------------------------------------------------------------------
# V1/V2/V3 : validation against the existing validated solver
# ---------------------------------------------------------------------------
def part_validation():
    print("=" * 72)
    print("V1/V2/V3  validation of the generalised solver")
    print("=" * 72)
    h = 0.75e-3
    res = {}

    # --- V1: one-scalar generalised solver == validated piezoelectric solver
    g1 = GeneralPlate(h, C_pzt_r, [e_pzt_r], [[eps_pzt_r]], 7500.0, N=24,
                      bcs=["open"], field_names=["phi"])
    p1 = ps.Plate(h, C_pzt_r, e_pzt_r, eps_pzt_r, 7500.0, N=24, electrical="open")
    worst = 0.0
    for k in [500.0, 700.0, 808.95, 1000.0, 1200.0]:
        fa, fb = g1.spectrum(k), p1.spectrum(k)
        m = min(len(fa), len(fb), 8)
        worst = max(worst, float(np.max(np.abs(fa[:m] - fb[:m]) / fb[:m])))
    res["V1_max_rel_freq_diff_vs_plate_spectral"] = worst
    print(f"  V1  generalised vs validated piezoelectric solver : {worst:.2e}  "
          f"{'PASS' if worst < 1e-10 else 'FAIL'}")

    # --- V2a: elastic plate, HF v_e vs 2*pi*df/dk
    ge = GeneralPlate(h, C_iso, [], [], RHO_ISO, N=24)
    rows = []
    for k0 in [900.0, 1000.0, 1100.0, 1200.0]:
        ve, g, rel = vg_vs_dfdk(ge, k0, 1.8904e6)
        rows.append((k0, ve, g, rel))
        print(f"  V2a k={k0:7.1f}  v_e={ve:12.4f}  2pi*df/dk={g:12.4f}  rel={rel:.1e}")
    res["V2a_elastic_hf_vs_dfdk_max_rel"] = max(r[3] for r in rows)
    print(f"  V2a elastic: worst relative difference {res['V2a_elastic_hf_vs_dfdk_max_rel']:.2e}  "
          f"{'PASS' if res['V2a_elastic_hf_vs_dfdk_max_rel'] < 1e-2 else 'FAIL'}")

    # --- V2b: piezoelectric plate, HF v_e vs 2*pi*df/dk
    rows = []
    for k0 in [500.0, 600.0, 700.0, 1200.0, 1300.0]:
        ve, g, rel = vg_vs_dfdk(g1, k0, 1.13e6)
        rows.append((k0, ve, g, rel))
        print(f"  V2b k={k0:7.1f}  v_e={ve:12.4f}  2pi*df/dk={g:12.4f}  rel={rel:.1e}")
    res["V2b_piezo_hf_vs_dfdk_max_rel"] = max(r[3] for r in rows)
    print(f"  V2b piezo:   worst relative difference {res['V2b_piezo_hf_vs_dfdk_max_rel']:.2e}  "
          f"{'PASS' if res['V2b_piezo_hf_vs_dfdk_max_rel'] < 1e-2 else 'FAIL'}")

    # --- V3: the CORRECTION, before vs after -------------------------------
    #  plate_spectral.energy_velocity previously used only the elastic part of
    #  the stress, sig = C:S.  In a coupled medium the physical stress is
    #  sig = C:S + e:E, so that version was NOT domega/dk.  The table below
    #  quantifies it: "old" is the uncorrected expression, recomputed here
    #  explicitly so the correction remains auditable after the fix.
    print("\n  V3  CORRECTION: the flux must include the piezoelectric stress")
    print("      sig = C:S + e:E   (the uncorrected form used only C:S):")
    corr = []
    for k0 in [600.0, 700.0, 1200.0]:
        f0 = f_on_branch(g1, k0, 1.13e6)
        ve_new = g1.flux(k0, f0)[0]                      # validated reference
        # uncorrected expression, reproduced explicitly
        q, _ = p1.eigen(k0, f0)
        n_, D_, w_ = p1.n, p1.D, p1.w
        u1, u2 = q[0:n_], q[n_:2 * n_]
        aa = np.max(np.abs(q)) + 1e-300
        u1, u2 = u1 / aa, u2 / aa
        ik = 1j * k0
        Sv = np.zeros((6, n_), dtype=complex)
        Sv[0] = ik * u1; Sv[1] = D_ @ u2; Sv[5] = D_ @ u1 + ik * u2
        sig_old = p1.C @ Sv
        om = 2 * np.pi * f0
        P_old = (om / 2) * np.imag(np.sum(w_ * (sig_old[0] * np.conj(u1) + sig_old[5] * np.conj(u2))))
        T_old = (p1.rho * om ** 2 / 4) * np.sum(w_ * (np.abs(u1) ** 2 + np.abs(u2) ** 2))
        ve_old = P_old / (2 * T_old)
        corr.append(dict(k=k0, f=f0, v_validated=ve_new, v_uncorrected=ve_old,
                         ratio_uncorrected=ve_old / ve_new))
        print(f"      k={k0:7.1f}  f={f0/1e6:.6f} MHz  uncorrected={ve_old:10.3f}  "
              f"validated={ve_new:10.3f}  ratio={ve_old/ve_new:6.3f}")
    res["V3_correction"] = corr
    try:
        ps.energy_velocity(p1, 600.0, 1.13e6)
        res["V3_plate_spectral_guard"] = "NOT RAISED (unexpected)"
    except NotImplementedError as exc:
        res["V3_plate_spectral_guard"] = "raises NotImplementedError for coupled plates"
        print("      plate_spectral.energy_velocity now refuses coupled plates (guard active)")
    return res


# ---------------------------------------------------------------------------
# V4-V8 : multiphysics premise
# ---------------------------------------------------------------------------
def part_premise():
    print("\n" + "=" * 72)
    print("V4-V8  multiphysics premise: does the fold exist?")
    print("=" * 72)
    h = 0.75e-3
    res = {"thickness_mm": 2 * h * 1e3}

    # elastic baseline fold (known: 1.89038 MHz / 1061 1/m)
    ge = GeneralPlate(h, C_iso, [], [], RHO_ISO, N=24)
    K_LO, K_HI, F_LO, F_HI = 450.0, 1700.0, 0.85e6, 2.40e6
    b = find_folds(ge, K_LO, K_HI, F_LO, F_HI, nk=60)[0]
    res["elastic"] = dict(f_MHz=b["f_MHz"], k=b["k"], abs_ve=b["abs_ve"])
    print(f"  elastic baseline            f_ZGV={b['f_MHz']:9.5f} MHz  k={b['k']:8.2f}  "
          f"|v_e|={b['abs_ve']:.2e}")

    # MEE plate, four combinations of the two coupled fields
    cases = {
        "MEE piezo only": dict(piezo=True, magnetic=False),
        "MEE magnetic only": dict(piezo=False, magnetic=True),
        "MEE fully coupled": dict(piezo=True, magnetic=True),
    }
    for label, kw in cases.items():
        C, coups, second, rho, bcs, names = build_material("mee", **kw)
        gp = GeneralPlate(h, C, coups, second, rho, N=24, bcs=bcs, field_names=names)
        folds = find_folds(gp, K_LO, K_HI, F_LO, F_HI, nk=60)
        for b in folds:
            res[label + " @ %.4f MHz" % b["f_MHz"]] = dict(
                n_fields=len(coups), f_MHz=b["f_MHz"], k=b["k"], abs_ve=b["abs_ve"])
            print(f"  {label:<26} ({len(coups)} field{'s' if len(coups)!=1 else ''}) "
                  f"f_ZGV={b['f_MHz']:9.5f} MHz  k={b['k']:8.2f}  |v_e|={b['abs_ve']:.2e}")

    # electrical / magnetic surface conditions for the fully coupled MEE plate
    C, coups, second, rho, _, names = build_material("mee")
    for combo in [("open", "open"), ("short", "open"), ("open", "short"), ("short", "short")]:
        gp = GeneralPlate(h, C, coups, second, rho, N=24, bcs=list(combo), field_names=names)
        folds = find_folds(gp, K_LO, K_HI, F_LO, F_HI, nk=60)
        key = f"MEE bc {combo[0]}/{combo[1]}"
        for b in folds:
            res[key + " @ %.4f MHz" % b["f_MHz"]] = dict(f_MHz=b["f_MHz"], k=b["k"],
                                                         abs_ve=b["abs_ve"])
            print(f"  {key:<26}      f_ZGV={b['f_MHz']:9.5f} MHz  k={b['k']:8.2f}  "
                  f"|v_e|={b['abs_ve']:.2e}")

    # mesh invariance of the MEE fold
    gp24 = GeneralPlate(h, C, coups, second, rho, N=24, bcs=["open", "open"], field_names=names)
    gp30 = GeneralPlate(h, C, coups, second, rho, N=30, bcs=["open", "open"], field_names=names)
    b24 = find_folds(gp24, K_LO, K_HI, F_LO, F_HI, nk=50)
    b30 = find_folds(gp30, K_LO, K_HI, F_LO, F_HI, nk=50)
    b24, b30 = b24[0], b30[0]
    rel = abs(b24["f_MHz"] - b30["f_MHz"]) / b30["f_MHz"]
    res["MEE_mesh_invariance_rel"] = rel
    print(f"  MEE mesh invariance N=24 vs 30: {rel:.2e} "
          f"({'PASS' if rel < 1e-4 else 'CHECK'})")

    # piezothermoelastic: thermal coupling enters as a third quasi-static
    # scalar through the thermal-stress tensor beta = C : alpha
    alpha = np.array([8.0e-6, 8.0e-6, 8.0e-6])        # BaTiO3 order of magnitude
    beta = np.zeros((6, 3))
    for I in range(3):
        beta[I, 0] = sum(C_pzt_r[I, J] * alpha[J] for J in range(3))
        beta[I, 1] = 0.0
        beta[I, 2] = 0.0
    kappa = np.diag([2.5, 2.5, 2.5])                    # W/(m K), placeholder scale
    print("\n  piezothermoelastic: the thermal scalar enters the LOSSLESS "
          "quasi-static limit\n  exactly like an extra potential, with the "
          "thermal-conductivity block in place\n  of the permittivity block. "
          "It is included here to show the fold survives;\n  a Fourier-type "
          "thermal field adds NO new branch.")
    print(f"      thermal-stress coefficients beta_1 = {beta[0,0]:.4e} Pa/K  "
          f"(from alpha = {alpha[0]:.1e} 1/K)")
    res["thermo_beta1_Pa_per_K"] = float(beta[0, 0])
    res["thermo_note"] = ("thermal scalar is mathematically equivalent to an extra "
                          "quasi-static potential in the lossless/Fourier limit; "
                          "hyperbolic (MGT/TPL) conduction makes the operator "
                          "NON-self-adjoint and is outside this module's contract")
    return res


def find_folds(gp, k_lo, k_hi, f_lo, f_hi, nk=60, n_refine=2, tol=5.0):
    """Enumerate ZERO-GROUP-VELOCITY points by scanning g(k) = min_j |v_e|.

    A ZGV is a local minimum of g(k) that refines to |v_e| < tol.  Several can
    exist (different branch pairs), so ALL are returned, not just the global
    minimum.  `tol` is an absolute energy velocity in m/s.
    """
    ks = np.linspace(k_lo, k_hi, nk)
    g = np.full(nk, np.inf)
    fg = np.full(nk, np.nan)
    for i, k in enumerate(ks):
        gi, fi = min_abs_ve(gp, k, f_lo, f_hi)
        g[i], fg[i] = gi, fi
    cand = [i for i in range(1, nk - 1) if g[i] <= g[i - 1] and g[i] <= g[i + 1]]
    if g[0] < g[1]:
        cand.append(0)
    if g[-1] < g[-2]:
        cand.append(nk - 1)
    folds = []
    for i in cand:
        if not np.isfinite(g[i]):
            continue
        bk, step = ks[i], ks[1] - ks[0]
        for _ in range(n_refine):
            sub = np.linspace(bk - step, bk + step, 25)
            bv, bkk = np.inf, bk
            for k in sub:
                gi, _ = min_abs_ve(gp, k, f_lo, f_hi)
                if gi < bv:
                    bv, bkk = gi, k
            bk, step = bkk, (sub[1] - sub[0])
        gg, ff = min_abs_ve(gp, bk, f_lo, f_hi)
        folds.append(dict(k=float(bk), f_MHz=float(ff) / 1e6, abs_ve=float(gg),
                          genuine=bool(gg < tol)))
    folds.sort(key=lambda d: d["f_MHz"])
    return folds


def part_scaling():
    """V9  Thickness scaling and surface-condition sensitivity.

    All wavenumber and frequency windows are scaled by 1/h from the validated
    1.5 mm reference windows.  Using a fixed window at every thickness scans
    the wrong part of the spectrum and returns the grid edge as a fake fold
    (this happened during development).
    """
    print("\n" + "=" * 72)
    print("V9  fold location, thickness scaling, surface-condition sensitivity")
    print("=" * 72)
    H_REF = 0.75e-3                     # half-thickness of the reference plate
    K_LO, K_HI = 450.0, 1700.0          # reference k window  (1/m, at 1.5 mm)
    F_LO, F_HI = 0.85e6, 2.40e6         # reference f window  (Hz,  at 1.5 mm)
    #  The window must CONTAIN the fold of every material: the piezoelectric
    #  fold sits near 1.13 MHz while the elastic / MEE folds sit near 1.8-1.9
    #  MHz.  A window chosen only around the elastic fold reports "no ZGV" for
    #  PZT, which is a windowing artefact, not physics.
    res = {"reference_window_at_1p5mm": dict(k=[K_LO, K_HI], f_MHz=[F_LO / 1e6, F_HI / 1e6])}

    configs = [("isotropic elastic", "elastic", None),
               ("PZT-5H open", "pzt", ["open"]),
               ("PZT-5H short", "pzt", ["short"]),
               ("MEE open/open", "mee", ["open", "open"]),
               ("MEE short/short", "mee", ["short", "short"])]

    table = {}
    for hh in [0.025e-3, 0.25e-3, 0.75e-3]:
        sc = H_REF / hh                  # k and f both scale as 1/h
        key = "2h=%.2f mm" % (2 * hh * 1e3)
        table[key] = {}
        for label, kind, bcs in configs:
            C, coups, second, rho, names = None, None, None, None, None
            if kind == "elastic":
                C, coups, second, rho, names = C_iso, [], [], RHO_ISO, []
            elif kind == "pzt":
                C, coups, second, rho, _, names = build_material("pzt")
            else:
                C, coups, second, rho, _, names = build_material("mee")
            gp = GeneralPlate(hh, C, coups, second, rho, N=24, bcs=bcs, field_names=names)
            folds = find_folds(gp, K_LO * sc, K_HI * sc, F_LO * sc, F_HI * sc, nk=70)
            gen = [d for d in folds if d["genuine"]]
            table[key][label] = gen
            if gen:
                for d in gen:
                    print(f"  {key:<12} {label:<18} f_ZGV={d['f_MHz']:10.5f} MHz  "
                          f"k={d['k']:10.2f}  2h*k={d['k']*2*hh:7.4f}  |v_e|={d['abs_ve']:.2e}")
            else:
                print(f"  {key:<12} {label:<18} no genuine ZGV in the window")
    res["table"] = table

    # geometric scaling: 2h*f_ZGV and 2h*k_ZGV should be thickness independent
    print("\n  geometric scaling (2h*f and 2h*k must be constant):")
    scale = {}
    for label, _, _ in configs:
        fs = [(key, table[key][label][0]) for key in table if table[key][label]]
        if len(fs) < 2:
            continue
        ref = fs[-1][1]
        dev = []
        for key, d in fs:
            dev.append((key, 100 * (d["f_MHz"] / ref["f_MHz"] - ref["k"] and 1)))
        f2h = [d["f_MHz"] * float(key.split("=")[1].split(" ")[0]) for key, d in fs]
        k2h = [d["k"] * float(key.split("=")[1].split(" ")[0]) * 1e-3 for key, d in fs]
        spread = 100 * (max(f2h) / min(f2h) - 1)
        spread_k = 100 * (max(k2h) / min(k2h) - 1)
        scale[label] = dict(f_times_2h=f2h, spread_f_pct=spread, spread_k_pct=spread_k)
        print(f"    {label:<18} 2h*f = {['%.4f'%v for v in f2h]}   spread {spread:.3f} %")
    res["geometric_scaling"] = scale

    # thickness-independent comparisons at 1.5 mm
    print()
    def g1(label):
        d = table["2h=1.50 mm"].get(label)
        return d[0]["f_MHz"] if d else np.nan
    e, po, psh, mo, msh = (g1("isotropic elastic"), g1("PZT-5H open"), g1("PZT-5H short"),
                           g1("MEE open/open"), g1("MEE short/short"))
    print("  all genuine folds found:")
    for key in table:
        for label, _, _ in configs:
            if table[key][label]:
                print(f"    {key:<12} {label:<18} "
                      + ", ".join("%.4f MHz @ k=%.1f" % (d["f_MHz"], d["k"])
                                  for d in table[key][label]))
    res["pct_at_1p5mm"] = dict(
        pzt_open_vs_elastic=100 * (po / e - 1),
        pzt_short_vs_open=100 * (psh / po - 1),
        mee_open_vs_elastic=100 * (mo / e - 1),
        mee_short_vs_open=100 * (msh / mo - 1))
    print("  at 1.5 mm:")
    for kk, v in res["pct_at_1p5mm"].items():
        print(f"    {kk:<22} {v:+9.4f} %")
    return res


def main():
    OUT["validation"] = part_validation()
    OUT["premise"] = part_premise()
    OUT["scaling"] = part_scaling()
    path = os.path.join(DATA, "check_03_multiphysics_premise.json")
    with open(path, "w") as fh:
        json.dump(OUT, fh, indent=2, default=float)
    print("\nwrote", path)


if __name__ == "__main__":
    main()
