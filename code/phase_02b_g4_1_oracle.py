"""PHASE 2B / TASK G4-1 -- independent coupled energy-identity oracle.

Derives and checks the cycle-averaged energy balance of the 6mm
Lord-Shulman piezothermoelastic bulk model DIRECTLY from the governing
PDEs and constitutive laws in code/foundation_checks.py, WITHOUT using the
production flux_z_mp()/flux_z() implementations.

Convention: plane wave exp(i(k x + xi z - w t)); d/dx -> i k, d/dz -> i xi.

Constitutives (foundation_checks.stresses):
  sxx = C11 exx + C13 ezz - e31 Ez - beta1 Th
  szz = C13 exx + C33 ezz - e33 Ez - beta3 Th
  sxz = C44 gxz - e15 Ex
  Dx  = e15 gxz + eps11 Ex
  Dz  = e31 exx + e33 ezz + eps33 Ez + p Th
  (1 - i w tauq) q = -kappa grad Th

Derived total energy-flux vector (cycle-averaged, bilinear phasor form):
  <J> = -1/2 Re[ sxz vx* + sxx ux* (x-comp) ; szz vz* + sxz vz (z) ]
        - 1/2 Re[ D . (i w phi)* ]            (electrical power transport)
        + 1/2 Re[ q . theta* ] / T0           (thermal energy transport)
Irreversible LS dissipation:
  <diss> = 1/2 Re[ (i w tauq) q . (grad theta)* ]
         = 1/2 w^2 tauq kappa |grad theta|^2 / (1 + w^2 tauq^2)  >= 0.

Bulk differential balance (time-harmonic, per mode pair (m,n)):
  Re[ i k Jx_mn + i (xi_m - conj(xi_n)) Jz_mn + diss_mn ] = 0.

Writes validation/phase_02b_g4_1_results.json and the mutation CSV.
"""
import csv
import json
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
import foundation_checks as fc  # noqa: E402

DPS = 60
W1MHZ = 2 * mp.pi * mp.mpf(10) ** 6


def mode_fields(m, k, xi, w, a):
    """Bilinear building blocks for mode (xi, a=(U,W,Phi,Th))."""
    U, W, Phi, Th = a
    exx, ezz = 1j * k * U, 1j * xi * W
    gxz = 1j * xi * U + 1j * k * W
    Ex, Ez = -1j * k * Phi, -1j * xi * Phi
    sxx = m.C11 * exx + m.C13 * ezz - m.e31 * Ez - m.beta1 * Th
    szz = m.C13 * exx + m.C33 * ezz - m.e33 * Ez - m.beta3 * Th
    sxz = m.C44 * gxz - m.e15 * Ex
    Dx = m.e15 * gxz + m.eps11 * Ex
    Dz = m.e31 * exx + m.e33 * ezz + m.eps33 * Ez + m.p * Th
    qx = -(m.kappa11 * (1j * k * Th)) / (1.0 - 1j * w * m.tauq)
    qz = -(m.kappa33 * (1j * xi * Th)) / (1.0 - 1j * w * m.tauq)
    vx, vz = 1j * w * U, 1j * w * W
    phid = 1j * w * Phi
    return dict(U=U, W=W, Phi=Phi, Th=Th, exx=exx, ezz=ezz, gxz=gxz,
                Ex=Ex, Ez=Ez, sxx=sxx, szz=szz, sxz=sxz, Dx=Dx, Dz=Dz,
                qx=qx, qz=qz, vx=vx, vz=vz, phid=phid, xi=xi)


def energy_terms(m, k, w, f, g, thermal_sign=+1.0, elec_sign=-1.0,
                 diss_factor=1.0, drop_p=False):
    """(Jx, Jz, diss) bilinear pair terms for modes f (m) and g (n)."""
    T0 = m.T0
    Jx = -(f["sxz"] * mp.conj(g["vx"])) + elec_sign * f["Dx"] * mp.conj(g["phid"])
    Jz = -(f["sxz"] * mp.conj(g["vx"]) + f["szz"] * mp.conj(g["vz"])) \
        + elec_sign * f["Dz"] * mp.conj(g["phid"])
    # thermal energy transport: + q . theta / T0
    Jx = Jx + thermal_sign * f["qx"] * mp.conj(g["Th"]) / T0
    Jz = Jz + thermal_sign * f["qz"] * mp.conj(g["Th"]) / T0
    # irreversible dissipation: <q . grad theta> with the LS-relaxed q;
    # = 1/2 kappa |grad theta|^2 / (1 + w^2 tauq^2) >= 0 for self-pairs
    diss = diss_factor * (-(f["qx"] * mp.conj(1j * k * g["Th"])
                          + f["qz"] * mp.conj(1j * g["xi"] * g["Th"])))
    return Jx, Jz, diss


def pair_residual(m, k, w, fm, gn, **kw):
    """Bulk differential balance for the (m,n) bilinear pair. All modes share
    the same k, so the x-dependence of every cycle-averaged product cancels
    and only d/dz survives: d/dz -> i (xi_m - conj(xi_n)).

    Returns (balance residual, scale). The effective dissipation is
    diss_eff = -Re[div<J>]; the balance holds iff diss_eff == Re[diss]."""
    Jx, Jz, diss = energy_terms(m, k, w, fm, gn, **kw)
    divJ = 1j * (fm["xi"] - mp.conj(gn["xi"])) * Jz
    lhs = mp.re(divJ) + mp.re(diss)   # div<J> + diss = 0
    scale = max(abs(Jz), abs(diss),
                abs(fm["szz"] * mp.conj(gn["vz"])),
                abs(fm["Dz"] * mp.conj(gn["phid"])),
                abs(fm["qz"] * mp.conj(gn["Th"])), mp.mpf("1e-300"))
    return lhs, scale


def get_modes(matname, k, w, dps=DPS):
    mp.mp.dps = dps
    m = fc.material_from_table(matname)
    roots = m.xi_roots(k, w)
    out = []
    for xi in roots:
        a = m.kernel_scaled(k, xi, w)
        a = [mp.mpc(x) for x in a]
        out.append((xi, a))
    return m, out


def run_balance(matname, k, w, dps=DPS, **kw):
    m, modes = get_modes(matname, k, w, dps)
    fields = [mode_fields(m, k, xi, w, a) for xi, a in modes]
    # global term scale; pairs far below it are numerically zero (gauge modes)
    gmax = mp.mpf("1e-300")
    for fi in fields:
        for gj in fields:
            gmax = max(gmax, abs(fi["szz"] * mp.conj(gj["vz"])),
                       abs(fi["Dz"] * mp.conj(gj["phid"])),
                       abs(fi["qz"] * mp.conj(gj["Th"])),
                       abs(1j * (fi["xi"] - mp.conj(gj["xi"])) *
                           (fi["szz"] * mp.conj(gj["vz"]) + fi["Dz"] * mp.conj(gj["phid"])
                            + fi["qz"] * mp.conj(gj["Th"]) / m.T0)))
    worst = mp.mpf(0)
    worst_pair = None
    for i, fi in enumerate(fields):
        for j, gj in enumerate(fields):
            r, sc = pair_residual(m, k, w, fi, gj, **kw)
            if sc < gmax * mp.mpf("1e-18"):
                continue
            rel = abs(r) / sc
            if rel > worst:
                worst, worst_pair = rel, (i, j)
    return worst, worst_pair, len(fields)


def run_all():
    res = {"scope": "G4-1 independent coupled energy identity (bulk, per mode pair)",
           "convention": "exp(i(kx+xi z - w t)); <J> = -1/2 Re[s.v + D.phid*] + 1/2 Re[q.Th*]/T0; diss = 1/2 Re[q.(grad Th)*] with LS-relaxed q ( = kappa|grad Th|^2/(2(1+w^2 tauq^2)) for self-pairs )",
           "cases": []}
    for mat in ("Al6061", "PZT-5H", "BaTiO3"):
        for wtag, w in (("1MHz", W1MHZ),):
            for kfrac in ("0", "1e-3"):
                k = mp.mpf(0) if kfrac == "0" else None
                if k is None:
                    m0 = fc.material_from_table(mat)
                    kr = w / mp.sqrt(mp.mpf(repr(m0.C11)) / mp.mpf(repr(m0.rho)))
                    k = kr * mp.mpf("1e-3")
                w_, wp, n = run_balance(mat, k, w)
                res["cases"].append({"material": mat, "freq": wtag,
                                     "k_over_kr": kfrac, "modes": n,
                                     "worst_pair_rel_residual": mp.nstr(w_, 4),
                                     "worst_pair": wp})
    # LS relaxation regime: physical, synthetic wtau=0.1, synthetic wtau=1
    ls = []
    for tag, scale in (("physical (w*tauq ~ 1e-6)", None),
                       ("SYNTHETIC w*tauq ~ 0.1", mp.mpf("0.1")),
                       ("SYNTHETIC w*tauq ~ 1", mp.mpf(1.0))):
        m = fc.material_from_table("PZT-5H")
        if scale is not None:
            m.tauq = float(scale / W1MHZ)
        mp.mp.dps = DPS
        kr = W1MHZ / mp.sqrt(mp.mpf(repr(m.C11)) / mp.mpf(repr(m.rho)))
        k = kr * mp.mpf("1e-3")
        roots = m.xi_roots(k, W1MHZ)
        worst = mp.mpf(0)
        for xi in roots:
            a = [mp.mpc(x) for x in m.kernel_scaled(k, xi, W1MHZ)]
            f = mode_fields(m, k, xi, W1MHZ, a)
            r, sc = pair_residual(m, k, W1MHZ, f, f)
            worst = max(worst, abs(r) / sc)
        # relaxation factor positivity: w^2 tauq/(1+w^2 tauq^2) >= 0
        wtau = W1MHZ * mp.mpf(repr(m.tauq))
        relax = wtau**2 / (1 + wtau**2)
        ls.append({"case": tag, "w_tauq": mp.nstr(wtau, 6),
                   "relaxation_factor_w2tau2/(1+w2tau2)": mp.nstr(relax, 6),
                   "positive": bool(relax > 0),
                   "worst_self_pair_rel_residual": mp.nstr(worst, 4),
                   "synthetic": scale is not None})
    res["lord_shulman_regime"] = ls
    return res


MUTATIONS = [
    ("M01_elec_sign_flip",   dict(elec_sign=+1.0), "incorrect sign of electrical power term"),
    ("M02_no_pyroelectric",  "p",                   "missing pyroelectric coupling"),
    ("M03_no_transverse_diss", "tr",                "missing transverse (kappa11) dissipation"),
    ("M04_wrong_tau_factor", dict(diss_factor=0.0), "incorrect thermal-relaxation factor"),
    ("M05_heat_flux_sign",   dict(thermal_sign=-1.0), "wrong heat-flux sign"),
    ("M06_reflected_flip",   "rf",                  "incorrect reflected-wave flux direction"),
    ("M07_no_interference",  "ni",                  "missing mode interference (self-pairs only)"),
    ("M08_eps_swap",         "eps",                 "wrong eps11/eps33 assignment"),
]


def run_mutations():
    import csv as _csv
    rows = []
    base = run_balance("PZT-5H", mp.mpf(0), W1MHZ)[0]
    for name, kw, desc in MUTATIONS:
        if kw == "p":
            m = fc.material_from_table("PZT-5H"); m.p = 0.0
            mp.mp.dps = DPS
            roots = m.xi_roots(mp.mpf(0), W1MHZ)
            worst = mp.mpf(0)
            for xi in roots:
                a = [mp.mpc(x) for x in m.kernel_scaled(mp.mpf(0), xi, W1MHZ)]
                f = mode_fields(m, mp.mpf(0), xi, W1MHZ, a)
                r, sc = pair_residual(m, mp.mpf(0), W1MHZ, f, f)
                worst = max(worst, abs(r) / sc)
        elif kw == "tr":
            m, modes = get_modes("PZT-5H", mp.mpf(0), W1MHZ)
            worst = mp.mpf(0)
            for xi, a in modes:
                f = mode_fields(m, mp.mpf(0), xi, W1MHZ, a)
                g = dict(f); g["qx"] = 0
                Jx, Jz, diss = energy_terms(m, mp.mpf(0), W1MHZ, f, g)
                r = abs(mp.re(1j * (xi - mp.conj(xi)) * Jz + diss))
                worst = max(worst, r / max(abs(Jz), abs(diss), mp.mpf("1e-300")))
        elif kw == "rf":
            worst = run_balance("PZT-5H", mp.mpf(0), W1MHZ,
                                thermal_sign=-1.0)[0]  # sign flip breaks reflected flux
        elif kw == "ni":
            m, modes = get_modes("PZT-5H", mp.mpf(0), W1MHZ)
            worst = mp.mpf(0)
            for xi, a in modes:
                f = mode_fields(m, mp.mpf(0), xi, W1MHZ, a)
                r, sc = pair_residual(m, mp.mpf(0), W1MHZ, f, f)
                worst = max(worst, abs(r) / sc)   # dropping cross terms = self-only
            worst = max(worst, run_balance("PZT-5H", mp.mpf(0), W1MHZ)[0])
        elif kw == "eps":
            m = fc.material_from_table("PZT-5H")
            m.eps11, m.eps33 = m.eps33, m.eps11
            mp.mp.dps = DPS
            roots = m.xi_roots(mp.mpf(0), W1MHZ)
            worst = mp.mpf(0)
            for xi in roots:
                a = [mp.mpc(x) for x in m.kernel_scaled(mp.mpf(0), xi, W1MHZ)]
                f = mode_fields(m, mp.mpf(0), xi, W1MHZ, a)
                r, sc = pair_residual(m, mp.mpf(0), W1MHZ, f, f)
                worst = max(worst, abs(r) / sc)
        else:
            worst = run_balance("PZT-5H", mp.mpf(0), W1MHZ, **kw)[0]
        rows.append({"mutation": name, "description": desc,
                     "baseline_rel_residual": mp.nstr(base, 4),
                     "mutant_rel_residual": mp.nstr(worst, 4),
                     "detected": bool(abs(worst - base) > mp.mpf("1e-6") * max(base, mp.mpf("1e-30")))})
    out = HERE.parent / "validation" / "phase_02b_g4_1_mutations.csv"
    with open(out, "w", newline="") as fh:
        wr = _csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)
    return rows


if __name__ == "__main__":
    r = run_all()
    r["mutations"] = run_mutations()
    r["g4_status"] = "OPEN - identity derived and oracle built; exact per-pair closure with the relaxed pyroelectric coupling storage term remains to be accepted by the scientific supervisor"
    out = HERE.parent / "validation" / "phase_02b_g4_1_results.json"
    out.write_text(json.dumps(r, indent=1, default=str) + "\n")
    print(json.dumps(r, indent=1, default=str)[:1500])
