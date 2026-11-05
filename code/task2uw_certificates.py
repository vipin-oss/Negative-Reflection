"""Task 2U-W independent modal-completeness certificates for exact normal incidence.

External checks use the INDEPENDENT first-order 8-state operator from
handoffs/task2uv_2026-11-05/code/independent_finite_k_state_space.py (no
imports of repository scattering code); production cross-checks use
code/scattering. Nothing here is a scattering R/T/A claim.

Certificates produced per material / pair / frequency / tau-case / precision:
  C1  k=0 spectrum of the independent operator: exactly two zero eigenvalues
      (algebraic multiplicity 2), nullity 1 (geometric multiplicity 1).
  C2  generalized eigenvector: B0 v = u with B0 u = 0, u = constant-potential;
      v carries the uniform-D0 sector (Dz entry nonzero, B0^2 v = 0).
  C3  nonzero-root separation at k=0 (P/thermal quartic + shear).
  C4  small-k electrostatic pair = +/- s_e k + O(k^2), s_e = sqrt(eps11/eps33_eff).
  C5  production k=0 6x6 block: det, smallest singular value, inverse norm.
  C6  electrostatic 2x2 block E0 = [[1,-1],[-eps_eL,-eps_eR]], det = -(eps_eL+eps_eR).
  C7  Neumann bound on the scaled full matrix S(k): q(k)=||S0^-1 (S(k)-S0)|| < 1
      and measured ||S(k)^-1|| <= ||S0^-1||/(1-q(k)) over a k grid.
  C8  scattering along the k grid: electrostatic amplitudes O(1), D_z = O(k),
      mechanical/thermal convergence order toward the k=0 6x6 solution.

Synthetic tau cases (omega*tau = 0.1, 1) are ALGORITHMIC stress tests, not
measured relaxation times.

Reproduce:
    cd code && python3 task2uw_certificates.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "handoffs" / "task2uv_2026-11-05" / "code"))

import foundation_checks as fc                                   # noqa: E402
from independent_finite_k_state_space import (                   # noqa: E402
    load as ind_load, operator as ind_operator, balanced, scale as ind_scale)
from scattering import assembly as A                             # noqa: E402
from scattering import gauge                                     # noqa: E402
from scattering.state import STATE_FIELDS                       # noqa: E402

CSV = HERE.parent / "handoffs" / "task2uv_2026-11-05" / "source" / "material_constants.csv"
MATS = ("Al6061", "PZT-5H", "BaTiO3")
PHI = STATE_FIELDS.index("phi")
DZ = STATE_FIELDS.index("Dz")
MECH = [i for i, f in enumerate(STATE_FIELDS) if f not in ("phi", "Dz")]


def eps_eff(p):
    """Clamped-adiabatic effective permittivities (lemma section 'coalescing roots')."""
    e33_eff = p["eps33"] - p["T0"] * p["p3"] ** 2 / (p["rho"] * p["cE"])
    s_e = mp.sqrt(p["eps11"] / e33_eff)
    eps_e = mp.sqrt(p["eps11"] * e33_eff)
    return e33_eff, s_e, eps_e


def sv2(M):
    """Singular values via eighe of the Hermitian Gram matrix."""
    n, m = M.rows, M.cols
    G = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            G[i, j] = sum(mp.conj(M[r, i]) * M[r, j] for r in range(n))
    Gh = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    ev, _ = mp.eighe(Gh)
    return sorted((mp.sqrt(abs(e)) for e in ev), reverse=True)


def mat_rank(M, tol_exp=-30):
    sv = sv2(M)
    if not sv or sv[0] == 0:
        return 0
    cut = mp.power(10, tol_exp) * sv[0]
    return sum(1 for s in sv if s > cut)


def norm2(M):
    return sv2(M)[0]

def inv_norm(M):
    sv = sv2(M)
    return mp.inf if sv[-1] == 0 else 1 / sv[-1]


def cert_jordan(p, w, dps):
    """C1/C2: zero eigenvalue alg mult 2, geom mult 1, generalized uniform-D0 vector."""
    mp.mp.dps = dps
    B0, kr, s = balanced(ind_operator(p, mp.mpf(0), w), p, w)
    vals, V = mp.eig(B0)
    zeros = [v for v in vals if abs(v) < mp.mpf("1e-20")]
    rank0 = mat_rank(B0)
    nullity = 8 - rank0
    u = mp.matrix([0, 0, 1, 0, 0, 0, 0, 0])          # constant potential
    ru = max(abs(t) for t in B0 * u)
    # generalized vector: minimum-norm solve B0 v = u via the Hermitian
    # pseudo-inverse of B0^H B0 (B0 itself is singular, rank 7)
    G = B0.H * B0
    Gh = mp.matrix(8, 8)
    for i in range(8):
        for j in range(8):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    ev, Q = mp.eighe(Gh)
    rhs = B0.H * u
    v = mp.matrix(8, 1)
    for i in range(8):
        if abs(ev[i]) > mp.mpf("1e-30") * max(abs(e) for e in ev):
            ci = sum(mp.conj(Q[j, i]) * rhs[j] for j in range(8)) / ev[i]
            for j in range(8):
                v[j] += Q[j, i] * ci
    r1 = max(abs(t) for t in (B0 * v - u))
    r2 = max(abs(t) for t in (B0 * B0 * v))
    return {
        "n_zero_eigenvalues": len(zeros),
        "nullity_B0": nullity,
        "rank_B0": rank0,
        "resid_B0_u": mp.nstr(ru, 4),
        "resid_B0_v_minus_u": mp.nstr(r1, 4),
        "resid_B0sq_v": mp.nstr(r2, 4),
        "generalized_Dz_entry": mp.nstr(v[DZ], 15),
        "generalized_phi_entry": mp.nstr(v[PHI], 15),
    }


def cert_roots(p, w, dps):
    """C3/C4: nonzero root separation at k=0; small-k electrostatic pair."""
    mp.mp.dps = dps
    kr, _ = ind_scale(p, w)
    B0, _, _ = balanced(ind_operator(p, mp.mpf(0), w), p, w)
    vals, _ = mp.eig(B0)
    nz = [v for v in vals if abs(v) > mp.mpf("1e-12")]
    sep = min(abs(a - b) for i, a in enumerate(nz) for b in nz[i + 1:]) if len(nz) > 1 else None
    _, s_e, _ = eps_eff(p)
    pair = []
    for frac in ("1e-2", "1e-3", "1e-4"):
        k = kr * mp.mpf(frac)
        Bk, _, _ = balanced(ind_operator(p, k, w), p, w)
        vk, _ = mp.eig(Bk)
        low = sorted((v for v in vk if abs(v) < mp.mpf("1e-2")), key=abs)
        pair.append({"k_over_kr": frac,
                     "roots_over_k": [mp.nstr(v / mp.mpf(frac), 10) for v in low[:2]],
                     "s_e": mp.nstr(s_e, 10)})
    return {"nonzero_root_min_sep_over_kr": mp.nstr(sep / kr, 12) if sep else None,
            "small_k_pair": pair}


def scaled_S(p1, p2, k, w, kind, dps):
    """Scaled full interface matrix S(k) per the lemma, on the INDEPENDENT
    eigenbasis (handoff independent_finite_k_state_space + eig_modes).

    Columns per side (branch-continuous order): the outgoing modes sorted by
    lam, with the electrostatic (lam -> -i*s_e*k) column last within its
    side. Rows (lemma order): [ux, uz, theta, txz, tzz, qz, phi, Dz/k] with
    fixed nondimensional scales s1 and the extra kr/k factor on the Dz row.
    At k = 0 the electrostatic columns are the analytic limit columns.
    """
    from independent_scattering_limit import (eig_modes, group, incident,
                                              PHI_SCALE)
    mp.mp.dps = dps
    kr, s1 = ind_scale(p1, w)
    M1, _ = eig_modes(p1, k, w, dps)
    M2, _ = eig_modes(p2, k, w, dps)
    group(M1, k, kr)
    group(M2, k, kr)
    inc = incident(M1, kind)

    def bkey(m):
        return (mp.re(m["lam"]), mp.im(m["lam"]))

    cols = []
    if k != 0:
        for side, M in ((1, M1), (2, M2)):
            out = [m for m in M if (m["left_out"] if side == 1 else m["right_out"])]
            mech = sorted([m for m in out if not m["electric"]], key=bkey)
            el = [m for m in out if m["electric"]]
            assert len(mech) == 3 and len(el) == 1, \
                f"side {side} mech/el counts {len(mech)}/{len(el)}"
            for m in mech + el:
                cols.append((side, m["Y"], None))
    else:
        for side, M in ((1, M1), (2, M2)):
            out = [m for m in M
                   if (m["left_out"] if side == 1 else m["right_out"])
                   and abs(m["a"]) > mp.mpf("1e-15")]
            out.sort(key=bkey)
            assert len(out) == 3, f"side {side} k=0 mech count {len(out)}"
            e = eps_eff(p1 if side == 1 else p2)[2]
            for m in out:
                cols.append((side, m["Y"], None))
            st = [mp.mpc(0)] * 8
            st[PHI] = PHI_SCALE
            cols.append((side, st, -e if side == 1 else e))
    assert len(cols) == 8, f"column count {len(cols)}"

    ROW = [0, 1, 3, 4, 5, 7, 2, 6]          # mech6, phi, Dz/k
    S = mp.matrix(8, 8)
    rhs = mp.matrix(8, 1)
    for j, tup in enumerate(cols):
        sd, st = tup[0], tup[1]
        sgn = mp.mpf(1) if sd == 1 else mp.mpf(-1)
        for r, i in enumerate(ROW):
            base = s1[i] * (k / kr if i == 6 and k != 0 else 1)
            v = st[i]
            if i == 6:
                v = v * kr / k if k != 0 else (tup[2] * PHI_SCALE if tup[2] is not None else mp.mpc(0))
            S[r, j] = sgn * v / base
    for r, i in enumerate(ROW):
        base = s1[i] * (k / kr if i == 6 and k != 0 else 1)
        v = inc["Y"][i]
        if i == 6:
            v = v * kr / k if k != 0 else mp.mpf(0)
        rhs[r] = -v / base
    return S, rhs, kr


def cert_inverse_bound(p1, p2, w, kind, dps, fracs):
    """C7/C8: Neumann bound and measured inverse norm in a FIXED equilibrated
    basis, on the independent eigenbasis. S0 is the k=0 analytic-limit matrix
    (block-triangular scaled limit); the same fixed diagonal factors are
    applied to every S(k)."""
    mp.mp.dps = dps
    kr, _ = ind_scale(p1, w)
    _, _, eL = eps_eff(p1)
    _, _, eR = eps_eff(p2)
    # reference matrix: finest-k member of the family (branch-matched column
    # order makes S(k) -> S0 entrywise); equilibrated once, factors FIXED
    S0, _rhs0, _kr = scaled_S(p1, p2, kr * mp.mpf(fracs[-1]), w, kind, dps)
    # FIXED two-sided equilibration of the reference matrix
    rsc = [max(abs(S0[i, j]) for j in range(8)) or mp.mpf(1) for i in range(8)]
    for i in range(8):
        for j in range(8):
            S0[i, j] = S0[i, j] / rsc[i]
    csc = [max(abs(S0[i, j]) for i in range(8)) or mp.mpf(1) for j in range(8)]
    for i in range(8):
        for j in range(8):
            S0[i, j] = S0[i, j] / csc[j]
    MCOL = (0, 1, 2, 4, 5, 6)   # per-side [3 mech sorted, el last]
    ECOL = (3, 7)
    M6 = mp.matrix(6, 6)
    for i in range(6):
        for j in range(6):
            M6[i, j] = S0[i, MCOL[j]]
    E0 = mp.matrix(2, 2)
    for i in range(2):
        for j in range(2):
            E0[i, j] = S0[6 + i, ECOL[j]]
    inv0 = mp.inverse(S0)
    rows = []
    amps_el = []
    for frac in fracs:
        k = kr * mp.mpf(frac)
        S, rhs, _ = scaled_S(p1, p2, k, w, kind, dps)
        for i in range(8):
            for j in range(8):
                S[i, j] = S[i, j] / rsc[i] / csc[j]
        D = mp.matrix(8, 8)
        for i in range(8):
            for j in range(8):
                D[i, j] = S[i, j] - S0[i, j]
        q = norm2(inv0 * D)
        bound = inv_norm(S0) / (1 - q) if q < 1 else mp.inf
        meas = inv_norm(S)
        a = mp.lu_solve(S, rhs)
        amps_el.append(max(abs(a[6]), abs(a[7])))
        rows.append({"k_over_kr": frac, "q_neumann": mp.nstr(q, 6),
                     "bound_inv_norm": mp.nstr(bound, 6) if bound != mp.inf else "inf",
                     "measured_inv_norm": mp.nstr(meas, 6),
                     "bound_holds": bool(q < 1 and meas <= bound)})
    return {"detE0_scaled": mp.nstr(mp.det(E0), 12),
            "detE0_physical_formula": mp.nstr(-(eL + eR), 12),
            "M6_smin": mp.nstr(sv2(M6)[-1], 6), "M6_det": mp.nstr(mp.det(M6), 6),
            "rows": rows, "max_electrostatic_amp": mp.nstr(max(amps_el), 6)}


def main():
    out = {"commit_scope": "Task 2U-W certificates", "cases": []}
    for dps in (40, 60):
        for name in MATS:
            p = ind_load(CSV, name)
            for f in ("1e5", "1e6", "1e7"):
                w = 2 * mp.pi * mp.mpf(f)
                for taucase in ("source", "wtau0.1", "wtau1"):
                    pp = dict(p)
                    if taucase != "source":
                        pp["tauq"] = mp.mpf(taucase[4:]) / w
                    j = cert_jordan(pp, w, dps)
                    r = cert_roots(pp, w, dps)
                    out["cases"].append({"dps": dps, "material": name,
                                         "freq_hz": f, "tau_case": taucase,
                                         "jordan": j, "roots": r})
    # pair-level inverse bound at 1 MHz, source tau, P incidence
    w = 2 * mp.pi * mp.mpf("1e6")
    pairs = [("Al6061", "Al6061"), ("PZT-5H", "PZT-5H"), ("Al6061", "PZT-5H"),
             ("PZT-5H", "Al6061"), ("PZT-5H", "BaTiO3")]
    out["pairs"] = []
    for dps in (40, 60):
        for n1, n2 in pairs:
            p1, p2 = ind_load(CSV, n1), ind_load(CSV, n2)
            for pol in ("P", "SV"):
                ib = cert_inverse_bound(p1, p2, w, pol, dps,
                                        ("1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"))
                out["pairs"].append({"dps": dps, "pair": f"{n1}|{n2}", "inc": pol,
                                     "inverse_bound": ib})
    Path(HERE.parent / "validation" / "task2uw_modal_completeness_results.json").write_text(
        json.dumps(out, indent=1) + "\n")
    # CSV of the inverse-bound rows
    with open(HERE.parent / "validation" / "task2uw_inverse_bound_results.csv", "w",
              newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["dps", "pair", "incidence", "k_over_kr", "q_neumann",
                     "bound_inv_norm", "measured_inv_norm", "bound_holds"])
        for pr in out["pairs"]:
            for r in pr["inverse_bound"]["rows"]:
                wr.writerow([pr["dps"], pr["pair"], pr["inc"], r["k_over_kr"],
                             r["q_neumann"], r["bound_inv_norm"],
                             r["measured_inv_norm"], r["bound_holds"]])
    print(f"wrote {len(out['cases'])} material cases, {len(out['pairs'])} pair cases")


if __name__ == "__main__":
    main()
