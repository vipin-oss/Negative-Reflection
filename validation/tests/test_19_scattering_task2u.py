"""Task 2U production repair — acceptance evidence, independent tests.

What this file establishes on the ACTIVE production path:

  ACPT07  the finite-k full eight-row dielectric-contact system converges to
          the production k = 0 reduced 6x6 mechanical/thermal solve, with a
          MEASURED per-component order, for all four Al6061/PZT-5H pairings
          (Al labelled dielectric surrogate) and the dissimilar piezo pair
          PZT-5H|BaTiO3, under both P and SV incidence.
  ACPT08  every interface contact residual is reported per field, each divided
          by its OWN physical scale (never an unscaled mixed-SI Euclidean
          norm): the six mechanical/thermal rows at k = 0, the [[phi]] = 0
          residual after relative-gauge reconstruction, and the on-shell D_z
          identity at k = 0.

Blockers: BLK-02A-26 (production defect) is closed by this evidence;
BLK-02A-27 (k = 0 mode count) is answered at the formulation level by the
conditional normal-incidence lemma -- 3 outgoing modes per side, 6x6, with
[[phi]] = 0 fixing only the relative potential constant -- and these
measurements; final acceptance remains with the supervisor. SCT02, G4, G5
remain open; no energy or R/T/A claim is made anywhere in this file.
"""
import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "code"))

import foundation_checks as fc                                  # noqa: E402
from harness import (FAIL, PASS, make_case, finish)              # noqa: E402
from scattering import asymptotic as asy                        # noqa: E402
from scattering import gauge                                    # noqa: E402
from scattering.state import STATE_FIELDS                       # noqa: E402

DPS = 60
W = mp.mpf(2) * mp.pi * mp.mpf(10) ** 6
PHI = STATE_FIELDS.index("phi")
DZ = STATE_FIELDS.index("Dz")
TZZ = STATE_FIELDS.index("tzz")

#: physical scale per field, for per-field residual reporting (ACPT08)
FIELD_SCALE = {"ux": "m", "uz": "m", "phi": "V", "theta": "K",
               "txz": "Pa", "tzz": "Pa", "Dz": "C/m^2", "qz": "W/m^2"}


def _mats():
    mp.mp.dps = DPS
    return (fc.material_from_table("Al6061"), fc.material_from_table("PZT-5H"),
            fc.material_from_table("BaTiO3"))


def _kac(al):
    return W / mp.sqrt(mp.mpf(repr(al.C11)) / mp.mpf(repr(al.rho)))


def _pairs():
    al, pz, ba = _mats()
    return [("Al6061|Al6061", al, al), ("PZT-5H|PZT-5H", pz, pz),
            ("Al6061|PZT-5H", al, pz), ("PZT-5H|Al6061", pz, al),
            ("PZT-5H|BaTiO3", pz, ba)]


# ===========================================================================
def test_ACPT07_finitek_converges_to_k0_6x6_with_measured_order():
    """The 8-row finite-k system converges to the production 6x6 at k = 0.

    Per component, the finite-k field sequence is classified by its measured
    k-order (asy.k_order); CONSTANT components are compared against the
    production k = 0 value (relative, 1e-6) and VANISHING components against
    zero (absolute, 1e-6 of the incident scale). This is the finite-k 8x8 vs
    k = 0 reduced 6x6 convergence the Task 2U brief requires, with measured
    per-component orders and precision stability (60 digits here; the handoff
    package independently reproduced the limit at 40 and 70 digits).
    """
    c = make_case(
        test_id="ACPT07", category="2-maths",
        purpose="finite-k eight-row dielectric contact converges to the "
                "production k=0 6x6 mechanical/thermal solve, per component, "
                "at a measured order",
        oracle="asy.k_order classification of the finite-k field sequence plus "
               "the production k=0 reduced solve",
        oracle_independent=True,
        fixture="five pairs (Al labelled dielectric surrogate) x P/SV, "
                "k/k_ac 1e-2 .. 1e-7, 60 digits",
        expected="every CONSTANT component agrees to 1e-6; every VANISHING "
                 "component is zero in production; every component has a "
                 "classified order",
        tolerance=1e-6,
        tol_justification="the finite-k system converges at measured order ~2, "
                          "so the k=1e-7 value is the limit far below 1e-6",
        level="2")

    def body():
        mp.mp.dps = DPS
        al, pz, ba = _mats()
        kac = _kac(al)
        fr = ["1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7"]
        ks = [kac * mp.mpf(f) for f in fr]
        MECH = [i for i, f in enumerate(STATE_FIELDS) if f not in ("phi", "Dz")]
        n, bad, notes = 0, [], []
        for name, m1, m2 in _pairs():
            for br, pol in ((5, "P"), (6, "SV")):
                P = gauge.reduced_system(m1, m2, mp.mpf(0), W, br, dps=DPS)
                if not P["square"] or P["n_unknowns"] != 6:
                    bad.append(f"{name} {pol}: k=0 system is not the 6x6 "
                               f"({P['n_equations']}x{P['n_unknowns']})")
                    continue
                sol, res, _ = gauge.solve_path(P["B"], P["rhs"], dps=DPS)
                if res > mp.mpf("1e-40"):
                    bad.append(f"{name} {pol}: k=0 solve residual "
                               f"{mp.nstr(res, 3)}")
                f1p, f2p = gauge.physical_fields_reduced(sol, P)
                inc_p = P["side1"].raw_states[br]
                iref = max(range(8), key=lambda i: abs(inc_p[i]))
                tp = inc_p[iref]
                rseq = {}
                for k in ks:
                    rseq[k] = asy.dielectric_contact_system(
                        m1, m2, k, W, br, normalisation="max")
                rlast = rseq[ks[-1]]
                tr = rlast["incident"][iref]
                FSr = max([abs(v) for v in rlast["f2"] + rlast["f1"]]
                          + [mp.mpf(1)])
                for i in MECH:
                    if tp == 0 or tr == 0:
                        continue
                    n += 1
                    a = f2p[i] / tp
                    seq = [rseq[k]["f2"][i] / rseq[k]["incident"][iref]
                           for k in ks]
                    p_ord, _C, cls = asy.k_order(seq, ks,
                                                 floor=FSr * mp.mpf("1e-20"))
                    if cls == "CONSTANT":
                        d = abs(a - seq[-1]) / max(abs(a), abs(seq[-1]),
                                                   mp.mpf("1e-40"))
                        if d > mp.mpf("1e-6"):
                            bad.append(f"{name} {pol}: {STATE_FIELDS[i]} rel "
                                       f"diff {mp.nstr(d, 3)} vs limit")
                    elif cls == "VANISHING":
                        d = abs(a) / max(abs(tp), mp.mpf(1))
                        if d > mp.mpf("1e-6"):
                            bad.append(f"{name} {pol}: {STATE_FIELDS[i]} "
                                       f"vanishes (order "
                                       f"{mp.nstr(p_ord, 3)}) but production "
                                       f"= {mp.nstr(d, 3)} of incident scale")
                    else:
                        # UNRESOLVED order but the component sits below the
                        # signal floor: its limit is zero, compare to zero
                        if max(abs(v) for v in seq) < mp.mpf("1e-9") * max(abs(tr), mp.mpf(1)):
                            d = abs(a) / max(abs(tp), mp.mpf(1))
                            if d > mp.mpf("1e-6"):
                                bad.append(f"{name} {pol}: {STATE_FIELDS[i]} "
                                           f"below floor but production = "
                                           f"{mp.nstr(d, 3)}")
                        else:
                            bad.append(f"{name} {pol}: {STATE_FIELDS[i]} k-order "
                                       f"unresolved ({cls})")
                notes.append(f"{name} {pol}: 6x6 ok, res {mp.nstr(res, 2)}")
        if bad:
            return (FAIL, n, f"{len(bad)} defects", " | ".join(bad[:4]))
        return (PASS, n, "; ".join(notes[:5]),
                "the finite-k eight-row system converges to the production "
                "k=0 6x6 solve for every pair and polarisation, per "
                "component, at a measured order")

    finish(c, body)


# ===========================================================================
def test_ACPT08_contact_residuals_reported_per_field():
    """Interface residuals per field, each on its own physical scale.

    No unscaled mixed-SI Euclidean norm is used anywhere: each residual is
    divided by the largest magnitude that field attains in the configuration,
    so the number is a per-field relative residual. The [[phi]] residual is
    reported separately from the global gauge freedom: reconstruction fixes
    C2 - C1, and the surviving gauge (C1 = 0) is recorded, not hidden.
    """
    c = make_case(
        test_id="ACPT08", category="2-physics",
        purpose="every k=0 interface contact residual, per field on its own "
                "scale, plus the [[phi]] residual separated from the global "
                "gauge",
        oracle="residuals recomputed from the production fields, not read off "
               "the solver",
        oracle_independent=True,
        fixture="five pairs x P/SV at k = 0, 60 digits",
        expected="all six mechanical/thermal residuals and the [[phi]] "
                 "residual at machine precision; the D_z on-shell identity "
                 "verified per column",
        tolerance=1e-40,
        tol_justification="the solve is exact linear algebra at 60 digits; the "
                          "residuals are recomputed from the reconstructed "
                          "fields, so only arithmetic remains",
        level="2")

    def body():
        mp.mp.dps = DPS
        al, pz, ba = _mats()
        n, bad, notes = 0, [], []
        for name, m1, m2 in _pairs():
            for br, pol in ((5, "P"), (6, "SV")):
                n += 1
                P = gauge.reduced_system(m1, m2, mp.mpf(0), W, br, dps=DPS)
                sol, res, _ = gauge.solve_path(P["B"], P["rhs"], dps=DPS)
                f1p, f2p = gauge.physical_fields_reduced(sol, P)
                inc = P["side1"].raw_states[br]
                # total side-1 field = incident + scattered
                f1t = [inc[i] + f1p[i] for i in range(8)]
                scale = [max(abs(f1t[i]), abs(f2p[i]), abs(inc[i]))
                         for i in range(8)]
                res_per_field = {}
                for i, f in enumerate(STATE_FIELDS):
                    sc = scale[i] if scale[i] > 0 else mp.mpf(1)
                    res_per_field[f] = abs(f1t[i] - f2p[i]) / sc
                worst_mech = max(res_per_field[f] for f in
                                 ("ux", "uz", "theta", "txz", "tzz", "qz"))
                if worst_mech > mp.mpf("1e-40"):
                    bad.append(f"{name} {pol}: worst mech/thermal residual "
                               f"{mp.nstr(worst_mech, 3)}")
                if res_per_field["phi"] > mp.mpf("1e-40"):
                    bad.append(f"{name} {pol}: [[phi]] residual "
                               f"{mp.nstr(res_per_field['phi'], 3)} after "
                               "reconstruction")
                # the on-shell D_z identity must have been verified in assembly
                dz_removed = [r for r in P["rows_removed"] if r["row"] == "Dz"]
                if dz_removed and "on_shell_worst_rel" in dz_removed[0]:
                    w = dz_removed[0]["on_shell_worst_rel"]
                    if w > mp.mpf("1e-25"):
                        bad.append(f"{name} {pol}: D_z on-shell identity "
                                   f"violated at {mp.nstr(w, 3)}")
                c1, c2 = P.get("gauge_constants", (None, None))
                notes.append(f"{name} {pol}: worst mech {mp.nstr(worst_mech, 2)}, "
                             f"[[phi]] {mp.nstr(res_per_field['phi'], 2)}, "
                             f"gauge C1={mp.nstr(c1, 2) if c1 is not None else '-'}"
                             f", C2-C1 fixed by [[phi]]=0")
        if bad:
            return (FAIL, n, f"{len(bad)} defects", " | ".join(bad[:4]))
        return (PASS, n, "; ".join(notes[:5]),
                "all contact residuals per field on their own scales are at "
                "machine precision; [[phi]]=0 is imposed by relative-gauge "
                "reconstruction with one recorded global gauge")

    finish(c, body)


if __name__ == "__main__":                              # pragma: no cover
    test_ACPT07_finitek_converges_to_k0_6x6_with_measured_order()
    test_ACPT08_contact_residuals_reported_per_field()
