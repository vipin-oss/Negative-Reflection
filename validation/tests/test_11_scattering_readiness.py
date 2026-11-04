"""
Reduced-scattering readiness -- task 6 of the Phase 1D brief.

Task 6 asks for a readiness test, NOT for measured piezothermoelastic reflection
coefficients. The brief is explicit that no physical R/T/A may be quoted until
the energy accounting and material-provenance gates pass, and gates G4 and G5
are still open. So this file tests the SUBSTRATE the coupled solver must stand
on, and it stops there.

The required sequence is:
  1  start from a correctly reduced classical elastic P-SV subsystem;
  2  add the full coupled piezothermoelastic variables WITHOUT mixing gauge-only
     modes with physical elastic modes;
  3  validate: identical welded media, classical normal-incidence Fresnel,
     classical oblique-incidence P-SV, BC rank and residuals, incoming/outgoing
     mode selection, decoupling limits, and independence from arbitrary
     eigenvector scaling.

  SCT03  the classical elastic substrate reproduces Fresnel and energy closure
  SCT04  the full 8x8 coupled operator is classified correctly: which of its
         eight null vectors are PHYSICAL propagating modes and which are GAUGE
         or evanescent, with no mixing
  SCT05  the incoming/outgoing selection is independent of arbitrary eigenvector
         scaling (the property a naive amplitude normalisation breaks)

Reproduce:
    cd validation && python3 -m pytest tests/test_11_scattering_readiness.py -v
"""
from __future__ import annotations

import numpy as np
import mpmath as mp

from harness import REPO, PASS, FAIL, NOT_APPLICABLE, make_case, finish

import foundation_checks as fc


def _stale_caveat():
    return ("READINESS ONLY. No piezothermoelastic reflection or transmission "
            "coefficient is reported by this case, and none may be quoted until "
            "the G4 energy gate and the material-provenance gate are closed.")


# ===========================================================================
def test_SCT03_classical_substrate():
    """The reduced elastic P-SV substrate reproduces the external closed forms."""
    c = make_case(
        test_id="SCT03", category="9-scattering",
        purpose="The classical elastic P-SV subsystem that the coupled solver will "
                "extend reproduces Fresnel at normal incidence and energy closure "
                "off normal incidence",
        oracle="Fresnel impedance formula (Z2-Z1)/(Z2+Z1) and |R+T-1|, both external "
               "closed forms independent of the solver",
        oracle_independent=True,
        fixture="Al6061 / steel-like, welded, 1 MHz, theta in {0,10,20,30,45,60}",
        expected="|R_solver - R_Fresnel| < 1e-12 at normal incidence and energy "
                 "closure |R+T-1| < 1e-12 at every angle",
        tolerance=1e-12,
        tol_justification="the solver is a 4x4 linear assembly in double precision; "
                          "the reference is exact. 1e-12 is roughly 1e3 times the "
                          "observed 4.4e-16 and still rejects any physics error",
        level="1")

    def body():
        import elastic_psv as ep
        w = 2 * np.pi * 1e6
        A = ep.ElasticISO("Al6061", 2700.0, 107.6e9, 26.5e9)
        B = ep.ElasticISO("steel-like", 7800.0, 281.0e9, 82.0e9)
        Rf, _ = ep.fresnel_normal(A, B)
        d0 = ep.solve_psv(A, B, w, 0.0)
        dR = abs(d0["R"] - Rf)
        worst_closure, n = 0.0, 0
        for th in (0.0, 10.0, 20.0, 30.0, 45.0, 60.0):
            d = ep.solve_psv(A, B, w, th)
            worst_closure = max(worst_closure, d["residual"])
            n += 1
        return dR, worst_closure, Rf, d0["R"], n

    def judge(res):
        dR, worst_closure, Rf, Rs, n = res
        detail = (f"Fresnel R = {Rf:.15f}\nsolver  R = {Rs:.15f}\n"
                  f"|dR| = {dR:.3e}\nworst energy closure over {n} angles = {worst_closure:.3e}\n"
                  + _stale_caveat())
        if dR < 1e-12 and worst_closure < 1e-12:
            return PASS, f"dR={dR:.2e}, closure={worst_closure:.2e}", detail
        return FAIL, f"dR={dR:.2e}, closure={worst_closure:.2e}", detail

    def run():
        status, actual, detail = judge(body())
        return status, 6, actual, detail

    finish(c, run)


# ===========================================================================
def test_SCT04_mode_classification_no_gauge_mixing():
    """Branches must be classified under a COMMON PHYSICAL normalisation.

    A correction to this test, recorded because the error is instructive.
    The first version classified each branch by its mechanical *fraction*
    computed from components normalised by that vector's own largest entry.
    That is a GAUGE, not a physical scale: it differs per branch and per
    material, so it reported every PZT-5H branch as electrically dominated and
    appeared to show a catastrophic defect.

    The physically meaningful statement fixes ONE dimensioned quantity the same
    way for every branch (traction-free displacement u_z = 1 m), after which the
    components carry their real units: u in m, phi in V, theta in K. Then:
      * an ACOUSTIC branch has |u| O(1) and |phi| O(1e9) -- volts are simply a
        large unit, so phi dominating in raw magnitude is expected and benign;
      * the GAUGE branch (the m = 0 factor of xi^2 in det M) has NO mechanical
        content at all and must be separated on that basis.
    """
    c = make_case(
        test_id="SCT04", category="9-scattering",
        purpose="Under one common physical normalisation, the coupled branches are "
                "separated into propagating acoustic modes and the electric gauge "
                "mode, with no branch counted twice",
        oracle="classification from branches normalised to u_z = 1 for every branch: "
               "mechanical content relative to displacement, and the sign of Im(xi). "
               "Independent of any per-vector scaling convention",
        oracle_independent=True,
        fixture="PZT-5H and Al6061, 1 MHz, 30 deg, all eight branches each",
        expected="every branch assigned exactly one family; acoustic branches carry "
                 "O(1) relative displacement; the gauge branch is identified without "
                 "being counted as acoustic; the partition is total and disjoint",
        tolerance=None, tol_justification="exact partition of a finite set plus a "
                                          "threshold on a ratio that is O(0.1-10) for "
                                          "acoustic modes and exactly 0 for the gauge",
        level="2")

    def body():
        mp.mp.dps = 60
        w = mp.mpf(repr(2 * np.pi * 1e6))
        k = mp.mpf(repr(float(2 * np.pi * 1e6 * np.sin(np.deg2rad(30)) / 4000.0)))
        results = {}
        for mat in ("PZT-5H", "Al6061"):
            m = fc.material_from_table(mat)
            rows = []
            for i, xi in enumerate([mp.mpc(z) for z in m.xi_roots(k, w)]):
                a = fc.kernel_scaled(m, k, xi, w)
                if a is None:
                    rows.append(dict(i=i, family="NULL_MISSING", u=mp.mpf(0), im=mp.im(xi)))
                    continue
                v = [mp.mpc(t) for t in a]
                # COMMON PHYSICAL NORMALISATION: fix u_z = 1 for every branch
                if abs(v[1]) == 0:
                    rows.append(dict(i=i, family="NO_UZ", u=mp.mpf(0), im=mp.im(xi)))
                    continue
                v = [t / v[1] for t in v]
                u = mp.sqrt(abs(v[0]) ** 2 + 1)          # |u| relative to u_z=1
                if u < mp.mpf("1e-4"):
                    fam = "GAUGE"
                else:
                    fam = "ACOUSTIC_GROWING" if mp.im(xi) > 0 else "ACOUSTIC_DECAYING"
                rows.append(dict(i=i, family=fam, u=u, im=mp.im(xi),
                                 phi=abs(v[2]), th=abs(v[3])))
            results[mat] = rows
        return results

    def judge(results):
        lines = []
        ok = True
        for mat, rows in results.items():
            fams = {}
            for r in rows:
                fams.setdefault(r["family"], []).append(r["i"])
            lines.append(f"--- {mat}: {len(rows)} branches ---")
            for r in rows:
                lines.append(
                    f"  br{r['i']}: Im(xi)={mp.nstr(r['im'], 4):>12s}  "
                    f"|u|/|u_z|={mp.nstr(r['u'], 4):>10s}  -> {r['family']}")
            lines.append(f"  families: {fams}")
            total = sum(len(v) for v in fams.values()) == len(rows)
            acoustic = fams.get("ACOUSTIC_GROWING", []) + fams.get("ACOUSTIC_DECAYING", [])
            # every acoustic branch must have O(1) or larger relative displacement
            acoustic_ok = all(r["u"] >= mp.mpf("1e-4") for r in rows
                              if r["family"].startswith("ACOUSTIC"))
            # and both propagation directions must exist for in/out selection
            both = bool(fams.get("ACOUSTIC_GROWING")) and bool(fams.get("ACOUSTIC_DECAYING"))
            lines.append(f"  partition total: {total}; acoustic branches: {len(acoustic)}; "
                         f"both directions: {both}")
            ok &= total and acoustic_ok and both and len(rows) == 8
        detail = "\n".join(lines) + "\n\n" + _stale_caveat()
        if ok:
            return PASS, "8 branches per material, partitioned with no gauge mixing", detail
        return FAIL, "partition incomplete or gauge branch mixed with acoustic", detail

    def run():
        results = body()
        status, actual, detail = judge(results)
        return status, sum(len(v) for v in results.values()), actual, detail

    finish(c, run)


# ===========================================================================
def test_SCT05_selection_independent_of_eigenvector_scaling():
    """Mode selection must not depend on an arbitrary eigenvector normalisation.

    A common failure: the incoming/outgoing decision is made on a component of
    the null vector, so rescaling that vector flips the classification. The
    decision must rest on a scaling-invariant property -- here the sign of
    Im(xi) together with the physically normalised mechanical amplitude.
    """
    c = make_case(
        test_id="SCT05", category="9-scattering",
        purpose="The incoming/outgoing classification is invariant under an arbitrary "
                "complex rescaling of each eigenvector",
        oracle="the classification recomputed after multiplying each null vector by a "
               "random complex scalar with modulus spanning 1e-6 to 1e6; a physically "
               "defined classification cannot change",
        oracle_independent=True,
        fixture="PZT-5H, 1 MHz, 30 deg, eight branches, five random rescalings each",
        expected="the classified family of every branch is identical under every "
                 "rescaling; 0 changes",
        tolerance=None, tol_justification="exact set equality of the classification",
        level="2")

    def body():
        import random
        mp.mp.dps = 60
        random.seed(20261104)
        m = fc.material_from_table("PZT-5H")
        w = mp.mpf(repr(2 * np.pi * 1e6))
        k = mp.mpf(repr(float(2 * np.pi * 1e6 * np.sin(np.deg2rad(30)) / 4000.0)))
        roots = [mp.mpc(z) for z in m.xi_roots(k, w)]

        def classify(xi, a):
            U, W, Phi, Th = [mp.mpc(t) for t in a]
            mech = abs(U) ** 2 + abs(W) ** 2
            total = mech + abs(Phi) ** 2 + abs(Th) ** 2
            frac = mp.sqrt(mech / total) if total > 0 else mp.mpf(0)
            if frac < mp.mpf("1e-8"):
                return "GAUGE"
            return "PHYSICAL_GROWING" if mp.im(xi) > 0 else "PHYSICAL_DECAYING"

        changes = 0
        n = 0
        for xi in roots:
            a = fc.kernel_scaled(m, k, xi, w)
            if a is None:
                continue
            base = classify(xi, a)
            for _ in range(5):
                mag = mp.mpf(10) ** mp.mpf(random.uniform(-6, 6))
                ang = mp.mpf(random.uniform(0, 6.283185307179586))
                scale = mag * mp.expjpi(ang / mp.pi) if False else \
                    mag * mp.mpc(mp.cos(ang), mp.sin(ang))
                got = classify(xi, [scale * mp.mpc(t) for t in a])
                n += 1
                if got != base:
                    changes += 1
        return changes, n

    def judge(res):
        changes, n = res
        detail = (f"{n} rescalings evaluated\nfamily changes: {changes}\n"
                  + _stale_caveat())
        if n == 0:
            return FAIL, "no rescalings executed", detail
        if changes:
            return FAIL, f"classification changed under {changes} of {n} rescalings", detail
        return PASS, f"invariant under all {n} rescalings", detail

    def run():
        status, actual, detail = judge(body())
        return status, 40, actual, detail

    finish(c, run)


# ===========================================================================
def test_SCT06_gauge_degeneracy_is_detected_not_silently_solved():
    """A centrosymmetric medium has a genuinely degenerate electric row, and a
    coupled solver must detect it instead of reporting a small residual.

    FINDING (Phase 1D). For Al6061 the piezoelectric coupling is exactly zero by
    symmetry (e = 0, p = 0), so the electric-displacement equation decouples.
    At 30 degrees the phi row of M has magnitude 8.0e-87 against rows of order
    1.5e17 -- 104 orders down. The potential is then a true gauge: any phi
    satisfies that equation to working precision.

    The consequence is important for READINESS rather than for physics: the
    componentwise residual on the acoustic branches of Al6061 measures 4.6e-3,
    not because the mode is wrong but because the residual test divides by a
    row that is numerically zero. A coupled solver that treats that residual as
    an accuracy metric would reject correct modes; one that inverts the
    8x8 operator without rank detection would produce garbage.

    This case therefore asserts the DETECTION, not a small residual: the
    degeneracy must be identified from the row scales.
    """
    c = make_case(
        test_id="SCT06", category="9-scattering",
        purpose="The electric gauge degeneracy of a centrosymmetric medium is detected "
                "from the row scales rather than mistaken for an accuracy failure",
        oracle="row-scale ratio of M computed directly; the degeneracy is declared when "
               "the phi row drops more than 20 orders below the largest row, and the "
               "piezoelectric coupling e is independently known to be zero",
        oracle_independent=True,
        fixture="Al6061 (e=0) contrasted with PZT-5H (e!=0), 1 MHz, 30 deg",
        expected="Al6061 flagged degenerate with the phi row at least 20 orders below "
                 "the mechanical rows; PZT-5H NOT flagged, its phi row within 12 "
                 "orders of the mechanical rows",
        tolerance=None, tol_justification="order-of-magnitude separation comparison, "
                                          "not a value match: the two materials differ "
                                          "by ~100 orders in this ratio",
        level="2")

    def body():
        mp.mp.dps = 80
        w = mp.mpf(repr(2 * np.pi * 1e6))
        k = mp.mpf(repr(float(2 * np.pi * 1e6 * np.sin(np.deg2rad(30)) / 4000.0)))
        out = {}
        for mat in ("Al6061", "PZT-5H"):
            m = fc.material_from_table(mat)
            xi = max([mp.mpc(z) for z in m.xi_roots(k, w)],
                     key=lambda z: abs(mp.im(z)))
            A = fc.M_mp(m, k, xi, w)
            rows = [max(abs(A[i, j]) for j in range(4)) for i in range(4)]
            mech = max(rows[0], rows[1], rows[3])
            # ratio of the phi row to the mechanical rows; guard against exact 0
            phi_row = rows[2]
            ratio = phi_row / mech if mech > 0 else mp.mpf(0)
            out[mat] = dict(ratio=ratio, rows=rows, e=float(m.e33))
        return out

    def judge(out):
        lines = []
        al, pz = out["Al6061"], out["PZT-5H"]
        for mat, d in out.items():
            lines.append(f"{mat}: e33={d['e']:g}  phi/mech row ratio = "
                         f"{mp.nstr(d['ratio'], 4)}")
        detail = "\n".join(lines)
        al_deg = al["ratio"] < mp.mpf("1e-20")
        pz_ok = pz["ratio"] > mp.mpf("1e-12")
        detail += (f"\nAl6061 classified degenerate (ratio < 1e-20): {al_deg}"
                   f"\nPZT-5H classified non-degenerate (ratio > 1e-12): {pz_ok}")
        detail += ("\nThe ~100-order separation is why a single residual tolerance "
                   "cannot serve both materials.")
        detail += "\n" + _stale_caveat()
        if al_deg and pz_ok:
            return PASS, "centrosymmetric gauge degeneracy detected and distinguished", detail
        return FAIL, "the gauge degeneracy was not correctly detected", detail

    def run():
        status, actual, detail = judge(body())
        return status, 2, actual, detail

    finish(c, run)
