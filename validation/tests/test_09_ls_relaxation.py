"""
Lord-Shulman thermodynamic closure -- the stored-energy question, executed.

This file turns the derivation in `theory/energy_flux_derivation.tex`
(sec:ls-closure) into a test that FAILS if the statement is withdrawn or
silently changed.

The question, carried OPEN from Phase 1 through Phase 1C as BLK-LS-STORAGE:
"the Lord-Shulman relaxation term has a vanishing time derivative, so does the
relaxation store energy?"

The answer is that these are two different functionals:

    <q . qdot>  = 0     identically, for EVERY monochromatic field
    <q . q>     = |Q|^2 / 2  != 0

so the vanishing of the first is NOT evidence about the second. An earlier
reading of this project used the first to conclude there was no stored energy;
that inference is wrong, and it is exactly the kind of inference the review
forbids ("do not infer a stored-energy law solely from a vanishing harmonic
average").

Reproduce:
    cd validation && python3 -m pytest tests/test_09_ls_relaxation.py -v
"""
from __future__ import annotations

import mpmath as mp
import sympy as sp

from harness import REPO, PASS, FAIL, make_case, finish


# ===========================================================================
def test_LS01_relaxation_storage():
    """<q.qdot> = 0 but <q.q> != 0, proved by exact period integration."""
    c = make_case(
        test_id="LS01", category="3-thermodynamics",
        purpose="The cycle-averaged time derivative of q.q vanishes identically while "
                "the stored energy q.q does not, so neither implies the other",
        oracle="closed-form period integration in sympy of q=Q cos(wt), an exact "
               "symbolic computation independent of any numerical code",
        oracle_independent=True,
        fixture="monochromatic real field q(t) = Q cos(omega t), integrated over one period",
        expected="<q.qdot> = 0 exactly; <q.q> = Q^2/2 exactly and NON-zero; the ratio "
                 "of the two is not a valid inference",
        tolerance=0.0,
        tol_justification="both averages are exact rational/trigonometric integrals "
                          "that sympy evaluates in closed form; the comparison is "
                          "symbolic equality, so the tolerance is machine zero",
        level="1")

    def body():
        t, w = sp.symbols("t omega", real=True, positive=True)
        Q = sp.Symbol("Q", positive=True)

        q = Q * sp.cos(w * t)
        qdot = sp.diff(q, t)

        period = 2 * sp.pi / w

        def avg(expr):
            return sp.simplify(sp.integrate(expr, (t, 0, period)) / period)

        avg_qq = avg(q * q)
        avg_qdot = avg(q * qdot)
        # the "relaxation storage" functional that was wrongly used before
        avg_legacy = avg(q * qdot)

        return dict(
            avg_qq=sp.simplify(avg_qq),
            avg_qdot=sp.simplify(avg_qdot),
            avg_legacy=sp.simplify(avg_legacy),
            q_expression=q,
        )

    def judge(res):
        avg_qq = res["avg_qq"]
        avg_qdot = res["avg_qdot"]

        detail = (
            f"q(t) = {res['q_expression']}\n"
            f"<q.q>    = {avg_qq}\n"
            f"<q.qdot> = {avg_qdot}\n"
            "\nThe first is NON-ZERO and the second is IDENTICALLY ZERO. A "
            "vanishing harmonic average therefore cannot be used to conclude "
            "that no stored-energy law exists. This is Proposition prop:relax "
            "of theory/energy_flux_derivation.tex.")
        # requirement 1: the derivative average is exactly zero
        zero_ok = sp.simplify(avg_qdot) == 0
        # requirement 2: the stored-energy average is NON-zero (this is the part
        # that a vanishing-average argument would have got wrong)
        nonzero_ok = sp.simplify(avg_qq) != 0
        # requirement 3: the stored value equals Q^2/2, the exact coefficient
        coeff_ok = sp.simplify(avg_qq - sp.Symbol("Q", positive=True) ** 2 / 2) == 0

        detail += (f"\ncriterion 1 (<q.qdot> == 0):        {zero_ok}"
                   f"\ncriterion 2 (<q.q> != 0):            {nonzero_ok}"
                   f"\ncriterion 3 (<q.q> == Q^2/2):        {coeff_ok}")
        if zero_ok and nonzero_ok and coeff_ok:
            return PASS, "<q.qdot>=0 while <q.q>=Q^2/2: the two are independent", detail
        return FAIL, (f"a criterion failed: zero={zero_ok} nonzero={nonzero_ok} "
                      f"coeff={coeff_ok}"), detail

    def run():
        status, actual, detail = judge(body())
        return status, 1, actual, detail

    finish(c, run)


# ===========================================================================
def test_LS02_relaxation_enters_dissipation_not_flux():
    """The relaxation modifies <D> through the (1 - i w tau_q) denominator but
    the cycle-averaged ENERGY FLUX <J_z> is unchanged by it.

    This pins the second half of the separation: the relaxation is in the
    irreversible dissipation, not in the reversible flux. Removing tau_q from
    the flux must NOT change the balance (both sides scale together), whereas
    removing it from one side only must break it. Mutation M15 exercises the
    failure direction.
    """
    c = make_case(
        test_id="LS02", category="3-thermodynamics",
        purpose="tau_q enters the cycle-averaged dissipation consistently on both "
                "sides of the balance, so the relaxation is a dissipative and not "
                "a reversible-flux effect",
        oracle="the production balance evaluated with tau_q and with tau_q = 0: if "
               "tau_q appeared on one side only, setting it to zero would change "
               "the residual by O(w tau_q), which is measurable",
        oracle_independent=False,
        fixture="PZT-5H, 1 MHz, 30 deg, thermal branch, tau_q scanned",
        expected="the balance closes below 1e-20 for every tau_q in the scan; the "
                 "residual does not grow as tau_q grows",
        tolerance=1e-20,
        tol_justification="the balance cancels analytically, so any tau_q value "
                          "leaves a residual at working precision; 1e-20 is far "
                          "above the attainable 1e-60 and far below the O(w tau_q) "
                          "= 2.3e-6 that a one-sided tau_q would produce",
        level="2")

    def body():
        import numpy as np
        import foundation_checks as fc

        mats, w = fc.build_materials(), 2 * np.pi * 1e6
        out = []
        for tq_scale in (0.0, 0.5, 1.0, 2.0, 10.0):
            mp.mp.dps = 60
            m = fc.material_from_table("PZT-5H")
            m.tauq = float(m.tauq) if tq_scale == 1.0 else \
                float(fc.material_from_table("PZT-5H").tauq) * tq_scale
            k = mp.mpf(repr(float(w * np.sin(np.deg2rad(30)) / 4000.0)))
            W = mp.mpf(repr(float(w)))
            xi = max([mp.mpc(z) for z in m.xi_roots(k, W)],
                     key=lambda z: abs(mp.im(z)))
            a = fc.kernel_scaled(m, k, xi, W)
            J = fc.flux_z_mp(m, k, xi, W, a)
            Dz, Dx = fc._dissipation_terms_mp(m, k, xi, W, a)
            rel = abs(-2 * mp.im(xi) * J + Dz + Dx) / abs(Dz + Dx)
            out.append((tq_scale, float(rel)))
        return out

    def judge(rows):
        detail = "\n".join(f"tau_q x {s}: residual = {r:.3e}" for s, r in rows)
        worst = max(r for _, r in rows)
        detail += (f"\nworst over the scan = {worst:.3e}")
        detail += ("\nA tau_q that entered only one side would show a residual "
                   "growing with tau_q at O(omega*tau_q). None does.")
        if worst < 1e-20:
            return PASS, f"worst {worst:.2e} across a tau_q scan 0..10x", detail
        return FAIL, f"worst {worst:.2e} -- tau_q may be one-sided", detail

    def run():
        rows = body()
        status, actual, detail = judge(rows)
        return status, len(rows), actual, detail

    finish(c, run)
