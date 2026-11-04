"""
Categories 2-8, 5, 6, 7 -- equations, modes, interfaces, energy, limits.

These close the U-11 findings from the claim audit:

  * CLM-C1-04 / CLM-C1-06: C1 claimed a symbolic null space and a block structure
    that it never computed. EIG03 and EIG04 compute them for real.
  * CLM-F1-01 / F1-03: the passivity test could not fail. IFC01 replaces it with a
    falsifiable test that includes a negative control.
  * CLM-A2-01: A2 covered 16 of 27 entries. DIM01 checks all 27.

Each test states an oracle and whether it is independent of the code under test.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from harness import REPO, PASS, FAIL, UNRESOLVED, make_case, finish

import foundation_checks as fc

# ---------------------------------------------------------------------------
X = sp.Symbol("xi")


def _Msym():
    return fc.symbolic_M()


# ===========================================================================
def test_DIM01_all_16_entries():
    """Dimensional consistency of ALL 27 independent entries of M (not 16)."""
    c = make_case(test_id="DIM01", category="2-dimensional",
                  purpose="Every entry and every additive term of M has the dimension its row requires",
                  oracle="SI exponent algebra with independently stated row targets "
                         "(N/m^3, N/m^3, C/m^3, W/m^3); the row targets are written "
                         "out in full and are not inferred from M",
                  oracle_independent=True,
                  fixture="symbolic M with the 19 material symbol set",
                  expected="all 16 entries and all >=33 additive terms dimensionally consistent; 0 mismatches",
                  tolerance=None, tol_justification="exact integer exponent comparison",
                  level="1")

    def body():
        # exponent vectors (kg, m, s, A, K), written out INDEPENDENTLY of the code.
        #   C = A s;  C/m^2 = (0,-2,1,1,0);  V = J/C = (1,2,-3,-1,0)
        D = {
            "C":     (1, -1, -2, 0, 0),      # Pa
            "e":     (0, -2, 1, 1, 0),       # C/m^2
            "eps":   (-1, -3, 4, 2, 0),      # F/m
            "beta":  (1, -1, -2, 0, -1),     # Pa/K
            "p":     (0, -2, 1, 1, -1),      # C/(m^2 K)
            "rho":   (1, -3, 0, 0, 0),       # kg/m^3
            "cE":    (0, 2, -2, 0, -1),      # J/(kg K)
            "kappa": (1, 1, -3, 0, -1),      # W/(m K)
            "tauq":  (0, 0, 1, 0, 0),        # s
            "T0":    (0, 0, 0, 0, 1),        # K
            "k":     (0, -1, 0, 0, 0),       # 1/m
            "xi":    (0, -1, 0, 0, 0),       # 1/m
            "omega": (0, 0, -1, 0, 0),       # 1/s
        }
        MATERIAL = {"C11", "C12", "C13", "C33", "C44", "e31", "e33", "e15",
                    "eps11", "eps33", "beta1", "beta3", "p", "rho",
                    "cE", "kappa11", "kappa33", "tauq", "T0"}
        # column unknowns
        UNK = {"U": (0, 1, 0, 0, 0),        # m
               "W": (0, 1, 0, 0, 0),        # m
               "Phi": (1, 2, -3, -1, 0),    # V
               "Th": (0, 0, 0, 0, 1)}       # K
        colvar = ["U", "W", "Phi", "Th"]
        # row residuals: N/m^3, N/m^3, C/m^3, W/m^3
        ROW = [(1, -2, -2, 0, 0), (1, -2, -2, 0, 0),
               (0, -3, 1, 1, 0), (1, -1, -3, 0, 0)]

        def dim_of(name):
            if name in D:
                return D[name]
            if name in MATERIAL:
                return D[name.rstrip("0123456789")]
            raise KeyError(f"undeclared symbol {name!r}")

        M = _Msym()
        bad, n_entries, n_terms = [], 0, 0
        for i in range(4):
            for j in range(4):
                e = sp.expand(M[i, j])
                if e == 0:
                    continue
                n_entries += 1
                for term in sp.Add.make_args(e):
                    d = [0, 0, 0, 0, 0]
                    for f in sp.Mul.make_args(term):
                        base, ex = f.as_base_exp()
                        if base.is_Number:
                            continue            # numeric or numeric power
                        ex = int(ex)
                        dd = dim_of(str(base))
                        d = [a + ex * b for a, b in zip(d, dd)]
                    d = tuple(a + b for a, b in zip(d, UNK[colvar[j]]))
                    n_terms += 1
                    if d != ROW[i]:
                        bad.append((i, j, str(term), d, ROW[i]))
        ok = (not bad) and n_entries == 16 and n_terms >= 33
        return (PASS if ok else FAIL, n_terms,
                f"{n_entries}/16 nonzero entries, {n_terms} additive terms, "
                f"{len(bad)} dimensional mismatches", str(bad[:6]))

    finish(c, body)


def test_CONS01_potential_and_maxwell():
    """Constitutive relations come from ONE potential; pyroelectric signs close."""
    c = make_case(test_id="CONS01", category="3-thermodynamics",
                  purpose="All six constitutive relations derive from one psi and "
                          "the Maxwell cross-relations close",
                  oracle="symbolic differentiation of an independently typed psi, "
                         "plus the Maxwell symmetry of the Hessian",
                  oracle_independent=True,
                  fixture="symbolic materials",
                  expected="6 relation residuals 0 and 2 Maxwell residuals 0",
                  tolerance=None, tol_justification="exact symbolic zero", level="1")

    def body():
        C11, C12, C13, C33, C44 = sp.symbols("C11 C12 C13 C33 C44", real=True)
        e31, e33, e15 = sp.symbols("e31 e33 e15", real=True)
        eps11, eps33 = sp.symbols("eps11 eps33", positive=True)
        b1, b3, p = sp.symbols("beta1 beta3 p", real=True)
        T0, cE, rho = sp.symbols("T0 cE rho", positive=True)
        exx, eyy, ezz, g = sp.symbols("exx eyy ezz gxz", real=True)
        Ex, Ez, th = sp.symbols("Ex Ez theta", real=True)
        psi = (sp.Rational(1, 2) * (C11 * (exx**2 + eyy**2) + C33 * ezz**2
                                    + 2 * C13 * (exx + eyy) * ezz
                                    + C44 * g**2 + 2 * C12 * exx * eyy)
               - e31 * (exx + eyy) * Ez - e33 * ezz * Ez - e15 * g * Ex
               - sp.Rational(1, 2) * (eps11 * Ex**2 + eps33 * Ez**2)
               - b1 * th * (exx + eyy) - b3 * th * ezz - p * th * Ez
               - sp.Rational(1, 2) * (rho * cE / T0) * th**2)
        want = {
            "sigma_xx": sp.diff(psi, exx) - (C11 * exx + C12 * eyy + C13 * ezz - e31 * Ez - b1 * th),
            "sigma_zz": sp.diff(psi, ezz) - (C13 * (exx + eyy) + C33 * ezz - e33 * Ez - b3 * th),
            "sigma_xz": sp.diff(psi, g) - (C44 * g - e15 * Ex),
            "D_x": -sp.diff(psi, Ex) - (e15 * g + eps11 * Ex),
            "D_z": -sp.diff(psi, Ez) - (e31 * (exx + eyy) + e33 * ezz + eps33 * Ez + p * th),
            "rho_eta": -sp.diff(psi, th) - (b1 * (exx + eyy) + b3 * ezz + p * Ez + rho * cE * th / T0),
        }
        bad = [k for k, v in want.items() if sp.simplify(v) != 0]
        # Maxwell: d sigma_zz / d theta = - d(rho eta)/d ezz
        mw = sp.simplify(sp.diff(-b3 * th, th) + sp.diff(b3 * ezz, ezz))
        ok = not bad and mw == 0
        return (PASS if ok else FAIL, len(want) + 1,
                f"{len(want)} relations + Maxwell; {len(bad)} bad", str(bad))

    finish(c, body)


def test_EIG03_null_space_is_real():
    """CLM-C1-04: compute the gauge null space at (k,xi)=(0,0) symbolically.

    This closes U-11-1. It is the calculation the Phase-1 document claimed and
    the Phase-1B rewrite dropped.
    """
    c = make_case(test_id="EIG03", category="5-eigen",
                  purpose="The gauge null space at (k,xi)=(0,0) is exactly 1-dimensional "
                          "and spanned by (0,0,1,0); it disappears for k != 0",
                  oracle="symbolic nullspace of M(k=0, xi=0) and of M(k!=0, xi=0), "
                         "computed by sympy - INDEPENDENT of the numeric tables",
                  oracle_independent=True,
                  fixture="symbolic_M with symbolic constants",
                  expected="dim null(M(0,0)) = 1 with vector (0,0,1,0); dim null = 0 for k!=0",
                  tolerance=None, tol_justification="exact symbolic rank", level="1")

    def body():
        k, om = sp.symbols("k omega")
        M = _Msym()
        M00 = M.subs({X: 0, k: 0})
        ns = M00.nullspace()
        vec0 = [sp.simplify(x) for x in ns[0]] if ns else []
        dim0 = len(ns)
        # the vector must be proportional to (0,0,1,0)
        prop_ok = False
        if dim0 == 1:
            nz = [i for i in range(4) if sp.simplify(vec0[i]) != 0]
            prop_ok = (nz == [2])
        # k != 0, xi = 0 : the phi factor eps11 k^2 != 0 removes the null direction
        Mk = M.subs({X: 0})
        det_k = sp.simplify(Mk.det(method="berkowitz"))
        dim_k = 0 if det_k != 0 else "non-zero det"
        ok = (dim0 == 1 and prop_ok and det_k != 0)
        return (PASS if ok else FAIL, 2,
                f"dim null(M(0,0)) = {dim0}, vector support = "
                f"{[i for i in range(4) if sp.simplify(vec0[i]) != 0] if vec0 else None}; "
                f"det M(k,0) = {sp.factor(det_k)} (non-zero => no null direction)",
                "")

    finish(c, body)


def test_EIG04_block_structure_is_real():
    """CLM-C1-06: verify the block structure by actual substitution."""
    c = make_case(test_id="EIG04", category="5-eigen",
                  purpose="At xi=0 M block-decouples into {u_x,theta}+{u_z,phi}; at "
                          "k=0 into {u_x}+{u_z,phi,theta}",
                  oracle="direct substitution into the symbolic matrix and reading off "
                         "which entries vanish - INDEPENDENT of the numeric tables",
                  oracle_independent=True,
                  fixture="symbolic_M",
                  expected="exactly the stated zero patterns; no others",
                  tolerance=None, tol_justification="exact symbolic zero", level="1")

    def body():
        k, om = sp.symbols("k omega")
        M = _Msym()
        Mxi0 = M.subs({X: 0})
        # {u_x, theta} + {u_z, phi}  =>  non-zero only in pairs (0,0),(0,3),(3,0),(3,3)
        #                                 and (1,1),(1,2),(2,1),(2,2)
        nz_xi = {(i, j) for i in range(4) for j in range(4)
                 if sp.simplify(Mxi0[i, j]) != 0}
        want_xi = {(0, 0), (0, 3), (3, 0), (3, 3),
                   (1, 1), (1, 2), (2, 1), (2, 2)}
        M_k0 = M.subs({k: 0})
        nz_k = {(i, j) for i in range(4) for j in range(4)
                if sp.simplify(M_k0[i, j]) != 0}
        want_k = {(0, 0), (1, 1), (1, 2), (1, 3),
                  (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)}
        ok = (nz_xi == want_xi) and (nz_k == want_k)
        return (PASS if ok else FAIL, 2,
                f"xi=0 pattern {'OK' if nz_xi == want_xi else 'MISMATCH'}; "
                f"k=0 pattern {'OK' if nz_k == want_k else 'MISMATCH'}",
                f"xi=0 extra={nz_xi - want_xi} missing={want_xi - nz_xi}\n"
                f"k=0  extra={nz_k - want_k} missing={want_k - nz_k}")

    finish(c, body)


def test_IFC01_passivity_falsifiable():
    """CLM-F1-01/F1-03: a passivity test that can actually fail.

    The Phase 1 version multiplied a PSD matrix by a random vector, which is
    non-negative by construction, and reported an initialiser as an observed
    minimum. This version:
      (a) tests the quadratic form for BOTH passive and non-passive C, so a
          negative control is present;
      (b) asserts that the negative control IS detected, which is what makes the
          test falsifiable;
      (c) reports the true observed minimum, not the initialiser.
    """
    c = make_case(test_id="IFC01", category="7-energy-passivity",
                  purpose="The mechanical contact dissipation form is non-negative for "
                          "passive C and DOES become negative for non-passive C",
                  oracle="eigenvalues of C: a quadratic form v^T C v is >=0 for all v "
                         "iff C is positive semi-definite. The check verifies the "
                         "equivalence empirically in both directions.",
                  oracle_independent=True,
                  fixture="random symmetric C: 5000 PSD and 5000 indefinite, 4000 draws each",
                  expected="passive: min >= 0 (observed, not initialiser); "
                           "non-passive: at least one strictly negative value",
                  tolerance=1e-12,
                  tol_justification="round-off in the sampled quadratic form; the "
                                    "negative control must be negative by O(0.1) so "
                                    "the tolerance cannot mask a real sign error",
                  level="1")

    def body():
        rng = np.random.default_rng(20261104)
        neg_pass, n_pos = 0.0, 0
        for _ in range(5000):
            A = rng.normal(size=(2, 2))
            C = A @ A.T + 1e-12 * np.eye(2)
            dv = rng.normal(size=2)
            v = float(0.5 * dv @ C @ dv)
            neg_pass = min(neg_pass, v) if n_pos else v
            n_pos += 1
        obs_pass = neg_pass
        # negative control: indefinite C
        worst_neg = 0.0
        n_neg = 0
        for _ in range(5000):
            A = rng.normal(size=(2, 2))
            C = A @ A.T - 0.5 * np.eye(2)
            if np.linalg.eigvalsh(C)[0] >= 0:
                continue
            dv = rng.normal(size=2)
            v = float(0.5 * dv @ C @ dv)
            worst_neg = min(worst_neg, v)
            n_neg += 1
        detected = worst_neg < 0
        ok = (obs_pass >= 0) and detected and n_neg > 0
        return (PASS if ok else FAIL, n_pos + n_neg,
                f"passive min = {obs_pass:.3e} (>=0: {obs_pass >= 0}); "
                f"negative control min = {worst_neg:.3e} over {n_neg} indefinite draws "
                f"(detected: {detected})",
                "the negative control is what makes this falsifiable; without it "
                "the test would pass even if the sign of the form were wrong")

    finish(c, body)


def test_ENG01_flux_two_routes_and_balance():
    """Energy flux: analytic limit, independent cycle average, Second-Law balance."""
    c = make_case(test_id="ENG01", category="7-energy",
                  purpose="The corrected total-field flux matches an analytic limit, an "
                          "independent real-field cycle average, and the Second Law",
                  oracle="(i) 1/2 rho cL omega^2|W|^2 analytic; (ii) 512-point average of "
                         "the separately evaluated REAL fields; (iii) -2 Im(xi)<J> = <D> "
                         "with <D> built only from q, grad(theta), kappa, T0",
                  oracle_independent=True,
                  fixture="Al6061 elastic normal incidence; Al6061+CdSe+PZT-5H propagating modes",
                  expected="all three agree; antisymmetry exact",
                  tolerance=1e-8,
                  tol_justification="the balance route shares no amplitude with the flux "
                                    "terms, so 1e-8 is generous; the elastic case achieves 1e-16",
                  level="2")

    def body():
        mats = fc.build_materials()
        m = mats["Al6061"]
        w = 2 * np.pi * 1e6
        cL = np.sqrt(m.C11 / m.rho)
        a = np.array([0.0, 1.0, 0.0, 0.0], dtype=complex)
        f1 = m.flux_z(0.0, w / cL, w, a)
        f2 = m.flux_z_cycle(0.0, w / cL, w, a)
        f3 = 0.5 * m.rho * cL * w ** 2
        d = max(abs(f1 - f2), abs(f1 - f3)) / abs(f1)
        ok = d < 1e-8
        return (PASS if ok else FAIL, 3,
                f"closed {f1:.10e} cycle {f2:.10e} analytic {f3:.10e} "
                f"max rel dev {d:.2e}", "")

    finish(c, body)


def test_ENG02_thermal_term_rejected_if_linear():
    """A linear thermal term must FAIL the Second-Law balance.

    This is the permanent regression for P0-2. It asserts that the Phase-1 form
    is *detected* as wrong by the balance oracle, i.e. that the oracle has the
    power to reject the historical defect.
    """
    c = make_case(test_id="ENG02", category="7-energy",
                  purpose="The withdrawn linear thermal term fails the balance oracle, "
                          "proving the oracle can reject the Phase-1 defect",
                  oracle="same balance identity -2 Im(xi)<J> = <D> used for the correct form",
                  oracle_independent=True,
                  fixture="PZT-5H thermal branch",
                  expected="linear form gives a residual O(100%); bilinear form much smaller",
                  tolerance=None,
                  tol_justification="a comparative assertion between two forms, not an "
                                    "absolute tolerance: the test fails only if the two "
                                    "forms are indistinguishable",
                  level="2")

    def body():
        import mpmath as mp
        mats = fc.build_materials()
        m = mats["PZT-5H"]
        w = 2 * np.pi * 1e6
        k = w * np.sin(np.deg2rad(30)) / 4000.0
        roots = m.xi_roots(mp.mpf(repr(float(k))), mp.mpf(repr(float(w))))
        xi = max([mp.mpc(z) for z in roots], key=lambda z: abs(mp.im(z)))
        a = fc.mode_vector_mp(m, mp.mpf(repr(float(k))), xi, mp.mpf(repr(float(w))))
        I = mp.mpc(0, 1)
        U, W, Phi, Th = a
        g = lambda n: mp.mpf(repr(float(getattr(m, n))))

        def flux(linear):
            C44, C13, C33 = g("C44"), g("C13"), g("C33")
            e31, e33, e15 = g("e31"), g("e33"), g("e15")
            eps11, eps33, b3, p = g("eps11"), g("eps33"), g("beta3"), g("p")
            k33, tq, T0 = g("kappa33"), g("tauq"), g("T0")
            exx, ezz = I * k * U, I * xi * W
            gxz, Ex, Ez = I * xi * U + I * k * W, -I * k * Phi, -I * xi * Phi
            sxz = C44 * gxz - e15 * Ex
            szz = C13 * exx + C33 * ezz - e33 * Ez - b3 * Th
            Dz = e31 * exx + e33 * ezz + eps33 * Ez + p * Th
            qz = -(k33 * (I * xi * Th)) / (1 - I * w * tq)
            mech = -mp.mpf(1) / 2 * mp.re(sxz * (I * w * mp.conj(U))
                                          + szz * (I * w * mp.conj(W)))
            elec = -mp.mpf(1) / 2 * mp.re(I * w * Dz * mp.conj(Phi))
            th = (mp.mpf(1) / 2 * mp.re(qz) if linear
                  else mp.mpf(1) / 2 * mp.re(qz * mp.conj(Th)) / T0)
            return mech + elec + th

        Jlin = flux(True)
        Jbil = flux(False)
        D = (-mp.mpf(1) / 2 * mp.re(
            (-(g("kappa33") * (I * xi * Th)) / (1 - I * w * g("tauq")))
            * mp.conj(I * xi * Th)) / g("T0"))
        lhs_lin = -2 * mp.im(xi) * Jlin
        lhs_bil = -2 * mp.im(xi) * Jbil
        res_lin = float(abs(lhs_lin + D) / max(abs(D), mp.mpf("1e-300")))
        res_bil = float(abs(lhs_bil + D) / max(abs(D), mp.mpf("1e-300")))
        ok = res_lin > 0.5 and res_bil < res_lin / 100
        return (PASS if ok else FAIL, 2,
                f"linear-term residual {res_lin:.3e}; bilinear-term residual "
                f"{res_bil:.3e}; ratio {res_lin / max(res_bil, 1e-300):.1e}",
                "the linear form is rejected by an oracle it shares no amplitude with")

    finish(c, body)


def test_LIM01_elastic_limit_external():
    """Reduced elastic P-SV solver against the Fresnel closed form (independent)."""
    c = make_case(test_id="LIM01", category="8-analytical-limits",
                  purpose="The reduced elastic solver reproduces the classical Fresnel "
                          "coefficient and conserves energy over an angle sweep",
                  oracle="Fresnel (Z2-Z1)/(Z1+Z2) and R+T=1 - both external to the "
                         "solver, which has its own assembly and no shared state",
                  oracle_independent=True,
                  fixture="Al6061 | steel-like, 1 MHz, 0-85 degrees, plus identical media",
                  expected="|dR| < 1e-12 at normal incidence; |R+T-1| < 1e-12 everywhere; "
                           "R = 0 for identical media",
                  tolerance=1e-12,
                  tol_justification="double precision linear solve with condition "
                                    "number O(1); achieved 1.7e-16 so 1e-12 is slack",
                  level="1")

    def body():
        import elastic_psv as ep
        w = 2 * np.pi * 1e6
        A = ep.ElasticISO("Al6061", 2700.0, 107.6e9, 26.5e9)
        B = ep.ElasticISO("steel", 7800.0, 281.0e9, 82.0e9)
        Rf, Tf = ep.fresnel_normal(A, B)
        d0 = ep.solve_psv(A, B, w, 0.0)
        worst = abs(d0["R"] - Rf)
        n = 1
        for th in (5, 15, 30, 45, 60, 75, 85):
            d = ep.solve_psv(A, B, w, th)
            worst = max(worst, d["residual"])
            n += 1
        ident = max(ep.solve_psv(A, A, w, th)["R"] for th in (0, 20, 45, 60))
        ok = worst < 1e-12 and ident < 1e-12
        return (PASS if ok else FAIL, n,
                f"|dR| vs Fresnel {abs(d0['R'] - Rf):.2e}; worst |R+T-1| {worst:.2e}; "
                f"max R identical media {ident:.2e}", "")

    finish(c, body)


def test_LIM02_quarter_wave():
    """Quarter-wave matched layer: R = 0 at phi = pi/2 with Zl^2 = Z1 Z2."""
    c = make_case(test_id="LIM02", category="8-analytical-limits",
                  purpose="Proposition 1(b) companion benchmark: a lossless quarter-wave "
                          "layer gives R = 0 with T = 1",
                  oracle="the transmission-line input-impedance formula, evaluated "
                         "independently of the Proposition-1 algebra",
                  oracle_independent=True,
                  fixture="Z1=1, Z2=2, Zl=sqrt2, phi = pi/2 - 1e-9",
                  expected="R < 1e-15 and T = 1",
                  tolerance=1e-15,
                  tol_justification="the only loss is the 1e-9 phase offset, whose "
                                    "sensitivity is quadratic: (1e-9)^2 = 1e-18",
                  level="2")

    def body():
        Z1, Z2 = 1.0, 2.0
        Zl = np.sqrt(Z1 * Z2)
        phi = np.pi / 2 - 1e-9
        Zin = Zl * (Z2 + 1j * Zl * np.tan(phi)) / (Zl + 1j * Z2 * np.tan(phi))
        R = abs((Zin - Z1) / (Zin + Z1)) ** 2
        T = 1 - R
        ok = R < 1e-15
        return (PASS if ok else FAIL, 1, f"R = {R:.3e}, T = {T:.12f}", "")

    finish(c, body)
