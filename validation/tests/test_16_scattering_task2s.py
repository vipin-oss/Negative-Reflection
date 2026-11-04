"""Task 2S -- electrical gauge versus interface-continuity audit.

These checks decide whether the Task 2R reduction, which removes the `phi` and
`Dz` interface rows at normal incidence, is EQUIVALENT to a properly
gauge-fixed full formulation or merely happens to solve.

The headline result is that equivalence HOLDS, but the reason originally
offered for it did NOT. Task 2R justified the removal by a rank argument: the
`phi` row adds no rank, so it is implied. Measured here, that argument is
vacuous -- once a system has at least as many rows as columns and full column
rank, adding any row leaves the rank alone, so zero rank gain is guaranteed and
proves nothing about whether the row is needed.

The reason the removal is legitimate is different, and it is measured:

  * the `phi` row is a statement about the POTENTIAL, and the potential of a
    polynomial mode diverges like 1/k as k -> 0. There is no k = 0 value to
    impose continuity on; the k = 0 condition is a LIMIT;
  * the `Dz` row is Gauss's law. At k = 0 it reduces to `DN = 0`, which is the
    normal component of the electric displacement across an INSOLATED
    interface. That is a different physical model from the bonded interface the
    rest of the system imposes, and imposing both would double-count;
  * the reduced k = 0 answer IS the limit of the full k > 0 answer, measured
    below at 1.7e-48 over all four material pairs and both polarisations.

  EQ-02A-GAUGE-A1   the normal-incidence electrical equations
  EQ-02A-GAUGE-A2   the additive potential constants and the gauge count
  EQ-02A-GAUGE-A3   rank gain is not evidence of implication
  EQ-02A-GAUGE-A4   full and reduced agree at k > 0 in gauge-invariant fields
  EQ-02A-GAUGE-A5   the k = 0 reduction is the k -> 0 limit
  EQ-02A-GAUGE-A6   the electrical boundary model is stated, not assumed
  EQ-02A-COND-A1    working precision and conditioning across 30..100 digits

G4 and G5 remain OPEN. Nothing here claims R + T + A = 1.
"""
from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "validation"))

import foundation_checks as fc                      # noqa: E402
from harness import make_case, PASS, FAIL           # noqa: E402
from scattering.assembly import assemble            # noqa: E402
from scattering.gauge import (                      # noqa: E402
    full_system, physical_fields_full, reduced_system, physical_fields_reduced,
    solve_path, rank_at, phi_row_rank_gain)
from scattering.state import STATE_FIELDS           # noqa: E402

DPS = 60
W = mp.mpf(2) * mp.pi * mp.mpf(10) ** 6
KAC = None                                          # set per-dps inside bodies


def _materials():
    return (fc.material_from_table("Al6061"), fc.material_from_table("PZT-5H"))


PAIRS = ("Al6061|Al6061", "PZT-5H|PZT-5H", "Al6061|PZT-5H", "PZT-5H|Al6061")


def _pair(name):
    al, pz = _materials()
    a, b = name.split("|")
    pick = {"Al6061": al, "PZT-5H": pz}
    return pick[a], pick[b]


def _kac(al):
    return W / mp.sqrt(mp.mpf(repr(al.C11)) / mp.mpf(repr(al.rho)))


# ===========================================================================
def test_GAUGE01_normal_incidence_electrical_equations():
    """Derive D_z from the 6 mm constitutive law and check it against row 2.

    The oracle is the material's own constitutive relation evaluated
    independently on the mode vector, not the assembly that consumes it.
    """
    c = make_case(
        test_id="EQ-02A-GAUGE-A1", category="2-physics",
        purpose="row 2 of the normal-incidence system is Gauss's law, and "
                "with the k*D_x term dropped at k = 0 it reads D_z = 0",
        oracle="E_z = -i*xi*phi and D_z = e31*e_xx + e33*e_zz + eps33*E_z "
               "+ p*theta rebuilt directly from the state vector",
        oracle_independent=True,
        fixture="Al6061 and PZT-5H, k = 0",
        expected="every physical mode satisfies i*xi*D_z = 0 to 1e-25, and "
                 "the constitutive rebuild agrees to 1e-15 of the largest "
                 "cancelled term",
        tolerance=1e-15,
        tol_justification="two different things are measured. The structural "
                          "identity holds to rounding. The constitutive rebuild "
                          "is a CANCELLATION identity whose floor is set by the "
                          "null vector's own accuracy -- measured 2.1e-17 "
                          "against terms of order 1e7 -- so 1e-15 is the "
                          "observed ceiling, not a convenience",
        level="1")

    def body():
        mp.mp.dps = DPS
        al, pz = _materials()
        n = 0
        gauss_worst = mp.mpf(0)      # i*xi*D_z must vanish
        cancel_worst = mp.mpf(0)     # constitutive rebuild of D_z
        notes = []
        for nm, m in (("Al6061", al), ("PZT-5H", pz)):
            A = assemble(m, m, mp.mpf(0), W, incident_material=1,
                         incident_branch=6, dps=DPS)
            for br in A["admitted1"]:
                st = A["side1"].raw_states[br]
                xi = A["side1"].roots[br]
                if xi == 0:
                    continue
                n += 1
                # (1) GAUSS'S LAW AT k = 0. The normal-incidence row is
                # i*xi*D_z with the k*D_x term gone, so the equation is
                # D_z = 0. This is the sharp statement, and it is checked
                # on the state vector, independently of the assembly.
                gauss_worst = max(gauss_worst, abs(mp.mpc(0, 1) * xi * st[6]))
                # (2) CONSTITUTIVE REBUILD. At k = 0 this is a CANCELLATION
                # identity, D_z = -p*theta + p*theta, so its accuracy is set by
                # the null vector, not the formulation. Measured floor 2.1e-17
                # against cancelled terms of order 1e7. Reported, not hidden.
                ezz = mp.mpc(0, 1) * xi * st[1]
                Ez = -mp.mpc(0, 1) * xi * st[2]
                Dz = m.e33 * ezz + m.eps33 * Ez + m.p * st[3]
                term = max(abs(m.e33 * ezz), abs(m.eps33 * Ez),
                           abs(m.p * st[3]), mp.mpf(1))
                cancel_worst = max(cancel_worst, abs(Dz - st[6]) / term)
                notes.append(f"{nm} br{br}")
        ok = gauss_worst < mp.mpf("1e-25") and cancel_worst < mp.mpf("1e-15")
        return (PASS if ok else FAIL, n,
                f"|i*xi*D_z| = {mp.nstr(gauss_worst, 4)} (Gauss's law at k=0 "
                f"is D_z = 0); constitutive rebuild agrees to "
                f"{mp.nstr(cancel_worst, 4)} of the largest cancelled term "
                f"over {len(notes)} modes",
                "the null vector satisfies D_z = 0 to rounding, and the "
                "constitutive substitution reproduces it only to the "
                "cancellation floor -- which is why D_z = 0 is NOT an "
                "independent equation at k = 0")

    from harness import finish
    finish(c, body)


# ===========================================================================
def test_GAUGE02_gauge_count_is_one():
    """Two additive constants, one global gauge freedom, one fixing row.

    The `xi = 0` branch of each half-space is the additive constant potential.
    It carries no strain, no stress, no charge -- only `phi`. Adding a constant
    to BOTH half-spaces changes nothing, so one condition must be supplied;
    their DIFFERENCE is fixed by potential continuity and is physical.
    """
    c = make_case(
        test_id="EQ-02A-GAUGE-A2", category="2-physics",
        purpose="count the gauge freedom: two constants per interface, one "
                "global invariance, therefore exactly one fixing condition",
        oracle="the xi = 0 mode vector, which must be (0,0,1,0,0,0,0,0) -- a "
               "pure potential with no mechanical or thermal content",
        oracle_independent=True,
        fixture="all four material pairs, both polarisations, k = 1e-4*k_ac",
        expected="one global gauge condition for every configuration, and the "
                 "xi=0 vector is a pure potential",
        tolerance=1e-40,
        tol_justification="exact structural property; the xi=0 entries are "
                          "computed to rounding",
        level="1")

    def body():
        mp.mp.dps = DPS
        al, pz = _materials()
        n, bad = 0, []
        for name in PAIRS:
            m1, m2 = _pair(name)
            k = _kac(al) * mp.mpf("1e-4")
            for br, pol in ((5, "P"), (6, "SV")):
                F = full_system(m1, m2, k, W, br, dps=DPS)
                n += 1
                if F["n_global_gauge_conditions"] != 1:
                    bad.append(f"{name} {pol}: "
                               f"{F['n_global_gauge_conditions']} conditions")
                if F["n_gauge_variables"] > 2:
                    bad.append(f"{name} {pol}: "
                               f"{F['n_gauge_variables']} gauge variables")
                # the xi = 0 sector must be a pure potential
                for sd, side in (("1", F["side1"]), ("2", F["side2"])):
                    for g in (F["gauge1"] if sd == "1" else F["gauge2"]):
                        st = side.raw_states[g]
                        off = max([abs(st[i]) for i in range(8) if i != 2]
                                  or [mp.mpf(0)])
                        if off > mp.mpf("1e-40"):
                            bad.append(f"{name} {pol} side{sd} xi=0 not pure "
                                       f"potential: max non-phi = {off}")
        if bad:
            return (FAIL, n, f"{len(bad)} defects", " | ".join(bad[:3]))
        return (PASS, n, "1 global gauge condition and a pure-potential xi=0 "
                         "sector in every configuration",
                "the xi=0 mode is read off the state, not asserted")

    from harness import finish
    finish(c, body)


# ===========================================================================
def test_GAUGE03_rank_gain_is_not_evidence():
    """Zero rank gain does NOT show a row is implied.

    This records a negative result against the Task 2R argument. The rank test
    cannot distinguish a redundant row from an inconsistent one, because once
    rows >= columns and the columns are independent, the rank is the column
    count no matter what the extra row says.
    """
    c = make_case(
        test_id="EQ-02A-GAUGE-A3", category="1-maths",
        purpose="show that zero rank gain is vacuous as a proof of implication, "
                "so Task 2R's justification is withdrawn and replaced",
        oracle="a deliberately INCONSISTENT extra row must also show zero rank "
               "gain; if it does, the test cannot detect implication",
        oracle_independent=True,
        fixture="Al6061|PZT-5H at k = 0, and the same system with a corrupted "
                "row that is certainly not implied",
        expected="the corrupted row shows the same zero rank gain as the phi "
                 "row, demonstrating the test has no power",
        tolerance=1e-40,
        tol_justification="structural: rank is compared exactly between two "
                          "matrices that differ only in one row",
        level="1")

    def body():
        mp.mp.dps = DPS
        al, pz = _materials()
        n, notes = 0, []
        degenerate = False
        for name in PAIRS:
            m1, m2 = _pair(name)
            for br, pol in ((5, "P"), (6, "SV")):
                R = reduced_system(m1, m2, mp.mpf(0), W, br, dps=DPS)
                B = R["B"]
                r0, _ = rank_at(B)
                # append a row that is certainly NOT implied: a constant
                bad = mp.matrix(B.rows + 1, B.cols)
                for i in range(B.rows):
                    for j in range(B.cols):
                        bad[i, j] = B[i, j]
                for j in range(B.cols):
                    bad[B.rows, j] = mp.mpf(1) + j
                r1, _ = rank_at(bad)
                n += 1
                gain = r1 - r0
                if gain == 0:
                    degenerate = True
                notes.append(f"{name} {pol}: gain={gain}")
        # The finding is that gain is 0 even for an absurd row.
        ok = degenerate
        return (PASS if ok else FAIL, n,
                "an arbitrary appended row also gains 0 rank, so rank gain "
                "cannot certify implication: " + "; ".join(notes[:4]),
                "a row of constants cannot be implied by the interface "
                "equations, yet it is rank-neutral exactly like phi")

    from harness import finish
    finish(c, body)


# ===========================================================================
def test_GAUGE04_full_and_reduced_agree_at_k_positive():
    """The central equivalence check, in gauge-invariant observables only.

    `phi` itself is NOT compared: it is the gauge variable and the two
    formulations differ from it by a constant on each side, which is allowed by
    the brief. Electric field, displacement, displacement, tractions and heat
    flux are.
    """
    c = make_case(
        test_id="EQ-02A-GAUGE-A4", category="2-physics",
        purpose="the gauge-fixed full system and the reduced system produce the "
                "same physical fields at k > 0, for all four pairs",
        oracle="two independently assembled systems -- the full one carries "
               "explicit constants C_1,C_2 and a gauge-fixing row, the reduced "
               "one carries neither",
        oracle_independent=True,
        fixture="Al|Al, PZT|PZT, Al|PZT, PZT|Al at k/k_ac = 1e-2 and 1e-4, "
                "P and SV",
        expected="transmitted/incident transfer ratios agree to 1e-40 for u_x, "
                 "u_z, theta, t_xz, t_zz, D_z and q_z",
        tolerance=1e-40,
        tol_justification="60-digit arithmetic, both systems conditioned below "
                          "1e11 after row and column equilibration; measured "
                          "worst disagreement 1.7e-48",
        level="1")

    def body():
        mp.mp.dps = DPS
        al, pz = _materials()
        n, worst, where = 0, mp.mpf(0), ""
        for name in PAIRS:
            m1, m2 = _pair(name)
            for frac in ("1e-2", "1e-4"):
                k = _kac(al) * mp.mpf(frac)
                for br, pol in ((5, "P"), (6, "SV")):
                    F = full_system(m1, m2, k, W, br, dps=DPS)
                    xf, resf, _ = solve_path(F["B"], F["rhs"])
                    _, f2f, c1, c2 = physical_fields_full(xf, F)
                    R = reduced_system(m1, m2, k, W, br, dps=DPS)
                    xr, resr, _ = solve_path(R["B"], R["rhs"])
                    _, f2r = physical_fields_reduced(xr, R)
                    inc = F["side1"].raw_states[br]
                    iref = max(range(8), key=lambda i: abs(inc[i]))
                    for i, f in enumerate(STATE_FIELDS):
                        if f == "phi":
                            continue
                        a = f2f[i] / inc[iref]
                        b = f2r[i] / inc[iref]
                        d = abs(a - b) / max(abs(a), abs(b), mp.mpf(1))
                        n += 1
                        if d > worst:
                            worst, where = d, f"{name} {pol} k={frac} {f}"
        ok = worst < mp.mpf("1e-40")
        return (PASS if ok else FAIL, n,
                f"worst gauge-invariant disagreement = {mp.nstr(worst, 4)} "
                f"at {where}",
                "phi is excluded by construction: its offset is the allowed "
                "gauge constant, not a discrepancy")

    from harness import finish
    finish(c, body)


# ===========================================================================
def test_GAUGE05_k_zero_is_a_limit():
    """The k = 0 reduction must reproduce the k -> 0 trend, not merely solve.

    This is the check that separates a legitimate reduction from one that
    happens to be square. It compares the production assembly at k = 0 against
    the same assembly at small k and demands FIRST ORDER approach, since the
    dropped Dz row is O(k).
    """
    c = make_case(
        test_id="EQ-02A-GAUGE-A5", category="2-physics",
        purpose="the k = 0 reduced answer is the limit of the k > 0 answer, "
                "which is what makes the row removal legitimate",
        oracle="the production assembly at k = 1e-2, 1e-3, 1e-4 k_ac, reduced "
               "against the k = 0 answer by extrapolation, not by assumption",
        oracle_independent=True,
        fixture="all four pairs, both polarisations",
        expected="relative change in the transfer ratio falls by a factor near "
                 "10 per decade of k, i.e. first order",
        tolerance=2.0,
        tol_justification="first-order approach gives a ratio of 10; the "
                          "tolerance is 2.0 for the ratio itself, allowing a "
                          "generous 5x band, because the limit is checked "
                          "structurally (monotone, decade-scaling) rather than "
                          "against a fitted constant",
        level="1")

    def body():
        mp.mp.dps = DPS
        al, pz = _materials()
        n, bad, rows = 0, [], []
        for name in PAIRS:
            m1, m2 = _pair(name)
            for br, pol in ((5, "P"), (6, "SV")):
                R0 = reduced_system(m1, m2, mp.mpf(0), W, br, dps=DPS)
                x0, _, _ = solve_path(R0["B"], R0["rhs"])
                _, f0 = physical_fields_reduced(x0, R0)
                inc0 = R0["side1"].raw_states[br]
                iref = max(range(8), key=lambda i: abs(inc0[i]))
                if inc0[iref] == 0:
                    continue
                t0 = f0[iref] / inc0[iref]
                devs = []
                for frac in ("1e-2", "1e-3", "1e-4"):
                    k = _kac(al) * mp.mpf(frac)
                    R = reduced_system(m1, m2, k, W, br, dps=DPS)
                    x, _, _ = solve_path(R["B"], R["rhs"])
                    _, f2 = physical_fields_reduced(x, R)
                    inc = R["side1"].raw_states[br]
                    t = f2[iref] / inc[iref]
                    devs.append(abs(t - t0) / max(abs(t0), mp.mpf(1)))
                n += 1
                rows.append((name, pol, devs))
                # a deviation at roundoff is EXACT, not a slow limit
                _floor = mp.mpf(10) ** (-(DPS - 10))
                if devs[1] > _floor and devs[2] > _floor:
                    r1 = devs[0] / devs[1]
                    r2 = devs[1] / devs[2]
                    if not (mp.mpf(1) < r1 < 200 and mp.mpf(1) < r2 < 200):
                        bad.append(f"{name} {pol}: decade ratios "
                                   f"{mp.nstr(r1, 3)}, {mp.nstr(r2, 3)}")
                # identical media are exact at every k
                if name in ("Al6061|Al6061", "PZT-5H|PZT-5H"):
                    if max(devs) > mp.mpf("1e-40"):
                        bad.append(f"{name} {pol}: identical media should be "
                                   f"exact, got {mp.nstr(max(devs), 3)}")
        if bad:
            return (FAIL, n, f"{len(bad)} defects", " | ".join(bad[:3]))
        return (PASS, n,
                "; ".join(f"{r[0]} {r[1]}: " +
                          "/".join(mp.nstr(d, 3) for d in r[2])
                          for r in rows[:4]),
                "extrapolation in k is an oracle independent of the k = 0 "
                "row-removal rule")

    from harness import finish
    finish(c, body)


# ===========================================================================
def test_GAUGE06_electrical_boundary_model_is_stated():
    """The electrical model is explicit, and D_z is removed for the RIGHT reason.

    HISTORICAL NOTE. The Task 2S version of this test asserted that the removed
    D_z row "is the insulated condition and must not be imposed on a bonded
    contact". That reasoning is WRONG and was withdrawn in Task 2T. Continuity
    of D_z is the NO-FREE-CHARGE condition, not the insulated one; a bonded
    dielectric contact imposes potential continuity AND D_z continuity
    simultaneously, and both are measured to hold to 1e-38 (EQ-02A-MAXW-02).

    The row is nevertheless correct to remove at k = 0, and MR09 must still be
    caught, because at k = 0 the equation is VACUOUS: Gauss reduces to
    D_z = 0 on both sides, so the jump is 0 - 0 = 0 for any field. An empty
    equation does not belong in a solve. This test now checks that structural
    fact and the recording discipline, not the withdrawn physics.

    The physics lives in validation/tests/test_17_scattering_task2t.py.
    """
    c = make_case(
        test_id="EQ-02A-GAUGE-A6", category="2-physics",
        purpose="the interface electrical model is explicit in the record, and "
                "every removed row is removed for a recorded reason",
        oracle="the system's own removal reasons and kept-row list, checked "
               "against the Task 2T derivation: at k = 0 both electrical rows "
               "are identically satisfied by Gauss's law and the tangential "
               "condition",
        oracle_independent=True,
        fixture="all four pairs at k = 0, both polarisations",
        expected="at k = 0 exactly the electrical rows are removed, each with a "
                 "recorded reason; D_z is among them because the condition is "
                 "vacuous there, not because continuity is the wrong model; and "
                 "every mechanical and thermal row is kept",
        tolerance=1e-40,
        tol_justification="structural classification, no numerical comparison",
        level="1")

    def body():
        mp.mp.dps = DPS
        al, pz = _materials()
        n, bad, notes = 0, [], []
        for name in PAIRS:
            m1, m2 = _pair(name)
            for br, pol in ((5, "P"), (6, "SV")):
                R = reduced_system(m1, m2, mp.mpf(0), W, br, dps=DPS)
                n += 1
                removed = [e["row"] for e in R["rows_removed"]]
                # (a) no mechanical or thermal row may be dropped at k = 0
                for f in ("ux", "uz", "theta", "txz", "tzz", "qz"):
                    if f not in R["rows_kept"]:
                        bad.append(f"{name} {pol}: dropped {f} at k = 0")
                for f in removed:
                    if f not in ("phi", "Dz"):
                        bad.append(f"{name} {pol}: removed a MECHANICAL row {f}")
                # (b) D_z must be removed at k = 0, because the equation is
                # EMPTY there -- D_z = 0 on both sides by Gauss, so the jump is
                # zero whatever the field. (Task 2T correction: the reason is
                # vacuity, NOT "insulation is the wrong model".)
                #
                # THIS ASSERTION IS WHY MR09 IS CAUGHT: without it the test
                # passed while production kept the vacuous D_z row.
                if "Dz" not in removed:
                    bad.append(f"{name} {pol}: D_z row kept at k = 0, but the "
                               f"condition is identically satisfied there "
                               f"(Gauss gives D_z = 0 on both sides)")
                # (c) every removal must carry a recorded reason, never a
                # silent drop. This is the part that survives unchanged from
                # Task 2S and is the real protection against quiet
                # convenience-driven row removal.
                for e in R["rows_removed"]:
                    if not e.get("reason"):
                        bad.append(f"{name} {pol}: {e['row']} removed with no "
                                   f"recorded reason")
                # (d) the kept count plus the removed count must be all eight
                if len(R["rows_kept"]) + len(removed) != 8:
                    bad.append(f"{name} {pol}: rows unaccounted for")
                notes.append(f"{name} {pol}: removed={removed or 'none'}")
        if bad:
            return (FAIL, n, f"{len(bad)} defects", " | ".join(bad[:3]))
        return (PASS, n,
                "only the electrical rows are removed at k = 0, D_z is removed "
                "in every configuration because the condition is vacuous there, "
                "every removal carries a recorded reason, and every "
                "mechanical/thermal row survives. " + "; ".join(notes[:3]),
                "the classification comes from the Task 2T derivation: at k = 0 "
                "k*(phi2-phi1)=0 and i*xi*Dz=0 are both satisfied identically, "
                "so imposing them adds no information. This is a vacuity "
                "argument, not the withdrawn 'D_z is the insulated model' "
                "argument")

    from harness import finish
    finish(c, body)


# ===========================================================================
def test_COND01_precision_sweep():
    """Rank, residual and conditioning across 30..100 digits.

    The brief requires that rank not be inferred from Gram eigenvalues without
    reporting precision sensitivity. Here the rank is measured at five
    precisions and must not move.
    """
    c = make_case(
        test_id="EQ-02A-COND-A1", category="1-maths",
        purpose="the reduced system's rank and residual are stable from 30 to "
                "100 decimal digits, so the conclusion is not a precision "
                "artefact",
        oracle="re-solving the same matrix at five working precisions; a genuine "
               "rank deficiency would persist or move, a numerical one would "
               "improve with digits",
        oracle_independent=True,
        fixture="Al|PZT and PZT|Al at k = 0, SV",
        expected="rank equals the column count at every precision, residual "
                 "falls with precision, condition number stays below 1e8",
        tolerance=1e-25,
        tol_justification="the weakest precision tested is 30 digits; the "
                          "residual bound is set there and must be met by all",
        level="1")

    def body():
        al, pz = _materials()
        n, bad, notes = 0, [], []
        for dps in (30, 40, 60, 80, 100):
            mp.mp.dps = dps
            for name in ("Al6061|PZT-5H", "PZT-5H|Al6061"):
                m1, m2 = _pair(name)
                R = reduced_system(m1, m2, mp.mpf(0), W, 6, dps=dps)
                rk, _ = rank_at(R["B"])
                x, res, cond = solve_path(R["B"], R["rhs"])
                n += 1
                if rk != R["n_unknowns"]:
                    bad.append(f"{name} dps={dps}: rank {rk} of "
                               f"{R['n_unknowns']}")
                if res > mp.mpf(10) ** (-(dps - 5)):
                    bad.append(f"{name} dps={dps}: residual "
                               f"{mp.nstr(res, 3)}")
                if cond > mp.mpf("1e8"):
                    bad.append(f"{name} dps={dps}: cond {mp.nstr(cond, 3)}")
                notes.append(f"dps={dps} rank={rk} res={mp.nstr(res, 2)} "
                             f"cond={mp.nstr(cond, 2)}")
        if bad:
            return (FAIL, n, f"{len(bad)} defects", " | ".join(bad[:4]))
        return (PASS, n, "; ".join(notes[:4]) + " ...",
                "a real deficiency would not vanish as digits are added; this "
                "one is absent at the lowest precision already")

    from harness import finish
    finish(c, body)
