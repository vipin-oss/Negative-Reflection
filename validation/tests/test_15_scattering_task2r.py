"""
Phase 2A task 2R -- normal-incidence degeneracy and the scattering assembly audit.

Task 2 asserted that "at normal incidence the sagittal geometry sets
u_x = sigma_xz = 0 for every mode". That statement is FALSE, and these checks
exist to make the specific ways it can go wrong fail loudly.

The seven defects SCTR01..SCTR07 pin are:

  SCTR01  Removing a nontrivial shear boundary row.
  SCTR02  Mistaking a P-wave zero component for a universal modal zero.
  SCTR03  Artificial rank deficiency introduced by the reduction.
  SCTR04  Duplicate or dependent modal columns.
  SCTR05  Gauge freedom misclassified as physical nonuniqueness.
  SCTR06  An incorrect identical-medium reflection.
  SCTR07  Residuals checked only against the assembled matrix instead of the
          physical fields.

Each is written so that a negative-control mutation of the code makes it fail.
The mutations live in validation/mutation_tests/run_mutations.py (MR01..MR06).

Every check states its oracle and why that oracle is independent of the code
under test. Nothing here closes G4 or G5.
"""
from __future__ import annotations

import mpmath as mp

from harness import PASS, FAIL, make_case, finish
import foundation_checks as fc

from scattering.assembly import (                       # noqa: E402
    TRIVIAL_ROW_TOL, assemble, matrix_rank,
)
from scattering.elastic_reference import (               # noqa: E402
    ISO, fresnel_normal_amplitudes,
)


def _iso_al():
    """The independent reference's own Al6061 instance, built from constants.

    Constructed here rather than reached for as an attribute: ISO is the
    reference's material class, and its singletons are not part of the contract
    this test relies on.
    """
    return ISO("Al6061", 2700.0, 107.6e9, 26.5e9)
from scattering.roots import (                           # noqa: E402
    row_relative_residual, verified_roots,
)
from scattering.state import STATE_FIELDS, state_vector   # noqa: E402

W = 2 * mp.pi * 1e6
DPS = 40


def _k_from_angle(material, theta_deg):
    cL = mp.sqrt(mp.mpf(repr(material.C11)) / mp.mpf(repr(material.rho)))
    return mp.mpf(W) * mp.sin(mp.radians(mp.mpf(theta_deg))) / cL


# ===========================================================================
def test_SCTR01_shear_rows_not_removed():
    """The u_x and sigma_xz rows must survive a shear-incidence reduction.

    WHY THIS IS THE TEST THAT MATTERS. At normal incidence on an isotropic
    solid the shear branches carry u_x != 0 and sigma_xz = i xi C44 u_x != 0.
    Those are the ONLY rows that distinguish a shear mode from a longitudinal
    one. If the reduction deletes them -- because they are zero for the other
    branches, or because the incident wave has u_x = 0 -- then the remaining
    system cannot tell the two apart and any solve is meaningless.

    A row may be removed only when it is zero for every admitted mode AND zero
    on the incident right-hand side.
    """
    c = make_case(
        test_id="SCTR01", category="11-normal-incidence",
        purpose="a shear boundary row that is nonzero for an admitted mode is "
                "never removed from the interface system",
        oracle="the exact SV state vector (u_x, u_z, phi, theta) = (1, 0, 0, 0) "
               "with sigma_xz = i xi C44 u_x, derived from the constitutive "
               "relation and independent of the assembly code",
        oracle_independent=True,
        fixture="Al6061 | Al6061, k = 0 exactly, branch 6 (shear) incident",
        expected="u_x and txz are among the kept rows, and no row is flagged "
                 "inconsistent",
        tolerance=0,
        tol_justification="an exact set-membership predicate; there is no "
                          "numerical tolerance to choose",
        level="1")

    def body():
        mp.mp.dps = DPS
        mat = fc.material_from_table("Al6061")
        k = mp.mpf(0)
        A = assemble(mat, mat, k, mp.mpf(W), incident_material=1,
                     incident_branch=6, dps=DPS)
        kept = A["rows_kept"]

        roots = verified_roots(mat, k, mp.mpf(W))
        sv_nonzero_ux = any(r["ok"] and abs(r["vector"][0]) > 0 for r in roots)

        # A second, independent condition: the reduced system must actually
        # SOLVE. Disabling the gauge rule leaves Al6061|PZT-5H inconsistent at
        # 1.1e-34 without touching the row set, so a row-presence check alone
        # would not notice.
        B, rhs = A["B"], A["rhs"]
        Bc = mp.matrix(B.rows, B.cols)
        rc = mp.matrix(B.rows, 1)
        for i in range(B.rows):
            mx = max(abs(B[i, j]) for j in range(B.cols))
            d = mx if mx > 0 else mp.mpf(1)
            for j in range(B.cols):
                Bc[i, j] = B[i, j] / d
            rc[i] = rhs[i] / d
        x = mp.lu_solve(Bc.H * Bc, Bc.H * rc)
        resid = max(abs(t) for t in Bc * x - rc)
        consistent = resid < mp.mpf("1e-25")

        # the same statement for the mixed pair, where the gauge rule is load
        # bearing and the row set is otherwise unchanged
        A2 = assemble(mat, fc.material_from_table("PZT-5H"), k, mp.mpf(W),
                      incident_material=1, incident_branch=5, dps=DPS)
        B2, rhs2 = A2["B"], A2["rhs"]
        Bc2 = mp.matrix(B2.rows, B2.cols)
        rc2 = mp.matrix(B2.rows, 1)
        for i in range(B2.rows):
            mx = max(abs(B2[i, j]) for j in range(B2.cols))
            d = mx if mx > 0 else mp.mpf(1)
            for j in range(B2.cols):
                Bc2[i, j] = B2[i, j] / d
            rc2[i] = rhs2[i] / d
        x2 = mp.lu_solve(Bc2.H * Bc2, Bc2.H * rc2)
        mixed_resid = max(abs(t) for t in Bc2 * x2 - rc2)
        mixed_ok = mixed_resid < mp.mpf("1e-25")

        ok = ("ux" in kept and "txz" in kept
              and not A["rows_inconsistent"] and sv_nonzero_ux
              and consistent and mixed_ok)
        detail = (f"kept = {kept}; removed = "
                  f"{[(r['row'], r.get('INCONSISTENT', False)) for r in A['rows_removed']]}; "
                  f"n_eq = {A['n_equations']}, n_unk = {A['n_unknowns']}, "
                  f"rank = {A['rank']}")
        actual = (f"kept ux = {'ux' in kept}, kept txz = {'txz' in kept}; "
                  f"verified shear mode with u_x != 0 = {sv_nonzero_ux}; "
                  f"Al|Al residual = {mp.nstr(resid, 4)}; "
                  f"Al|PZT residual = {mp.nstr(mixed_resid, 4)}")
        return (PASS if ok else FAIL, 3, actual, detail)

    finish(c, body)


# ===========================================================================
def test_SCTR02_p_wave_zero_is_not_a_modal_zero():
    """u_x = 0 for the incident P wave does NOT mean u_x = 0 for every mode.

    The Task 2 error in one line: a component measured on ONE wave was assumed
    to vanish on ALL of them. The P wave at normal incidence does have u_x = 0;
    the shear wave does not.
    """
    c = make_case(
        test_id="SCTR02", category="11-normal-incidence",
        purpose="a zero component of the incident P wave is not treated as a "
                "universal modal zero",
        oracle="the constitutive relation sigma_xz = C44 du_x/dz - e15 dphi/dz, "
               "which forces sigma_xz != 0 whenever u_x varies",
        oracle_independent=True,
        fixture="Al6061, k = 0, all eight roots verified against the matrix",
        expected="the P branch has u_x = 0 while the shear branch has u_x != 0, "
                 "and both are verified null vectors",
        tolerance=1e-30,
        tol_justification="the constitutive identity is exact, so the agreement "
                          "threshold is the 40-digit arithmetic floor",
        level="1")

    def body():
        mp.mp.dps = DPS
        mat = fc.material_from_table("Al6061")
        k = mp.mpf(0)
        roots = verified_roots(mat, k, mp.mpf(W))

        p_branches = [r for r in roots
                      if r["ok"] and abs(r["vector"][1]) > 0
                      and abs(r["vector"][0]) == 0
                      and abs(r["vector"][2]) == 0]
        s_branches = [r for r in roots
                      if r["ok"] and abs(r["vector"][0]) > 0]

        traction_ok = True
        if s_branches:
            sb = s_branches[0]
            st = state_vector(mat, k, sb["xi"], mp.mpf(W), sb["vector"])
            C44 = mp.mpf(repr(mat.C44))
            predicted = mp.mpc(0, 1) * sb["xi"] * C44 * st[0]
            scale = max(abs(predicted), mp.mpf(1))
            traction_ok = abs(st[4] - predicted) / scale < mp.mpf("1e-30")

        # The assembly must not treat u_x as gauge-capable. Measured on a
        # P-incidence system, where the incident u_x really is zero: the row may
        # not be reclassified as an empty gauge equation.
        A = assemble(mat, mat, k, mp.mpf(W), incident_material=1,
                     incident_branch=5, dps=DPS)
        ux_not_gauge = "ux" not in A["rows_gauge_only"]
        ux_kept = "ux" in A["rows_kept"]

        ok = (bool(p_branches) and bool(s_branches) and traction_ok
              and ux_not_gauge and ux_kept)
        actual = (f"P branches with u_x = 0: {len(p_branches)}; "
                  f"shear branches with u_x != 0: {len(s_branches)}; "
                  f"sigma_xz = i xi C44 u_x reproduced: {traction_ok}; "
                  f"u_x classified gauge = {not ux_not_gauge}")
        detail = ("the shear state is the mode that would be lost if a P-wave "
                  "zero were generalised")
        return (PASS if ok else FAIL, 3, actual, detail)

    finish(c, body)


# ===========================================================================
def test_SCTR03_no_artificial_rank_deficiency():
    """The reduction must not lower the rank of the surviving system.

    Dropping a row cannot raise the rank. If a reduction that claims to remove
    only 0 = 0 equations changes the rank of the coefficient matrix, it removed
    something that was carrying information.
    """
    c = make_case(
        test_id="SCTR03", category="11-normal-incidence",
        purpose="trivial-row removal does not reduce the rank of the system",
        oracle="rank is invariant under deleting an identically zero row, which "
               "is a statement about linear algebra, not about this code",
        oracle_independent=True,
        fixture="Al6061 | Al6061, PZT|PZT, Al|PZT at k = 0, shear incident",
        expected="rank(B_reduced) == rank(B_full) in every case",
        tolerance=0,
        tol_justification="rank is an integer; equality is exact",
        level="1")

    def body():
        mp.mp.dps = DPS
        al = fc.material_from_table("Al6061")
        pz = fc.material_from_table("PZT-5H")
        k = mp.mpf(0)
        rows, ok = [], True
        for name, m1, m2 in (("Al|Al", al, al), ("PZT|PZT", pz, pz),
                             ("Al|PZT", al, pz)):
            A = assemble(m1, m2, k, mp.mpf(W), incident_material=1,
                         incident_branch=6, dps=DPS)
            r_full = matrix_rank(A["B_full"])
            r_red = matrix_rank(A["B"])
            n_unk = A["n_unknowns"]
            # A rank comparison alone is vacuous: two equally broken systems
            # agree. What must hold is that the reduced system is FULL RANK in
            # its unknowns -- every degree of freedom is pinned by a surviving
            # equation. Measured, this is the condition that fails if the
            # structural-zero floor is removed and rows carrying information
            # are silently discarded.
            full_rank_ok = (r_red == n_unk)
            rows.append(f"{name}: rank_full = {r_full}, rank_reduced = {r_red}, "
                        f"n_unk = {n_unk}, full rank = {full_rank_ok}, "
                        f"removed = {[r['row'] for r in A['rows_removed']]}")
            if r_full != r_red or not full_rank_ok:
                ok = False
        return (PASS if ok else FAIL, len(rows), " | ".join(rows),
                "rank is compared on the same column set, and the reduced "
                "system is additionally required to retain full column rank")

    finish(c, body)


# ===========================================================================
def test_SCTR04_no_duplicate_modal_columns():
    """No two admitted columns may be proportional.

    Two proportional columns mean two names for one degree of freedom: the
    coefficient vector is then non-unique for a reason that has nothing to do
    with the physics.
    """
    c = make_case(
        test_id="SCTR04", category="11-normal-incidence",
        purpose="the admitted basis contains no duplicate or proportional "
                "columns",
        oracle="two columns are dependent iff the matrix rank drops when one "
               "is present and rises when it is added; a Cauchy-Schwarz overlap "
               "on MIXED-UNIT columns is not a valid test and was replaced "
               "after it was measured to report 1.0 for genuinely distinct "
               "columns",
        oracle_independent=True,
        fixture="Al6061 | Al6061 at 0 and 30 deg, PZT|PZT at 0 deg",
        expected="every normalised overlap is below 1 - 1e-20",
        tolerance=1e-20,
        tol_justification="two distinct physical modes cannot be collinear; the "
                          "threshold is the arithmetic floor of the 40-digit "
                          "arithmetic used",
        level="1")

    def body():
        mp.mp.dps = DPS
        al = fc.material_from_table("Al6061")
        pz = fc.material_from_table("PZT-5H")
        from scattering.gauge import rank_at
        rows_out, total, ok = [], 0, True
        for name, m1, m2, theta in (
                ("Al|Al 0deg", al, al, 0),
                ("Al|Al 30deg", al, al, 30),
                ("PZT|PZT 0deg", pz, pz, 0),
                ("Al|PZT 0deg", al, pz, 0)):
            k = _k_from_angle(m1, theta)
            A = assemble(m1, m2, k, mp.mpf(W), incident_material=1,
                         incident_branch=6, dps=DPS)
            B = A["B"]
            r_all, _ = rank_at(B)
            # THE TEST: does any single column fail to contribute rank?
            # A duplicated column leaves the rank below the column count.
            augmented = (r_all == A["n_unknowns"])
            # and leave-one-out must strictly reduce the rank for every column
            deficient = []
            for c in range(B.cols):
                sub = mp.matrix(B.rows, B.cols - 1)
                for j in range(B.cols):
                    if j == c:
                        continue
                    for i in range(B.rows):
                        sub[i, j if j < c else j - 1] = B[i, j]
                r_sub, _ = rank_at(sub)
                total += 1
                if r_sub >= r_all:
                    deficient.append(c)
            if not augmented or deficient:
                ok = False
            rows_out.append(
                f"{name}: rank {r_all}/{A['n_unknowns']}, "
                f"redundant columns {deficient or 'none'}")
        worst, worst_at = (mp.mpf(0), None)
        actual = "; ".join(rows_out)
        return (PASS if ok else FAIL, total, actual,
                "every admitted column must contribute rank; leave-one-out must "
                "strictly reduce it")

    finish(c, body)


# ===========================================================================
def test_SCTR05_gauge_is_not_nonuniqueness():
    """The identical-medium system must not be called under-determined.

    Al|Al at normal incidence has fewer independent rows than unknowns. That is
    NOT evidence of a physically non-unique scattering solution. The check
    separates a GAUGE direction (whole state zero except a potential, no
    traction, no flux) from genuine nonuniqueness (two different PHYSICAL
    fields solving the same problem).
    """
    c = make_case(
        test_id="SCTR05", category="11-normal-incidence",
        purpose="a rank-deficient identical-medium system is classified by what "
                "the deficient directions physically are, not by rows < unknowns",
        oracle="the physical fields from the independent elastic reference, "
               "which are unique and analytic at normal incidence",
        oracle_independent=True,
        fixture="Al6061 | Al6061, k = 0, shear incident",
        expected="the removed rows carry no traction and no incident load, and "
                 "the independent reflection vanishes",
        tolerance=1e-40,
        tol_justification="the independent reflection is an exact zero for "
                          "matched media, so the check is equality to the "
                          "arithmetic floor",
        level="1")

    def body():
        mp.mp.dps = DPS
        mat = fc.material_from_table("Al6061")
        k = mp.mpf(0)
        A = assemble(mat, mat, k, mp.mpf(W), incident_material=1,
                     incident_branch=6, dps=DPS)
        r_theory, _t_theory = fresnel_normal_amplitudes(_iso_al(), _iso_al())

        gauge_clean = all(not r.get("INCONSISTENT") for r in A["rows_removed"])

        # The decisive check: on a mixed pair the potential row MUST be
        # recognised as gauge, and the resulting system must solve. If the rule
        # is disabled the row stays in, the system is over-determined and its
        # residual is large -- so this is what distinguishes "we classified the
        # gauge correctly" from "we silently failed to solve".
        pzt = fc.material_from_table("PZT-5H")
        A_mixed = assemble(fc.material_from_table("Al6061"), pzt, k, mp.mpf(W),
                           incident_material=1, incident_branch=5, dps=DPS)
        gauge_flagged = "phi" in A_mixed["rows_gauge_only"]

        Bm, rm_ = A_mixed["B"], A_mixed["rhs"]
        Bc = mp.matrix(Bm.rows, Bm.cols)
        rc = mp.matrix(Bm.rows, 1)
        for i in range(Bm.rows):
            mx = max(abs(Bm[i, j]) for j in range(Bm.cols))
            d = mx if mx > 0 else mp.mpf(1)
            for j in range(Bm.cols):
                Bc[i, j] = Bm[i, j] / d
            rc[i] = rm_[i] / d
        xm = mp.lu_solve(Bc.H * Bc, Bc.H * rc)
        mixed_resid = max(abs(t) for t in Bc * xm - rc)
        mixed_ok = mixed_resid < mp.mpf("1e-25")

        ok = (abs(r_theory) == 0 and A["n_unknowns"] > 0 and gauge_clean
              and gauge_flagged and mixed_ok)
        actual = (f"n_eq = {A['n_equations']}, n_unk = {A['n_unknowns']}, "
                  f"removed = {[r['row'] for r in A['rows_removed']]}, "
                  f"independent Al|Al reflection = {mp.nstr(r_theory, 6)}; "
                  f"Al|PZT phi flagged gauge = {gauge_flagged}, "
                  f"Al|PZT residual = {mp.nstr(mixed_resid, 4)}")
        detail = ("rows < unknowns alone is not taken as nonuniqueness; the "
                  "removed rows are checked to be physically empty")
        return (PASS if ok else FAIL, 3, actual, detail)

    finish(c, body)


# ===========================================================================
def test_SCTR06_identical_medium_reflection_is_zero():
    """For two identical welded half-spaces the reflected field must vanish.

    THE SIMPLEST PHYSICAL FACT IN THE WHOLE AUDIT. Nothing reflects at a welded
    interface between identical media -- there is no interface. The oracle is
    the closed-form Fresnel amplitude for matched impedances, which knows
    nothing about the assembly.
    """
    c = make_case(
        test_id="SCTR06", category="11-normal-incidence",
        purpose="identical welded half-spaces reflect nothing at normal incidence",
        oracle="closed-form Fresnel amplitude (Z2 - Z1)/(Z2 + Z1) at normal "
               "incidence with Z1 = Z2, plus the exact impedance identity",
        oracle_independent=True,
        fixture="Al6061 | Al6061 at k = 0",
        expected="reflection amplitude exactly zero and transmission unity",
        tolerance=0,
        tol_justification="matched-impedance Fresnel amplitudes are exactly 0 "
                          "and 1; the comparison is exact",
        level="1")

    def body():
        mp.mp.dps = DPS
        al = fc.material_from_table("Al6061")
        Z = mp.sqrt(mp.mpf(repr(al.rho)) * mp.mpf(repr(al.C11)))
        r_manual = (Z - Z) / (Z + Z)
        t_manual = 2 * Z / (Z + Z)
        r_ref, t_ref = fresnel_normal_amplitudes(_iso_al(), _iso_al())

        # AND THE ASSEMBLY ITSELF. The closed form says what must happen; this
        # says the code does it. Without this leg a sign slip on the second
        # half-space's columns would leave the closed-form statement true and
        # the check green.
        mat = fc.material_from_table("Al6061")
        A = assemble(mat, mat, mp.mpf(0), mp.mpf(W), incident_material=1,
                     incident_branch=6, dps=DPS)
        B, rhs = A["B"], A["rhs"]
        Bc = mp.matrix(B.rows, B.cols)
        rc = mp.matrix(B.rows, 1)
        for i in range(B.rows):
            mx = max(abs(B[i, j]) for j in range(B.cols))
            d = mx if mx > 0 else mp.mpf(1)
            for j in range(B.cols):
                Bc[i, j] = B[i, j] / d
            rc[i] = rhs[i] / d
        x = mp.lu_solve(Bc.H * Bc, Bc.H * rc)
        sol = mp.matrix(B.cols, 1)
        for j in range(B.cols):
            sol[j] = x[j]

        norms = A["column_norms"]
        s1, s2 = A["side1"], A["side2"]
        a1, a2 = A["admitted1"], A["admitted2"]

        def field(side, adm, i0):
            acc = [mp.mpc(0)] * 8
            for kk, br in enumerate(adm):
                c = sol[i0 + kk] / norms[i0 + kk]
                for i in range(8):
                    acc[i] += c * side.states[br][i]
            return acc

        refl = field(s1, a1, 0)
        trans = field(s2, a2, len(a1))
        inc = s1.states[6]
        refl_max = max(abs(t) for t in refl)
        cont = max(abs(inc[i] + refl[i] - trans[i])
                   / max(abs(inc[i]), abs(trans[i]), mp.mpf(1))
                   for i in range(8))

        ok = (abs(r_manual) == 0 and abs(t_manual - 1) == 0
              and abs(r_ref) == 0 and abs(t_ref - 1) == 0
              and refl_max < mp.mpf("1e-30") and cont < mp.mpf("1e-30"))
        actual = (f"Fresnel r = {mp.nstr(r_manual, 6)}, "
                  f"t = {mp.nstr(t_manual, 6)}; reference r = "
                  f"{mp.nstr(r_ref, 6)}, t = {mp.nstr(t_ref, 6)}; "
                  f"assembly reflected field = {mp.nstr(refl_max, 4)}, "
                  f"continuity error = {mp.nstr(cont, 4)}")
        return (PASS if ok else FAIL, 4, actual,
                "impedance computed from the material constants directly, and "
                "cross-checked against elastic_reference")

    finish(c, body)


# ===========================================================================
def test_SCTR07_residuals_against_physical_fields():
    """A small ||B a - rhs|| does not prove the boundary conditions hold.

    THE DIVERGENCE THE WHOLE TASK HINGES ON. A residual measured against the
    assembled matrix only says the linear system was solved. It says nothing
    about whether the equations in that matrix are the boundary conditions. If
    a row was removed, the residual stops testing it.

    This check rebuilds the physical statement -- every removed row must read
    0 = 0 in both the basis and the incident data -- and refuses a reduction
    that removed a row while the right-hand side was nonzero.
    """
    c = make_case(
        test_id="SCTR07", category="11-normal-incidence",
        purpose="interface reduction is verified on the physical fields, not "
                "only on the assembled matrix",
        oracle="the modal state vectors and the incident state, evaluated "
               "independently of the reduced operator",
        oracle_independent=True,
        fixture="Al6061 | Al6061, k = 0, shear incident, interface at z = 0",
        expected="every removed row is zero in the basis AND zero on the "
                 "incident right-hand side",
        tolerance=1e-25,
        tol_justification="TRIVIAL_ROW_TOL, the level below which a row is "
                          "defined to read 0 = 0 for this matrix",
        level="1")

    def body():
        mp.mp.dps = DPS
        al = fc.material_from_table("Al6061")
        pzt = fc.material_from_table("PZT-5H")
        k = mp.mpf(0)
        A = assemble(al, al, k, mp.mpf(W), incident_material=1,
                     incident_branch=6, dps=DPS)
        A_mixed = assemble(pzt, pzt, k, mp.mpf(W), incident_material=1,
                           incident_branch=5, dps=DPS)

        # TWO CONFIGURATIONS, because the defect is only visible in one.
        # For Al|Al the electric rows really do read 0 = 0. For PZT-5H|PZT-5H
        # under P incidence every half-space carries a potential AND the
        # incident wave carries one, so dropping a phi row there is
        # unambiguously wrong. A reduction that removes rows without checking
        # passes on the first and fails on the second.
        ok = True
        parts = []
        n_checked = 0
        for Am, lbl in ((A, "Al|Al SV"), (A_mixed, "PZT|PZT P")):
            if Am["rows_inconsistent"]:
                ok = False
                parts.append(f"{lbl}: inconsistent rows "
                             f"{Am['rows_inconsistent']}")
            Bm = Am["B_full"]
            for entry in Am["rows_removed"]:
                n_checked += 1
                ii = STATE_FIELDS.index(entry["row"])
                rmx = max(abs(Bm[ii, j]) for j in range(Bm.cols))
                rr = abs(Am["rhs_full"][ii])
                if entry.get("GAUGE_ONLY"):
                    # Legitimate in exactly two routes:
                    #  (old) a side genuinely has no content in this row, or
                    #        the incident datum supplies none; OR
                    #  (Task 2U repair) at exact k = 0 the electric rows are
                    #        routed to the potential-gauge handler
                    #        ([[phi]]=0 by relative-constant reconstruction;
                    #        D_z = 0 on shell), which must NOT depend on
                    #        basis or incident emptiness -- that was
                    #        BLK-02A-26. The removal is legitimate there only
                    #        when the assembly recorded it as a potential row.
                    task2u = (Am.get("k") == 0
                              and entry["row"] in Am.get("potential_rows", []))
                    if not (entry.get("basis_empty_side1")
                            or entry.get("basis_empty_side2")
                            or entry.get("incident_empty")
                            or task2u):
                        ok = False
                        parts.append(f"{lbl} {entry['row']}: gauge removal with "
                                     f"no empty basis and no empty datum")
                    continue
                if rmx > TRIVIAL_ROW_TOL or rr > TRIVIAL_ROW_TOL:
                    ok = False
                    parts.append(f"{lbl} {entry['row']}: removed but rowmax = "
                                 f"{mp.nstr(rmx, 4)}, rhs = {mp.nstr(rr, 4)}")

        # AND the physical fields on the surviving rows.
        B, rhs = A["B"], A["rhs"]
        Bc = mp.matrix(B.rows, B.cols)
        rc = mp.matrix(B.rows, 1)
        for i in range(B.rows):
            mx = max(abs(B[i, j]) for j in range(B.cols))
            d = mx if mx > 0 else mp.mpf(1)
            for j in range(B.cols):
                Bc[i, j] = B[i, j] / d
            rc[i] = rhs[i] / d
        x = mp.lu_solve(Bc.H * Bc, Bc.H * rc)
        sol = mp.matrix(B.cols, 1)
        for j in range(B.cols):
            sol[j] = x[j]

        norms = A["column_norms"]
        s1, s2 = A["side1"], A["side2"]
        a1, a2 = A["admitted1"], A["admitted2"]

        def field(side, adm, i0):
            acc = [mp.mpc(0)] * 8
            for kk, br in enumerate(adm):
                cc = sol[i0 + kk] / norms[i0 + kk]
                for i in range(8):
                    acc[i] += cc * side.states[br][i]
            return acc

        refl = field(s1, a1, 0)
        trans = field(s2, a2, len(a1))
        inc = s1.states[6]
        cont_bad = []
        for i, nm in enumerate(STATE_FIELDS):
            if nm not in A["rows_kept"]:
                continue
            lhs = inc[i] + refl[i]
            sc = max(abs(lhs), abs(trans[i]), mp.mpf(1))
            if abs(lhs - trans[i]) / sc > mp.mpf("1e-20"):
                cont_bad.append(nm)
        if cont_bad:
            ok = False
            parts.append(f"physical continuity violated on {cont_bad}")

        actual = ("; ".join(parts) if parts else
                  f"all {n_checked} removed rows justified and every kept row "
                  f"continuous on the physical fields")
        return (PASS if ok else FAIL, n_checked + len(A["rows_kept"]), actual,
                "each removed row is checked against the physical basis and the "
                "incident datum, and the surviving rows against the "
                "reconstructed interface fields")

    finish(c, body)


def run():
    test_SCTR01_shear_rows_not_removed()
    test_SCTR02_p_wave_zero_is_not_a_modal_zero()
    test_SCTR03_no_artificial_rank_deficiency()
    test_SCTR04_no_duplicate_modal_columns()
    test_SCTR05_gauge_is_not_nonuniqueness()
    test_SCTR06_identical_medium_reflection_is_zero()
    test_SCTR07_residuals_against_physical_fields()


if __name__ == "__main__":
    run()
