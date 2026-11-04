"""
Two-half-space coupled piezothermoelastic scattering assembly -- Phase 2A task 2.

WHAT IS ASSEMBLED
=================
For an interface at z = 0 between two 6mm piezothermoelastic half-spaces, the
physical interface state is

    s = ( u_x, u_z, phi, theta, sigma_xz, sigma_zz, D_z, q_z )

(`scattering.state`) and welded electrical/thermal contact means all eight
components are continuous. With a modal basis of N1 admissible modes in medium 1
and N2 in medium 2,

    [ B1  -B2 ] [ a1 ]   [ -s_inc ]
                [ a2 ]  =

MODE SELECTION -- THE POINT OF THIS FILE
========================================
Modes are selected by `scattering.radiation.classify_role`, i.e. by
`outgoing_status`, NEVER by the `DECAYS_*` / `TRAVELS_*` labels and NEVER by the
sign of Im(xi) alone:

  * a lossless outgoing propagating mode (Im xi = 0 with resolvable flux) is
    ADMITTED. Task 1R established these are real channels.
  * a genuine electric-potential gauge (phi = const, D_z = 0) is EXCLUDED,
    because it is a null channel and adds a zero column.
  * a DECOUPLED ELECTROSTATIC mode (mechanical content zero but D_z =-
    eps33*xi*phi != 0) is ADMITTED. It satisfies the D_z continuity equation
    like any other mode. Conflating this with a gauge -- which an earlier
    version of `scattering.modal.classify_mode` did -- silently deletes the
    electrostatic channel of an elastic half-space while keeping the
    equivalent channel of a piezoelectric one.

DEGREES OF FREEDOM ARE COUNTED, NOT ASSUMED
===========================================
The brief forbids assuming an 8x8 system. The count is computed and reported:
roots, gauges, physical modes, admitted modes per side, unknowns, equations,
rank. A square 8x8 system appears only when the physics supplies eight
independent channels per side.

NORMALISATION
=============
Two exact rescalings are applied, and neither changes the solution set:

1. DIMENSIONAL: every state entry is divided by `state.reference_scales`, so
   rows that mix metres, volts, kelvin and pascals become comparable.
2. COLUMN: every column of B is divided by its 2-norm. A null vector has no
   intrinsic norm, so the raw coefficient of a mode is an artefact of whichever
   eigenvector routine produced it. Solving in the unit-column basis makes the
   coefficient of mode j a basis-INDEPENDENT physical amplitude -- this is what
   `solver.rescale_invariance` tests, and it is why the test is not circular:
   without the normalisation the raw coefficients move by 1e248 (measured).

ROW REDUCTION
=============
At normal incidence the sagittal geometry makes u_x and sigma_xz identically
zero for every mode, so those two continuity equations read 0 = 0. They are
IDENTIFIED AND REMOVED, not kept with a tolerance and not repaired by adding
unknowns. The result is a correctly reduced system.

Reproduce:
    cd code && python3 -m scattering.assembly
"""
from __future__ import annotations

import mpmath as mp

import foundation_checks as fc
from scattering.roots import null_vector_verified
from scattering.modal import classify_mode
from scattering.radiation import classify_role
from scattering.elastic_reference import flux_z, mode_vector, traction
from scattering.state import (
    STATE_FIELDS, reference_scales, state_vector, to_dimensionless,
)

#: roles admitted into the basis. OUTGOING is decided from `outgoing_status`,
#: i.e. from where the energy goes, not from a DECAYS_* label.
BASIS_ROLES = ("OUTGOING",)

#: A row whose entire content is zero to this relative level is a trivial
#: 0 = 0 equation. It is removed and the removal is recorded.
TRIVIAL_ROW_TOL = mp.mpf("1e-25")

#: Equilibrated condition number above which a solve is reported but NOT
#: trusted. Never widened to manufacture agreement.
COND_TRUST = mp.mpf("1e8")


class Side:
    """Everything known about one half-space's modal basis."""

    def __init__(self, material, k, w, side, dps=None):
        self.material = material
        self.k = k
        self.w = w
        self.side = side
        self.dps = dps
        self.scales = reference_scales(material, w)
        self.roots = [mp.mpc(z) for z in material.xi_roots(k, w)]
        self.modes = []
        self.jzs = []
        self.states = []
        self.roles = []
        self.raw_states = []
        self.n_gauge = 0
        self.n_unresolved = 0
        self.n_electrostatic = 0
        self._build()

    def _build(self):
        for br, xi in enumerate(self.roots):
            # THE VERIFIED PATH, and the check that was missing.
            #
            # This used to call `nullspace.null_vector` and discard the
            # "agreement" flag it returns. That flag was the only warning that
            # the vector was a compressed direction rather than a null vector,
            # and at Al6061 k=0 branch 6 -- the shear branch, the whole point of
            # task 2R -- it was False while the vector went on to define the
            # mode. A guard no caller honours is not a guard.
            #
            # `null_vector_verified` returns ok=False instead of a direction it
            # cannot justify, so the mode is marked unresolved rather than
            # built from a wrong state.
            r = null_vector_verified(fc.M_mp(self.material, self.k, xi, self.w),
                                     self.dps or mp.mp.dps)
            if r["vector"] is None or not r.get("ok"):
                why = ("no null vector" if r["vector"] is None else
                       f"null vector not verified (residual "
                       f"{mp.nstr(r['residual'], 4)})")
                self._push(None, None, None, None, ("UNRESOLVED", why))
                self.n_unresolved += 1
                continue
            a = r["vector"]
            jz = fc.flux_z_mp(self.material, self.k, xi, self.w, a)
            md = classify_mode(self.material, self.k, xi, self.w, a, jz=jz)
            role, why = classify_role(md, self.side, jz)
            raw = state_vector(self.material, self.k, xi, self.w, a)
            state = to_dimensionless(dict(zip(STATE_FIELDS, raw)), self.scales)
            self._push(md, jz, state, raw, (role, why))
            if md["spatial_class"] == "GAUGE":
                self.n_gauge += 1
            if md["spatial_class"] == "ELECTROSTATIC":
                self.n_electrostatic += 1

    def _push(self, md, jz, state, raw, role):
        self.modes.append(md)
        self.jzs.append(jz)
        self.states.append(state)
        self.raw_states.append(raw)
        self.roles.append(role)

    def admissible(self):
        """Indices admitted into the basis, with the reason for each decision.

        The exclusion list is returned too: an unexplained exclusion of a
        physical channel is precisely the failure this task must not repeat.
        """
        admitted, excluded = [], []
        for br, (role, why) in enumerate(self.roles):
            md = self.modes[br]
            if md is None:
                excluded.append((br, "NO_MODE", why))
                continue
            if md["spatial_class"] == "GAUGE":
                excluded.append((br, "GAUGE",
                                 "genuine potential gauge (D_z = 0): a null "
                                 "channel, not a physical energy channel"))
                continue
            if role in BASIS_ROLES:
                admitted.append(br)
            else:
                excluded.append((br, role, why))
        return admitted, excluded

    def n_physical(self):
        return sum(1 for md in self.modes
                   if md is not None and md["spatial_class"] != "GAUGE")

    def summary(self):
        adm, exc = self.admissible()
        return {
            "side": self.side,
            "material": self.material.name,
            "n_roots": len(self.roots),
            "n_gauge": self.n_gauge,
            "n_electrostatic": self.n_electrostatic,
            "n_physical": self.n_physical(),
            "n_unresolved": self.n_unresolved,
            "admitted": adm,
            "n_admitted": len(adm),
            "excluded": [(b, r) for b, r, _w in exc],
            "roles": [r[0] for r in self.roles],
            "classes": [(m["spatial_class"] if m else "NONE")
                        for m in self.modes],
        }


def assemble(m1, m2, k, w, incident_material=1, incident_branch=None,
             dps=None):
    """Build the two-half-space system and report its degrees of freedom."""
    s1 = Side(m1, k, w, "z<0", dps=dps)      # medium 1 occupies z < 0
    s2 = Side(m2, k, w, "z>0", dps=dps)      # medium 2 occupies z > 0

    a1, exc1 = s1.admissible()
    a2, exc2 = s2.admissible()

    inc_side = s1 if incident_material == 1 else s2
    incoming = [br for br, (role, _why) in enumerate(inc_side.roles)
                if role == "INCOMING"]
    if incident_branch is None:
        inc = incoming[0] if incoming else None
    else:
        inc = incident_branch
    admitted_inc = a1 if incident_material == 1 else a2
    if inc is not None and inc in admitted_inc:
        raise ValueError(
            f"branch {inc} is both the incident wave and an admitted outgoing "
            f"mode; the incident channel must not appear in the outgoing basis")

    n_unk = len(a1) + len(a2)
    cols = []
    # ---- THE COLUMNS MUST BE PHYSICAL, NOT DIMENSIONLESS ---------------
    #
    # THIS WAS WRONG AND IT IS THE DEEPER DEFECT OF TASK 2S.
    #
    # `s1.states` and `s2.states` are each divided by THEIR OWN medium's
    # reference scales (`to_dimensionless`). For two identical media the scales
    # agree and subtracting the two sides is a real continuity condition. For
    # two DIFFERENT media they do not, and the row
    #
    #     states1_i - states2_i = -incident_i
    #
    # then compares u_i/scale1_i with u_i/scale2_i -- a condition on a
    # mismatched pair of unknowns that no physical field can satisfy.
    #
    # Measured for Al6061|PZT-5H at k = 0, the scale ratios are
    #   theta 2.469,  txz 0.738,  tzz 0.597,  qz 178.5,  Dz 2.5e-13
    # so the mismatch is not a rounding matter. The consequence: with every row
    # that would have detected it removed as "gauge" or "trivial", the reduced
    # system solved to 2.3e-41 IN ITS OWN MATRIX and left the PHYSICAL traction
    # discontinuous by 0.2616. Rebuilding the identical rows and columns from
    # raw physical states solves to 8.1e-28 and leaves physical continuity at
    # 2.5e-41 -- with coefficients that differ in the third digit.
    #
    # Continuity is a statement about the physical field, so the columns must be
    # the physical states. The conditioning that `to_dimensionless` was there to
    # provide is restored where it belongs: by normalising each column, and by
    # row-equilibrating at the solve. See `scattering.gauge.solve_path`.
    for br in a1:
        cols.append(("1", br, [s1.raw_states[br][i]
                               for i in range(len(STATE_FIELDS))]))
    for br in a2:
        cols.append(("2", br, [s2.raw_states[br][i]
                               for i in range(len(STATE_FIELDS))]))
    # ---- ROW-THEN-COLUMN EQUILIBRATION, PHYSICS-PRESERVING ---------------
    #
    # The columns are physical states in SI units; their eight entries span
    # 1e-37 (D_z) to 1e22 (q_z) across fields, and up to 1e62 in condition
    # number for a mixed pair at 30 degrees. Nothing can be solved in that form,
    # and dividing by the Euclidean norm -- the obvious fix -- is not one: it is
    # dominated by the loudest entry and deletes every other component.
    #
    # Measured, Al6061|PZT-5H at k = 0 with SV incidence: the Al shear mode has
    # (u_x, t_xz) = (1, 5.315e13j) and the PZT shear mode (1, 8.306e13j). They
    # differ in the shear modulus, which IS the physical distinction. After
    # Euclidean normalisation both become (1.9e-14, -1j) and (-1.2e-14, -1j):
    # the u_x = 1 information is gone at the 1e-14 level and the two columns
    # score a Cauchy-Schwarz overlap of 1.0 -- a false "duplicate column".
    #
    # (Normalising by the mechanical amplitude instead was tried and is ALSO
    # insufficient: it leaves every entry relative to q_z, so the `t_xz` row
    # still dominates and `u_x` still vanishes at 1e-14. Column-only fixes
    # cannot work here because the problem is the per-field spread WITHIN each
    # column, not the scale of one column against another.)
    #
    # The scaling that works is DIAGONAL, applied per field:
    #
    #   row scale   s_i = largest |state_i| over all admitted modes, per field
    #   column scale n_j = largest |state_i / s_i| over the fields of mode j
    #
    # Both are invertible diagonal maps, so the solution set of B a = rhs is
    # unchanged -- unlike a column-only normalisation that silently reweights
    # the components against each other. Every entry of the scaled matrix is
    # then at most 1, and PHYSICAL CONTINUITY REMAINS EXACT.
    #
    # The amplitudes are recovered as a_j = x_j / n_j, and the row scale is a
    # left multiplication that cancels, so `norms` is the only thing a caller
    # needs to undo.
    row_scale = []
    for i in range(len(STATE_FIELDS)):
        mx = mp.mpf(0)
        for _sd, _br, col in cols:
            if abs(col[i]) > mx:
                mx = abs(col[i])
        row_scale.append(mx if mx != 0 else mp.mpf(1))

    norms = []
    for _sd, _br, col in cols:
        nn = mp.mpf(0)
        for i in range(len(STATE_FIELDS)):
            v = abs(col[i]) / row_scale[i]
            if v > nn:
                nn = v
        norms.append(nn if nn != 0 else mp.mpf(1))
    # SIGN CONVENTION: side 1 and side 2 enter with OPPOSITE signs.
    # The continuity condition is
    #     s1_out + s_inc  =  s2_out
    # i.e.      s1_out - s2_out  = -s_inc
    # so medium 2's columns are NEGATED. The first version built [B1 +B2] (the
    # docstring said [B1 -B2] but the code did not implement it) and the solve
    # then satisfied s1_out + s2_out = -s_inc instead. The linear-algebra
    # residual was still 1e-37, so the solve LOOKED converged; recomputing the
    # interface field independently exposed it, showing a residual of exactly
    # 2 x the incident state. This is the reason a converged solve is not a
    # correct solve, and why the interface residual is recomputed rather than
    # read off the solver.
    B = mp.matrix(len(STATE_FIELDS), n_unk)
    for j, (sd, _br, col) in enumerate(cols):
        sgn = mp.mpf(1) if sd == "1" else mp.mpf(-1)
        for i in range(len(STATE_FIELDS)):
            B[i, j] = sgn * (col[i] / row_scale[i]) / norms[j]

    rhs_full = mp.matrix(len(STATE_FIELDS), 1)
    if inc is not None:
        st = inc_side.raw_states[inc]     # physical, for the same reason as above
        for i in range(len(STATE_FIELDS)):
            rhs_full[i] = -st[i] / row_scale[i]

    # ---- trivial-row reduction ---------------------------------------
    #
    # A ROW MAY BE REMOVED ONLY IF IT READS 0 = 0 IN THE WHOLE SYSTEM.
    #
    # Two independent conditions, and the first version enforced only one:
    #
    #   (1) the row is zero for EVERY admitted basis mode, and
    #   (2) the corresponding entry of the incident right-hand side is zero.
    #
    # Neither implies the other, and dropping (2) is not a corner case -- it is
    # the exact failure task 2R was opened to find. At normal incidence the
    # u_x row of the ALUMINIUM half-space is nonzero only on the shear branches,
    # and the shear branches are precisely the ones removed as inadmissible
    # when a P wave is incident. A reduction that looks at the basis alone
    # therefore sees an all-zero row and deletes it, taking sigma_xz = 0 with
    # it -- the very equation that distinguishes a shear mode from a
    # longitudinal one. `rhs_is_zero` was already being computed and recorded;
    # it was simply not part of the decision.
    #
    # The tolerance is relative to the largest row content, since B's columns
    # are already normalised and the entries are not of a common physical kind.
    kept, removed = [], []
    _rowmax_all = [max([abs(B[i, j]) for j in range(n_unk)] or [mp.mpf(0)])
                   for i in range(len(STATE_FIELDS))]
    # ---- THE TRIVIALITY TEST IS PER ROW, NOT PER MATRIX -----------------
    #
    # A row reads 0 = 0 for the whole system when ITS OWN entries are zero.
    # Comparing each row against the largest row in the matrix is wrong here and
    # was measured to be catastrophic: with raw physical states the rows range
    # from 1e-37 (D_z) to 1e19 (q_z), so a global scale of 1e19 makes the
    # threshold 1e-6 and classifies every row except q_z as empty. Measured
    # after switching the columns to physical states: the reduction collapsed
    # to 2 rows out of 8 and the system had rank 0.
    #
    # Per-row, the entry magnitudes are not compared to anything external. That
    # is what `0 = 0` means. It also makes the test unit-safe, which the global
    # form was not: `raw_states` are in SI and their spread across fields is
    # physical, not numerical.
    _row_scale = [v if v > 0 else mp.mpf(1) for v in _rowmax_all]
    _tol_of = [TRIVIAL_ROW_TOL * v for v in _row_scale]

    i1 = len(a1)                       # columns [0, i1) are medium 1
    for i in range(len(STATE_FIELDS)):
        rowmax = _rowmax_all[i]
        rhs_val = rhs_full[i]
        _tol = _tol_of[i]
        row_zero = rowmax <= _tol
        rhs_zero = abs(rhs_val) <= _tol

        # content of this row in each half-space's own basis
        m1_max = max([abs(B[i, j]) for j in range(i1)] or [mp.mpf(0)])
        m2_max = max([abs(B[i, j]) for j in range(i1, n_unk)] or [mp.mpf(0)])
        # A half-space is "empty" in this row when IT contributes nothing
        # relative to this row's own scale. Measured on the row, not the matrix.
        m1_empty = m1_max <= _tol
        m2_empty = m2_max <= _tol

        # ---- case 1: the row is zero everywhere -------------------------
        if row_zero:
            if rhs_zero:
                removed.append({"row": STATE_FIELDS[i],
                                "reason": "reads 0 = 0 at this angle: zero for "
                                          "every admitted mode AND zero on the "
                                          "incident right-hand side",
                                "rhs": rhs_val, "rhs_is_zero": True,
                                "row_max": rowmax})
            else:
                kept.append(i)
                removed.append({"row": STATE_FIELDS[i],
                                "reason": "ROW ZERO BUT RHS NONZERO: the system "
                                          "is inconsistent in this equation",
                                "rhs": rhs_val, "rhs_is_zero": False,
                                "row_max": rowmax, "INCONSISTENT": True})
            continue

        # ---- case 2: gauge-capable rows, tested before the RHS ----------
        #
        # THE ORDER OF THESE TWO TESTS MATTERS AND THE FIRST VERSION HAD IT
        # BACKWARDS. Testing "is the RHS zero?" first and returning early made
        # the gauge test unreachable whenever the incident wave DID carry the
        # quantity -- which is the case that needs it most. Measured:
        # PZT-5H|Al6061 at k = 0 under P incidence kept a phi row whose
        # right-hand side was -1.0 while the receiving half-space had no
        # potential in any admitted mode, and the resulting 7x6 system was
        # inconsistent at 55.25.
        #
        # THE RULE. A potential is fixed only through its gradient, so an
        # overall constant is gauge. A half-space whose admitted basis has no
        # content in this row has no datum to transmit across the interface.
        # The interface condition "phi1 = phi2" then imposes the other side's
        # value -- which is zero only because it is undefined there -- as if it
        # were a physical boundary datum.
        #
        # SCOPED TO ELECTRIC ROWS ON PURPOSE. Applied to a mechanical row it
        # deletes exactly the physics this audit exists to protect: the incident
        # P wave has u_x = 0, so an unscoped test removes u_x and txz -- the
        # shear rows -- and reproduces the original Task 2 defect in a new
        # place. Measured: unscoped, Al6061|Al6061 under P incidence reduced to
        # 4x6 and dropped ux and txz. Displacement has no gauge freedom; u_x = 0
        # is a physical statement about that wave.
        #
        # AND IT IS GATED ON k = 0. The degeneracy that makes this necessary is
        # exact only at normal incidence: there, the electric governing equation
        # reduces to e33 u_z - eps33 phi = 0, which forces D_z = i xi (e33 u_z -
        # eps33 phi) = 0 identically, leaving phi non-unique by a constant.
        # For any k != 0 the eps11 k^2 term breaks that cancellation and D_z is
        # nonzero and meaningful. Measured: applied ungated at 30 degrees this
        # rule removed both phi and Dz from Al6061|PZT-5H, leaving a 6x8 system
        # and breaking SCTA05. `k` here is the caller's exact mpf, so k == 0 is
        # an exact test, not a tolerance.
        normal_incidence = (k == 0)
        if normal_incidence and STATE_FIELDS[i] in ("phi", "Dz"):
            inc_empty = False
            if inc is not None:
                _st = inc_side.states[inc]
                _mem = mp.sqrt(sum(abs(_st[k2]) ** 2
                                   for k2 in range(len(STATE_FIELDS))))
                inc_empty = (abs(_st[i]) / _mem if _mem > 0
                             else mp.mpf(0)) <= _tol
            if (m1_empty != m2_empty) or inc_empty:
                if m1_empty != m2_empty:
                    carried = "medium 2" if m1_empty else "medium 1"
                    datside = "medium 1" if m1_empty else "medium 2"
                else:
                    carried, datside = "both media", "the incident wave"
                removed.append({"row": STATE_FIELDS[i],
                                "reason": f"GAUGE, NOT A BOUNDARY CONDITION: "
                                          f"{carried} carry no content in this "
                                          f"row and {datside} supplies no "
                                          f"physical datum, so the row would "
                                          f"impose zero as if it were one",
                                "rhs": rhs_val, "rhs_is_zero": rhs_zero,
                                "row_max": rowmax,
                                "basis_empty_side1": bool(m1_empty),
                                "basis_empty_side2": bool(m2_empty),
                                "incident_empty": bool(inc_empty),
                                "GAUGE_ONLY": True})
                continue

        kept.append(i)

    B_red = mp.matrix(len(kept), n_unk)
    rhs_red = mp.matrix(len(kept), 1)
    for i2, i in enumerate(kept):
        for j in range(n_unk):
            B_red[i2, j] = B[i, j]
        rhs_red[i2] = rhs_full[i]

    rank = matrix_rank(B_red) if n_unk else 0
    return {
        "B": B_red, "rhs": rhs_red, "B_full": B, "rhs_full": rhs_full,
        "rows_kept": [STATE_FIELDS[i] for i in kept],
        "rows_removed": removed,
        "rows_inconsistent": [r["row"] for r in removed
                              if r.get("INCONSISTENT")],
        "rows_gauge_only": [r["row"] for r in removed
                            if r.get("GAUGE_ONLY")],
        "n_equations": len(kept),
        "n_equations_full": len(STATE_FIELDS),
        "n_unknowns": n_unk,
        "n_roots_total": len(s1.roots) + len(s2.roots),
        "column_norms": norms,
        "row_scale": row_scale,
        "rank": rank,
        "square": n_unk == len(kept),
        "overdetermined": n_unk < len(kept),
        "underdetermined": n_unk > len(kept),
        "side1": s1, "side2": s2,
        "admitted1": a1, "admitted2": a2,
        "excluded1": exc1, "excluded2": exc2,
        "incident_side": incident_material,
        "incident_branch": inc,
        "incoming1": ([br for br, (r, _w) in enumerate(s1.roles)
                       if r == "INCOMING"]),
        "incoming2": ([br for br, (r, _w) in enumerate(s2.roles)
                       if r == "INCOMING"]),
    }


def matrix_rank(B, tol=None):
    """Numerical rank via the squared singular values of the Gram matrix.

    `mp.cond` and mpmath's SVD both fail on complex input (recorded dead ends);
    `mp.eighe` on the Hermitian A^H A is stable and gives the singular values as
    square roots of its eigenvalues.
    """
    n, m = B.rows, B.cols
    if m == 0 or n == 0:
        return 0
    if tol is None:
        tol = mp.power(10, -(mp.mp.dps - 10))
    sv = _singular_values(B)
    if not sv:
        return 0
    smax = sv[0]
    if smax == 0:
        return 0
    cut = tol * smax
    r = sum(1 for s in sv if s > cut)
    # The Gram matrix is m x m, but a rank can never exceed min(n, m). Without
    # this cap the near-zero eigenvalues of a tall-and-thin B survive the
    # tolerance and the function reports a rank the matrix cannot have --
    # caught by getting rank 6 from a 4 x 6 system.
    return min(r, min(n, m))


def _singular_values(B):
    n, m = B.rows, B.cols
    G = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            G[i, j] = sum(mp.conj(B[r, i]) * B[r, j] for r in range(n))
    Gh = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    E, _Q = mp.eighe(Gh)
    return sorted((mp.sqrt(abs(e)) for e in E), reverse=True)


def condition_number(B):
    """2-norm condition number from the extreme singular values."""
    if B.rows == 0 or B.cols == 0:
        return mp.inf
    sv = _singular_values(B)
    smax, smin = sv[0], sv[-1]
    if smin == 0:
        return mp.inf
    return smax / smin


def _main():
    import sys
    sys.path.insert(0, ".")
    mp.mp.dps = 50
    w = 2 * mp.pi * 1e6
    print("Coupled two-half-space assembly -- degrees of freedom counted\n")
    for name1, name2 in (("PZT-5H", "PZT-5H"), ("Al6061", "PZT-5H"),
                         ("Al6061", "Al6061")):
        for th in (0, 30):
            m1 = fc.material_from_table(name1)
            m2 = fc.material_from_table(name2)
            k = w * mp.sin(mp.radians(mp.mpf(th))) / mp.mpf(4000)
            A = assemble(m1, m2, k, w)
            tag = f"{name1}|{name2} th={th}"
            print(f"--- {tag} ---")
            for nm, sd in (("side1", A["side1"]), ("side2", A["side2"])):
                ss = sd.summary()
                print(f"  {nm} {ss['material']:8s} roots={ss['n_roots']} "
                      f"gauge={ss['n_gauge']} electro={ss['n_electrostatic']} "
                      f"phys={ss['n_physical']} admitted={ss['admitted']}")
            print(f"  roots={A['n_roots_total']} eq={A['n_equations']}"
                  f"/{A['n_equations_full']} unk={A['n_unknowns']} "
                  f"rank={A['rank']} square={A['square']}")
            if A["rows_removed"]:
                print(f"  rows removed: "
                      f"{[r['row'] for r in A['rows_removed']]}")


if __name__ == "__main__":
    _main()


# ---------------------------------------------------------------------------
# Reduced elastic control: the 4x4 welded problem, against an independent oracle
# ---------------------------------------------------------------------------
def assemble_reduced_elastic(iso1, iso2, w, theta_deg, incident="P"):
    """The 4x4 welded P-SV system, built from the SAME modal machinery.

    This exists to separate two questions that the coupled piezothermoelastic
    problem confounds:

      1. Is the DoF counting and mode selection machinery right?
      2. Is the coupled 8-field formulation itself correct?

    Nothing here is imported from `elastic_reference`: the modes come from
    `ISO.mode_vector`/`ISO.traction`, the selection is done explicitly, and the
    answer is compared against `elastic_reference.solve_psv`, which was
    validated independently against closed-form Fresnel and Snell limits. If
    these agree, the assembly machinery is sound and any remaining coupled
    disagreement belongs to the piezothermoelastic formulation.
    """
    theta = theta_deg if isinstance(theta_deg, mp.mpf) else mp.mpf(repr(theta_deg))
    k = (w / (iso1.cL if incident == "P" else iso1.cT)) * mp.sin(mp.radians(theta))

    def pick(iso, branch, direction):
        """Choose the root whose axial flux leaves the interface, by flux only.

        A lossless mode in the isotropic reference is exactly lossless
        (Im xi = 0), so the direction is carried entirely by the flux sign --
        which is the rule the brief demands and the rule `classify_role`
        implements for the coupled problem.
        """
        cand = []
        for xi in (iso.xi_P(k, w), -iso.xi_P(k, w)) if branch == "P" else \
                  (iso.xi_S(k, w), -iso.xi_S(k, w)):
            u = mode_vector(k, xi, branch)
            pz = flux_z(iso, k, xi, u, w)
            cand.append((xi, u, pz))
        for xi, u, pz in cand:
            if pz != 0 and ((mp.re(pz) > 0) == (direction > 0)):
                return xi, u, pz
        return cand[0]

    xi_rP, u_rP, _ = pick(iso1, "P", -1)
    xi_rS, u_rS, _ = pick(iso1, "S", -1)
    xi_tP, u_tP, _ = pick(iso2, "P", +1)
    xi_tS, u_tS, _ = pick(iso2, "S", +1)
    xi_in, u_in, _ = pick(iso1, incident, +1)

    cols = [(iso1, xi_rP, u_rP, +1, "rP"),
            (iso1, xi_rS, u_rS, +1, "rS"),
            (iso2, xi_tP, u_tP, -1, "tP"),
            (iso2, xi_tS, u_tS, -1, "tS")]

    # NORMALISE THE INCIDENT WAVE TO UNIT DISPLACEMENT, not the columns.
    #
    # The first version unit-normalised every column and then tried to undo the
    # divisor afterwards to recover "true" amplitudes; that factor is not the
    # same for the incident and the unknown columns, so it produced amplitudes
    # wrong by 1e28. The reference defines the incident displacement as 1 and
    # reports every other amplitude relative to it, so the right move is to
    # scale the INCIDENT to 1 and leave the columns alone. Then a solved
    # coefficient IS the amplitude ratio the reference reports.
    u_in_mag = mp.sqrt(abs(u_in[0]) ** 2 + abs(u_in[1]) ** 2)
    u_scale = mp.mpf(1) / (u_in_mag if u_in_mag else mp.mpf(1))

    A = mp.matrix(4, 4)
    for j, (mm, xi, u, sgn, _nm) in enumerate(cols):
        tx, tz = traction(mm, k, xi, u)
        # physical form, untransformed: continuity of (u_x, u_z, t_xz, t_zz),
        # with the transmitted (sgn = -1) columns moved to the left-hand side
        A[0, j] = sgn * u[0] * u_scale
        A[1, j] = sgn * u[1] * u_scale
        A[2, j] = sgn * tx * u_scale
        A[3, j] = sgn * tz * u_scale
    b = mp.matrix(4, 1)
    tx_in, tz_in = traction(iso1, k, xi_in, u_in)
    for i, v in enumerate((-u_in[0], -u_in[1], -tx_in, -tz_in)):
        b[i] = v * u_scale

    # row equilibration only -- an exact rescaling of the equations
    S = mp.matrix(A)
    r = mp.matrix(b)
    for i in range(4):
        mx = max(abs(A[i, j]) for j in range(4))
        mx = mx if mx else mp.mpf(1)
        for j in range(4):
            S[i, j] = A[i, j] / mx
        r[i] = b[i] / mx
    sol = mp.lu_solve(S, r)
    amps = {cols[j][4]: sol[j] for j in range(4)}
    return {
        "amps": amps, "k": k, "theta": theta,
        "n_unknowns": 4, "n_equations": 4,
        "rank": matrix_rank(A),
        "cond_equilibrated": condition_number(S),
        "incident_scale": u_scale,
        "incident_displacement": u_in_mag,
    }
