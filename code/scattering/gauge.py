"""
code/scattering/gauge.py -- electrical gauge versus electrical boundary condition.

WHY THIS EXISTS
===============
Task 2R removed the `phi` and `Dz` interface rows at normal incidence and
recorded the reason as "the basis has no content in this row, so the incident
zero would be imposed as a physical datum". That reasoning reaches a defensible
conclusion by an argument that is not the real one.

THE REAL ARGUMENT, DERIVED
==========================
Under exp(i(kx + xi z)) with plane strain, from the 6mm constitutive relations:

    D_z = e31 e_xx + e33 e_zz + eps33 E_z + p theta
    E_z = -i xi phi,   e_zz = i xi u_z,   e_xx = i k u_x

Gauss's law with no free charge is

    div D = i k D_x + i xi D_z = 0

At k = 0 the first term carries its own factor of k and vanishes, leaving

    i xi D_z = 0   ->   D_z = 0    for xi != 0.

Combining the two expressions for D_z gives, at k = 0,

    -xi^2 (e33 u_z - eps33 phi) + i p xi theta = 0
    =>  phi = e33 u_z / eps33 - i p theta / (eps33 xi)

and substituting back:

    D_z = i xi (e33 u_z - eps33 phi) + p theta = -p theta + p theta = 0

So D_z = 0 is an IDENTITY on shell at k = 0, not a constraint. It follows from
Gauss's law after phi has been used to eliminate it. Two consequences:

  1. The `Dz` row of the interface system reads 0 = 0 for every mode and for
     every incident datum. Removing it removes no condition. [correct]

  2. `phi` is NOT a free field at k = 0. It is SLAVED to (u_z, theta) by the
     relation above. It is not an independent unknown.

WHY THAT MAKES THE `phi` ROW REDUNDANT RATHER THAN MERELY EMPTY
---------------------------------------------------------------
Because phi is slaved to (u_z, theta), the statement "phi_1 = phi_2" is a
statement about (u_z, theta) across the interface. The mechanical and thermal
rows already impose that. So the phi row is linearly dependent on rows that are
kept -- which is measured directly here as a rank that does not increase when
the row is added.

THE ADDITIVE CONSTANT, AND WHAT IT IS FOR
-----------------------------------------
The dispersion spectrum contains the xi = 0 root, whose state is (0,0,1,0,...):
a constant potential with no displacement, no traction, no flux. It is the
plane-wave representative of the additive gauge phi -> phi + c. Each half-space
has one, because phi is continuous only up to a constant in each.

For a bonded interface the physical statement is phi_1(0) = phi_2(0). Writing
phi_i = phi_i^modes + C_i, that is one equation relating C_1 and C_2. It fixes
their DIFFERENCE. Adding the same constant to both leaves every physical field
invariant, so exactly ONE global gauge freedom survives.

A half-space with no piezoelectric coupling (e31 = e33 = e15 = 0) has phi = 0
for every xi != 0 mode, so its entire potential is the constant. The potential
continuity condition then constrains only that constant and carries no
mechanical information. THAT is why dropping the row loses nothing for an
Al6061|PZT-5H pair -- not because a zero "would be imposed as a datum".

WHAT WOULD BREAK IT
-------------------
For two DIFFERENT piezoelectric materials at k = 0, phi is slaved to
(u_z, theta) with material-dependent coefficients, so phi_1(0) and phi_2(0)
differ even when u_z and theta are continuous. There the row is NOT redundant
and carries real information. `phi_row_rank_gain` measures exactly this; a
nonzero value means the current reduction would be wrong for that pair. See
BLK-02A-19.

WHAT THIS MODULE PROVIDES
-------------------------
`full_system`   -- every physical mode as an unknown, the xi = 0 gauge modes
                   kept as EXPLICIT gauge variables, all eight interface rows
                   imposed, and one global gauge condition.
`reduced_system`-- what the production assembly builds.
`compare`       -- physical fields from both, compared component by component.

The comparison is on PHYSICAL fields, never on the coefficient vector: a
successful reduced solve is not a proof of equivalence, and this module is
written so that it cannot be mistaken for one.
"""
from __future__ import annotations

import mpmath as mp

import foundation_checks as fc
from scattering.assembly import Side, assemble
from scattering.state import STATE_FIELDS

#: The two rows the electrical-branch question is about.
ELECTRIC_ROWS = ("phi", "Dz")


def characteristic_rows():
    """Rows of the interface system in STATE_FIELDS order."""
    return STATE_FIELDS


# ---------------------------------------------------------------------------
# the full formulation
# ---------------------------------------------------------------------------
def full_system(m1, m2, k, w, incident_branch, dps=None):
    """All eight interface rows, with the xi = 0 sector kept as gauge variables.

    Unknowns, in this order:
        a_1j   physical modes of medium 1 (xi != 0, admissible)
        a_2j   physical modes of medium 2
        C_1    the constant potential of medium 1
        C_2    the constant potential of medium 2

    Rows, in STATE_FIELDS order, all eight of them, plus ONE global gauge
    condition C_1 = 0. The gauge condition is needed because (C_1, C_2) ->
    (C_1 + c, C_2 + c) changes nothing physical; any single fixing removes it.

    Everything is built from RAW physical states -- no dimensionless scaling --
    so that the fields this returns are directly comparable with the reduced
    solve. Ill-conditioning from the mixed units is handled at the solve, not by
    changing the equations.
    """
    if dps is not None:
        mp.mp.dps = dps

    # ---- THE FULL SYSTEM IS NOT DEFINED AT EXACTLY k = 0 ----------------
    #
    # Polynomial mode vectors are amplitudes multiplying exp(i k x + i xi z).
    # As k -> 0 the potential entry of the PZT branches grows like 1/k, measured
    # at 2.63e11 (k/k_ac = 1e-2), 2.63e13 (1e-4), 2.63e15 (1e-6), 2.63e17
    # (1e-8) -- a clean 1/k. There is therefore no k = 0 value of the potential
    # to impose continuity on; the k = 0 answer is a LIMIT.
    #
    # MEASURED CONSEQUENCE, and why this guard exists: assembling the full
    # system at exactly k = 0 anyway (the additive constants C_1, C_2 absorb the
    # divergent potential, so the system still "solves") gave a transmitted
    # t_zz / incident t_zz of 1.3316249147 + 0.00115753873j, with C_2 = 1.196e9.
    # The same ratio at k/k_ac = 1e-4, 1e-5 is 1.3317978949 and 1.33179789532,
    # and the REDUCED solve at k = 0 gives 1.33179789533. The k = 0 full answer
    # is wrong in the fourth digit and invents an imaginary part; the reduced
    # answer is the limit. The full formulation is only valid for k > 0.
    if k == 0:
        raise ValueError(
            "full_system requires k != 0: the potential of a polynomial mode "
            "diverges like 1/k as k -> 0, so potential continuity is a "
            "statement about the LIMIT, not about k = 0. Use reduced_system at "
            "k = 0, or full_system at small k and take the limit.")

    s1 = Side(m1, k, w, "z<0", dps=dps)
    s2 = Side(m2, k, w, "z>0", dps=dps)

    def physical_unknowns(side):
        """Modes admissible as outgoing on this side, excluding the gauge sector."""
        out = []
        for br in range(8):
            if side.raw_states[br] is None:
                continue
            if side.roots[br] == 0:
                continue                       # xi = 0: handled as a gauge variable
            role = side.roles[br][0]
            if role in ("INCOMING", "EVANESCENT_NO_TRANSPORT",
                        "GAUGE_EXCLUDED", "UNRESOLVED"):
                continue
            out.append(br)
        return out

    def gauge_branches(side):
        return [br for br in range(8)
                if side.raw_states[br] is not None and side.roots[br] == 0]

    a1 = physical_unknowns(s1)
    a2 = physical_unknowns(s2)
    g1 = gauge_branches(s1)
    g2 = gauge_branches(s2)

    # a gauge variable is the constant potential; it exists iff the xi = 0
    # branch is realisable in that medium
    has_c1 = bool(g1)
    has_c2 = bool(g2)

    cols = []
    for br in a1:
        cols.append(("1", br, [s1.raw_states[br][i] for i in range(8)]))
    for br in a2:
        cols.append(("2", br, [s2.raw_states[br][i] for i in range(8)]))
    if has_c1:
        st = [mp.mpc(0)] * 8
        st[2] = mp.mpc(1)
        cols.append(("1", "C", st))
    if has_c2:
        st = [mp.mpc(0)] * 8
        st[2] = mp.mpc(1)
        cols.append(("2", "C", st))

    n = len(cols)
    n_rows = len(STATE_FIELDS) + 1          # eight conditions + one gauge fixing
    B = mp.matrix(n_rows, n)
    rhs = mp.matrix(n_rows, 1)

    for j, (sd, br, col) in enumerate(cols):
        sgn = mp.mpf(1) if sd == "1" else mp.mpf(-1)
        for i in range(8):
            B[i, j] = sgn * col[i]

    st_inc = s1.raw_states[incident_branch]
    for i in range(8):
        rhs[i] = -st_inc[i]

    # ---- one global gauge condition -----------------------------------
    #
    # C_1 = 0. Any single condition that breaks (C_1, C_2) -> (C_1+c, C_2+c)
    # will do; fixing the first constant is the least arbitrary choice because
    # it leaves medium 2's constant free to be DETERMINED by potential
    # continuity, which is the physically meaningful combination.
    gauge_row = len(STATE_FIELDS)
    if has_c1:
        B[gauge_row, n - (2 if has_c2 else 1)] = mp.mpf(1)
        rhs[gauge_row] = mp.mpf(0)
    elif has_c2:
        B[gauge_row, n - 1] = mp.mpf(1)
        rhs[gauge_row] = mp.mpf(0)

    return {
        "B": B, "rhs": rhs, "columns": cols,
        "side1": s1, "side2": s2,
        "a1": a1, "a2": a2, "gauge1": g1, "gauge2": g2,
        "has_c1": has_c1, "has_c2": has_c2,
        "incident_branch": incident_branch,
        "n_unknowns": n, "n_rows": n_rows,
        "n_gauge_variables": (1 if has_c1 else 0) + (1 if has_c2 else 0),
        "n_global_gauge_conditions": 1,
    }


def physical_fields_full(sol, system):
    """The eight-component fields on each side, from the FULL solve."""
    coeffs = sol
    s1, s2 = system["side1"], system["side2"]
    n1 = len(system["a1"])
    n2 = len(system["a2"])
    f1 = [mp.mpc(0)] * 8
    f2 = [mp.mpc(0)] * 8
    for j, br in enumerate(system["a1"]):
        for i in range(8):
            f1[i] += coeffs[j] * s1.raw_states[br][i]
    for j, br in enumerate(system["a2"]):
        for i in range(8):
            f2[i] += coeffs[n1 + j] * s2.raw_states[br][i]
    c1 = c2 = mp.mpc(0)
    idx = n1 + n2
    if system["has_c1"]:
        c1 = coeffs[idx]
        idx += 1
    if system["has_c2"]:
        c2 = coeffs[idx]
    f1[2] += c1
    f2[2] += c2
    return f1, f2, c1, c2


# ---------------------------------------------------------------------------
# the reduced formulation
# ---------------------------------------------------------------------------
def reduced_system(m1, m2, k, w, incident_branch, dps=None):
    """What the production assembly builds, with its own bookkeeping."""
    if dps is not None:
        mp.mp.dps = dps
    return assemble(m1, m2, k, w, incident_material=1,
                    incident_branch=incident_branch, dps=dps)


def physical_fields_reduced(sol, system):
    """The eight-component fields on each side, from the REDUCED solve.

    The columns are unit-normalised, so the physical contribution of column j is
    coeff_j / ||col_j||. Getting that divisor wrong understates the field by the
    column norm -- a mistake made once already in `solver.boundary_residuals`.

    At exact k = 0 the potential row is not a scattering equation (see
    assembly.assemble case 2): the potential is reconstructed here with the
    relative gauge fixed by [[phi]] = 0, leaving one recorded global gauge.
    """
    s1, s2 = system["side1"], system["side2"]
    a1, a2 = system["admitted1"], system["admitted2"]
    norms = system["column_norms"]
    f1 = [mp.mpc(0)] * 8
    f2 = [mp.mpc(0)] * 8
    for j, br in enumerate(a1):
        c = sol[j] / norms[j]
        for i in range(8):
            f1[i] += c * s1.raw_states[br][i]
    for j, br in enumerate(a2):
        c = sol[len(a1) + j] / norms[len(a1) + j]
        for i in range(8):
            f2[i] += c * s2.raw_states[br][i]
    if system.get("k") == 0 and "phi" in system.get("potential_rows", []):
        # [[phi]] = 0 is a statement about the TOTAL fields, and on the
        # incident side the total includes the incident wave's potential.
        inc_phi = mp.mpc(0)
        if system.get("incident_side", 1) == 1 and system.get("incident_branch") is not None:
            inc_phi = s1.raw_states[system["incident_branch"]][2]
        f1t = list(f1)
        f1t[2] = f1[2] + inc_phi
        f1, f2, c1, c2 = reconstruct_potential(f1t, f2)
        f1[2] = f1[2] - inc_phi          # report the scattered part, as before
        system["gauge_constants"] = (c1, c2)
        system["incident_phi"] = inc_phi
    return f1, f2


def reconstruct_potential(f1, f2, global_gauge="side1_zero"):
    """Impose [[phi]] = 0 on the k = 0 modal fields, modulo one global gauge.

    At k = 0 each half-space's potential is phi_i = phi_i^modes + C_i. The
    bonded-contact condition phi_1(0) = phi_2(0) fixes only the RELATIVE
    offset C_2 - C_1 = phi_1^modes - phi_2^modes; adding the same constant to
    both leaves every physical field invariant, so exactly one global gauge
    freedom survives. Here it is fixed by C_1 = 0 (`global_gauge="side1_zero"`),
    and the choice is returned, not hidden.

    Returns (f1, f2, C1, C2) with the constants added to the potential entry.
    """
    phi = STATE_FIELDS.index("phi")
    c1 = mp.mpc(0) if global_gauge == "side1_zero" else None
    if c1 is None:
        raise ValueError(f"unknown global gauge fixing {global_gauge!r}")
    c2 = f1[phi] - f2[phi]
    f1 = list(f1)
    f2 = list(f2)
    f1[phi] = f1[phi] + c1
    f2[phi] = f2[phi] + c2
    return f1, f2, c1, c2


# ---------------------------------------------------------------------------
# solving
# ---------------------------------------------------------------------------
def solve_path(B, rhs, dps=None):
    """Least-squares solve of B a = rhs, in ROW-EQUILIBRATED form.

    ROW EQUILIBRATION IS MANDATORY HERE, NOT A REFINEMENT. The rows of B are
    equations in different physical units and their entries span 1e-37 (D_z) to
    1e17 (q_z). Measured: the raw normal matrix has eigenvalues from 5.6e-5 to
    9.2e36, a condition number of 1.6e41 at 40 digits, and mp.lu_solve fails
    outright. Dividing each row by its own largest entry leaves the solution set
    unchanged -- it is a diagonal left multiplication by an invertible matrix --
    and brings the condition number down to something 40-digit arithmetic can
    handle.

    Returns (coeffs, residual_max, cond_estimate). The residual is measured on
    the ORIGINAL system so it is a statement about the equations, not about the
    scaling used to solve them.
    """
    if dps is not None:
        mp.mp.dps = dps
    n, m = B.rows, B.cols
    S = mp.matrix(n, m)
    t = mp.matrix(n, 1)
    scale = []
    for i in range(n):
        mx = max(abs(B[i, j]) for j in range(m))
        d = mx if mx != 0 else mp.mpf(1)
        scale.append(d)
        for j in range(m):
            S[i, j] = B[i, j] / d
        t[i] = rhs[i] / d

    # ---- COLUMN EQUILIBRATION, then the normal equations ---------------
    #
    # Row equilibration alone is not enough. `full_system` builds columns from
    # raw states at k > 0, where the potential entry grows like 1/k and reaches
    # 1e14 at k/k_ac = 1e-2 and 1e20 at 1e-8. The normal matrix S^H S is then
    # singular at any working precision and `mp.lu_solve` raises. Column scaling
    # is a diagonal RIGHT multiplication, so it changes the unknowns
    # (a_j = x_j / c_j) and not the solution set -- the caller is returned the
    # true amplitudes, already unscaled.
    cscale = []
    for j in range(m):
        mx = max(abs(S[i, j]) for i in range(n)) if n else mp.mpf(0)
        cscale.append(mx if mx != 0 else mp.mpf(1))
    S2 = mp.matrix(n, m)
    for i in range(n):
        for j in range(m):
            S2[i, j] = S[i, j] / cscale[j]

    G = S2.H * S2
    Gh = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    ev, _ = mp.eighe(Gh)
    mags = sorted(abs(e) for e in ev)
    cond = mags[-1] / mags[0] if mags[0] != 0 else mp.inf
    x = mp.lu_solve(Gh, S2.H * t)
    a = mp.matrix(m, 1)
    for j in range(m):
        a[j] = x[j] / cscale[j]
    resid = max(abs(v) for v in (B * a - rhs)) if n else mp.mpf(0)
    return a, resid, cond


# ---------------------------------------------------------------------------
# rank bookkeeping
# ---------------------------------------------------------------------------
def rank_at(A, tol_exp=-20):
    """Numerical rank by the singular values of the row-equilibrated matrix.

    Reported alongside the tolerance scaling used, because a rank is only
    meaningful together with the tolerance it was computed at. mp.svd and
    mp.cond both fail on complex input, so mp.eighe on A^H A is used.
    """
    n, m = A.rows, A.cols
    Ack = mp.matrix(n, m)
    for i in range(n):
        mx = max(abs(A[i, j]) for j in range(m))
        d = mx if mx != 0 else mp.mpf(1)
        for j in range(m):
            Ack[i, j] = A[i, j] / d
    G = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            G[i, j] = sum(mp.conj(Ack[r, i]) * Ack[r, j] for r in range(n))
    Gh = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    ev, _ = mp.eighe(Gh)
    sv = sorted((mp.sqrt(abs(e)) for e in ev), reverse=True)
    smax = sv[0] if sv else mp.mpf(0)
    if smax == 0:
        return 0, sv
    tol = mp.power(10, tol_exp) * smax
    return sum(1 for s in sv if s > tol), sv


def phi_row_rank_gain(full):
    """Does the phi row raise the rank BEYOND the modes alone?

    TWO RANKS, AND THE DISTINCTION IS THE WHOLE POINT.

    On the FULL system the phi row trivially adds rank, because it is the only
    row in which the gauge constants C_1 and C_2 appear. That measurement says
    nothing about whether potential continuity carries physical information; it
    only says the constants need somewhere to live.

    The question that matters is the one asked of the MODE-ONLY subsystem -- the
    system you get after the constants are eliminated, which is the form both
    the production reduction and the "gauge-invariant reduced system" use. There
    the derivation predicts zero gain: phi is slaved to (u_z, theta), and those
    are already continuous, so "phi_1 = phi_2" is implied. Zero gain means
    dropping the row discards nothing. Nonzero gain would mean it discards a
    real condition, and the reduction would be wrong for that material pair.

    Returns (gain, rank_without, rank_with, rank_with_constants).
    """
    B = full["B"]
    n_rows = B.rows - 1                      # drop the gauge-fixing row
    # columns of the physical modes only; the trailing ones are the constants
    n_modes = len(full["a1"]) + len(full["a2"])
    phi_idx = STATE_FIELDS.index("phi")

    def sub(rows, ncols):
        A = mp.matrix(len(rows), ncols)
        for r, i in enumerate(rows):
            for j in range(ncols):
                A[r, j] = B[i, j]
        return A

    mode_rows = [i for i in range(n_rows)]
    rows_wo = [i for i in mode_rows if i != phi_idx]
    r_wo, _ = rank_at(sub(rows_wo, n_modes))
    r_w, _ = rank_at(sub(mode_rows, n_modes))
    r_full, _ = rank_at(B)
    return r_w - r_wo, r_wo, r_w, r_full
