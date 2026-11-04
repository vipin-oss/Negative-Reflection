"""
code/scattering/roots.py -- dispersion roots VERIFIED against the matrix.

WHY THIS EXISTS
===============
Task 2R was opened to investigate a specific claim: that `Medium6mm.xi_roots`,
which obtains the eight roots from the characteristic POLYNOMIAL via
`mp.polyroots`, returns inaccurate roots, and that the inaccuracy leaves the
characteristic matrix M non-singular so that every nullspace routine silently
returns the most-compressed direction as a spurious "null vector".

THE FIRST HYPOTHESIS WAS WRONG, AND THE WAY IT WAS WRONG MATTERS
----------------------------------------------------------------
The probe that produced the motivating numbers -- xi_roots off by 4.93e-14,
det M(xi_roots) = 0.009219 - 279.759j, M^H M eigenvalues [0.42, 1.0, 1.0, 2.58]
-- built its "exact" reference as

    cT = sqrt(m.C44 / m.rho)

with m.C44 and m.rho both FLOAT64. That ratio is itself only accurate to 1.2e-13
as a distance in wavenumber. Measured at the working precision instead:

    dps   min |xi_roots - w/cT|
     30                2.0e-28
     40                0.0
     50                2.7e-48
     60                0.0

`xi_roots` is EXACT to working precision. There is no root error. The
"4.93e-14" was the error in the reference, not in the root, and the
"det M = 0.0092 - 279.759j" was that same reference error amplified by C44.

Consequence: a Newton polish against det M, built to repair a fault that was
never there, is not neutral -- it is destructive. From the exact root it moves
xi by +2.0e-17 imaginary and takes

    M[0,0]:  1.65e-24  ->  -2.13e-3

turning a structurally-zero row into a loud one and destroying the nullspace it
was meant to protect. That polish has been removed. Do not reintroduce it
without re-reading this section.

THE ACTUAL DEFECTS (both real, both measured)
============================================
1. THE METRIC. `nullspace._rel_residual` divides by the GLOBAL matrix max.
   M's entries span 1e-5 to 1e19 because rows are equations in different units.
   A candidate that violates the electric-potential row by 3.6e-5 therefore
   scores 3.6e-5/1.99e19 = 1.8e-24 and passes. That is how the pure gauge
   vector (phi = 1, u = 0) came back as "the mode" on a branch whose physical
   content is (u_x = 1, u_z = 0).

2. THE EQUILIBRATION. `_two_sided_equilibrate` rescales rows, which lifts a
   roundoff-level row to unit size and lets it dominate the selection. And an
   eigenvector search on A^H A, however scaled, picks the direction favoured by
   the numerically loudest row.

WHAT THIS MODULE DOES INSTEAD
=============================
* A null vector is found by ENUMERATION over (dropped row, free column) and
  judged by a ROW-RELATIVE residual on the ORIGINAL matrix: every row is
  normalised by its own largest entry, and rows whose largest entry is at the
  arithmetic floor are skipped as structural zeros. That last clause is load
  bearing -- without it, the shear root's row 0 (value 1.65e-24, itself the
  roundoff of a 1.07e17 subtraction) is scored as 1.0 and REJECTS the correct
  physical answer.

* Nothing is returned unverified. `null_vector_verified` reports `ok=False`
  rather than a compressed direction, and `verified_roots` marks a root usable
  only if a null vector is actually demonstrated.

* The floor that separates "structurally zero" from "small but real" is tied to
  the WORKING PRECISION, not to a physics scale. Measured: the shear branch's
  row-0 entry equals the roundoff floor and tracks 10^-dps (2.8e-14 at dps 30,
  1.7e-24 at 40, 1.1e-44 at 60).

Reproduce:
    cd code && python3 -m scattering.roots
"""
from __future__ import annotations

import mpmath as mp

#: A null residual below this is accepted. Measured floor across every branch
#: of both media at 40 dps is 1.8e-24; this sits four orders above the floor and
#: sixteen orders below any physically meaningful contamination. It answers
#: "is this vector in the nullspace", which is a question about arithmetic.
NULL_TOL = mp.mpf("1e-20")

#: Digits of headroom above the working-precision floor when deciding that a
#: row is structurally zero. Five digits is generous: the measured roundoff
#: tracks 10^-dps to within one or two digits.
FLOOR_DIGITS = 5


# ---------------------------------------------------------------------------
# the metric
# ---------------------------------------------------------------------------
def row_relative_residual(A, v):
    """max_i |(A v)_i| / (max_j |A_ij| * ||v||_inf), skipping structural zeros.

    THE METRIC THAT DECIDES WHETHER A VECTOR IS A NULL VECTOR.

    Dividing by the global matrix max -- the natural first choice, and what
    `nullspace._rel_residual` does -- is wrong for M. Its rows are equations in
    different units and span 1e-5 to 1e19, so a violation of the quietest row
    is invisible against the loudest. Normalising each ROW by itself asks the
    only question that matters: can this candidate satisfy EVERY equation?

    Rows whose largest entry is at the working-precision floor are skipped.
    They read 0 = 0 for every vector and cannot reject anything. This clause is
    not a convenience: at k = 0 on the shear branch, row 0 evaluates to
    1.65e-24 against a row max of 1.65e-24 -- a ratio of 1.0 -- and treating it
    as a real equation rejects the physically correct shear vector.
    """
    n, m = A.rows, A.cols
    vn = max(abs(t) for t in v)
    if vn == 0:
        return mp.inf
    entries = [abs(A[i, j]) for i in range(n) for j in range(m)]
    gmax = max(entries) if entries else mp.mpf(0)
    if gmax == 0:
        return mp.mpf(0)
    floor = mp.power(10, -(mp.mp.dps - FLOOR_DIGITS)) * gmax
    worst = mp.mpf(0)
    for i in range(n):
        rowmax = max(abs(A[i, j]) for j in range(m))
        if rowmax <= floor:
            continue
        val = sum(A[i, j] * v[j] for j in range(m))
        worst = max(worst, abs(val) / (rowmax * vn))
    return worst


def structural_zero_rows(A):
    """Indices of rows that are zero at the working precision."""
    n, m = A.rows, A.cols
    entries = [abs(A[i, j]) for i in range(n) for j in range(m)]
    gmax = max(entries) if entries else mp.mpf(0)
    if gmax == 0:
        return list(range(n))
    floor = mp.power(10, -(mp.mp.dps - FLOOR_DIGITS)) * gmax
    return [i for i in range(n)
            if max(abs(A[i, j]) for j in range(m)) <= floor]


# ---------------------------------------------------------------------------
# rank
# ---------------------------------------------------------------------------
def matrix_rank_nullity(A, tol=None):
    """Rank and nullity of a COMPLEX matrix, reported with their tolerance.

    Two things make this harder than it looks, and both were got wrong first.

    NOT on the raw matrix. M mixes rows whose natural sizes differ by 30+
    orders of magnitude -- a displacement equation is O(1e19), the electric
    potential equation of a centrosymmetric medium is literally 0, the thermal
    row is O(1e31). A relative singular-value test on that says nothing; every
    entry of a perfectly good matrix is tiny next to the largest one. The first
    version of this function reported nullity 4 for all eight branches at 30
    degrees for exactly that reason.

    NOT by mp.svd or mp.cond. Both fail or return +inf on complex input
    (recorded dead ends). mp.eighe on the Hermitian A^H A is used instead, and
    A^H A is explicitly symmetrised because mpmath's is not guaranteed Hermitian
    to the last digit.

    Rows are equilibrated to unit max before the singular values are taken.
    That preserves singularity exactly -- the scaling is by invertible diagonal
    matrices -- and it is what makes a relative tolerance meaningful.
    """
    n, m = A.rows, A.cols
    if tol is None:
        tol = mp.power(10, -(mp.mp.dps - 10))
    B = _row_equilibrated(A)
    G = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            G[i, j] = sum(mp.conj(B[r, i]) * B[r, j] for r in range(n))
    Gh = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    try:
        ev, _Q = mp.eighe(Gh)
    except Exception:                                          # noqa: BLE001
        return -1, -1, [], tol
    sv = sorted((mp.sqrt(abs(e)) for e in ev), reverse=True)
    smax = sv[0] if sv else mp.mpf(0)
    if smax == 0:
        return 0, m, sv, tol
    rank = min(sum(1 for s in sv if s > tol * smax), min(n, m))
    return rank, m - rank, sv, tol


def _row_equilibrated(A):
    """A with every row divided by its largest entry.

    Singularity is preserved exactly: the scaling is by an invertible diagonal
    matrix. Rows that are identically zero are left alone, which is correct --
    they are genuinely structural and register as exact null directions.
    """
    n, m = A.rows, A.cols
    B = mp.matrix(n, m)
    for i in range(n):
        mx = max(abs(A[i, j]) for j in range(m))
        d = mx if mx != 0 else mp.mpf(1)
        for j in range(m):
            B[i, j] = A[i, j] / d
    return B


# ---------------------------------------------------------------------------
# the null vector
# ---------------------------------------------------------------------------
def null_vector_verified(A, dps=None):
    """The null vector of A, or `ok=False` -- never a compressed direction.

    FOUND BY ENUMERATION, NOT BY AN EIGENVECTOR. Every eigenvector route failed
    on this matrix for a structural reason: whichever way A is scaled, the
    smallest singular direction of A^H A is dominated by the numerically
    loudest row, and at k = 0 on the shear branch that direction is the pure
    gauge (phi = 1), not the physical shear state (u_x = 1).

    So the candidates are enumerated. For a rank-(n-1) matrix, dropping row i
    and solving for the remaining entries with column j pinned to 1 gives a
    candidate for every (i, j) pair -- sixteen for a 4x4. Each is scored by
    `row_relative_residual` on the ORIGINAL matrix, and the best is returned
    WITH its score. A wrong answer is then visible as a poor score rather than
    silently accepted.

    Returns a dict:
        vector, ok, residual, rank, nullity, singular_values,
        structural_zero_rows, normalisation_basis, note
    """
    if dps is None:
        dps = mp.mp.dps
    n, m = A.rows, A.cols
    rank, nullity, sv, tol = matrix_rank_nullity(A)

    best = None
    for drop in range(n):
        rows = [i for i in range(n) if i != drop]
        for free in range(m):
            cols = [c for c in range(m) if c != free]
            if len(cols) != n - 1 or len(rows) != n - 1:
                continue
            sub = mp.matrix(n - 1, n - 1)
            for r_i, r in enumerate(rows):
                for c_i, c in enumerate(cols):
                    sub[r_i, c_i] = A[r, c]
            try:
                xs = mp.lu_solve(sub, mp.matrix([-A[r, free] for r in rows]))
            except Exception:                                  # noqa: BLE001
                continue
            cand = [mp.mpc(0)] * m
            cand[free] = mp.mpc(1)
            for c_i, c in enumerate(cols):
                cand[c] = mp.mpc(xs[c_i])
            sc = row_relative_residual(A, cand)
            if best is None or sc < best[0]:
                best = (sc, cand)

    if best is None:
        return {"vector": None, "ok": False, "residual": None,
                "rank": rank, "nullity": nullity, "singular_values": sv,
                "structural_zero_rows": structural_zero_rows(A),
                "normalisation_basis": None,
                "note": "no candidate subsystem was solvable"}

    residual, v = best
    v = _normalise(v)
    return {"vector": v, "ok": residual <= NULL_TOL, "residual": residual,
            "rank": rank, "nullity": nullity, "singular_values": sv,
            "structural_zero_rows": structural_zero_rows(A),
            "normalisation_basis": "mechanical",
            "note": ("verified null vector" if residual <= NULL_TOL else
                     "RESIDUAL ABOVE FLOOR: this is not a null vector")}


def _normalise(v):
    """Scale so the largest MECHANICAL component is 1.

    Displacement components (indices 0, 1) are preferred because the physical
    state has natural units there; psi and theta are the fallbacks. An arbitrary
    component would make the reported amplitudes unit-dependent.
    """
    for i in (0, 1):
        if i < len(v) and abs(v[i]) > 0:
            d = v[i]
            return [t / d for t in v]
    for i in range(len(v)):
        if abs(v[i]) > 0:
            d = v[i]
            return [t / d for t in v]
    return v


# ---------------------------------------------------------------------------
# roots
# ---------------------------------------------------------------------------
def verified_roots(material, k, w, dps=None):
    """All eight roots of the dispersion relation, each verified against M.

    NO POLISHING. See the module docstring: the root error that motivated a
    polish was an artifact of a float64 reference, and polishing against det M
    moves an exact root away from exactness. `xi_roots` is used as given and the
    result is VERIFIED, which is the part that was actually missing before.

    A root whose matrix admits no demonstrated null vector is returned with
    `ok=False` and must not be used to build a mode.
    """
    if dps is None:
        dps = mp.mp.dps
    out = []
    for br, xi in enumerate(material.xi_roots(k, w)):
        A = _matrix_at(material, k, xi, w)
        nv = null_vector_verified(A, dps=dps)
        nv = dict(nv)
        nv["branch"] = br
        nv["xi"] = xi
        out.append(nv)
    return out


def _matrix_at(material, k, xi, w):
    """The 4x4 characteristic matrix. Imported lazily to avoid a cycle."""
    import foundation_checks as fc
    return fc.M_mp(material, k, xi, w)


# ---------------------------------------------------------------------------
# self-check
# ---------------------------------------------------------------------------
def _main():
    import foundation_checks as fc
    mp.mp.dps = 40
    w = mp.mpf(2) * mp.pi * mp.mpf(10) ** 6

    print("Dispersion roots verified against the MATRIX (no polishing)")
    print(f"dps = {mp.mp.dps}, w = 2*pi*1e6\n")

    for name in ("Al6061", "PZT-5H"):
        mat = fc.material_from_table(name)
        for theta in (0, 30):
            k = _k_from_angle(mat, w, theta)
            roots = verified_roots(mat, k, w)
            print(f"--- {name}, theta = {theta} deg ---")
            print("   br                    xi            ok  resid  "
                  "rank null   ux    uz  phi  th")
            for r in roots:
                v = r["vector"] or [mp.mpc(0)] * 4
                print(f"   {r['branch']:2d}  {mp.nstr(r['xi'], 18):>22s}  "
                      f"{'Y' if r['ok'] else 'n':>3s}  "
                      f"{mp.nstr(r['residual'], 3):>7s}  "
                      f"{r['rank']:4d} {r['nullity']:4d}  "
                      + "  ".join(f"{mp.nstr(abs(t), 3):>4s}" for t in v))
            print()


def _k_from_angle(material, w, theta_deg):
    """Horizontal wavenumber for a plane wave at `theta_deg` from the normal."""
    cL = mp.sqrt(mp.mpf(repr(material.C11)) / mp.mpf(repr(material.rho)))
    return w * mp.sin(mp.radians(theta_deg)) / cL


if __name__ == "__main__":
    _main()
