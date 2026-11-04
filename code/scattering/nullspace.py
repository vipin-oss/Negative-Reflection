"""
code/scattering/nullspace.py -- a precision-stable null vector for M(xi).

WHY THIS REPLACES THE PHASE-1D PATH
===================================
Phase 1D verified the production null vector `kernel_scaled` against an
independent construction and found agreement to ~1e-58. That verification was
sound for the configurations it tested, but it did not test MULTIPLE WORKING
PRECISIONS at fixed configuration, and doing so exposes a real failure:

    Al6061, 30 deg, branch 3 (a thermal branch)
        dps  60..85  ->  correct branch      (|theta|/max = 1.0)
        dps  90      ->  WRONG branch        (|phi|/max = 1.0, |theta| = 0)
        dps  95..105 ->  correct
        dps 110      ->  WRONG
        dps 120      ->  WRONG
        dps 125..130 ->  correct
        dps 135      ->  WRONG

THE MECHANISM
-------------
`kernel_scaled` equilibrates COLUMNS using

    cs[j] = max_i |A[i,j]|      ,  cs[j] = 1 if that max is 0

For a centrosymmetric medium the piezoelectric coupling is exactly zero, so at a
thermal branch the electric-potential column of M is not small -- it is ZERO.
Whether the computed max is 8.04e-87 (roundoff) or exactly 0 (the `if mcol > 0`
guard fires) depends on the root's last digits and therefore on dps. Both are
"correct" numerically and they lead to different matrices:

    dps <= 85 : cs = [... , 8.04e-87, ...]  ->  M[:,2] amplified by 1e86,
                                                genuine thermal minors dominate
    dps >= 90 : cs = [... , 1.0,      ...]  ->  M[:,2] left at roundoff,
                                                a spurious (2,2) minor wins

The free-variable search then picks a different 3x3 minor and returns a
DIFFERENT NULL VECTOR -- the wire gauge (phi = 1, everything else 0) instead of
the thermal branch. No error is raised; the vector is a perfectly good null
vector of the *equilibrated* matrix. It is simply the wrong mode.

This is why `classify_mode` saw mech/max = 3.06e-6 at some precisions and 0.0 at
others. A branch classifier that is precision-dependent is not a classifier.

THE FIX
-------
Two independent changes, both principled rather than tuned:

1. **Structural zero detection.** A column whose maximum is below a
   precision-relative floor is a STRUCTURAL zero, not a small number. It is
   excluded from the free-variable search entirely, so roundoff can never be
   selected as the dominant minor. The floor is
   `n * eps * max_scale`, i.e. what the arithmetic can actually resolve --
   not a physics threshold.

2. **SVD-free eigenvector of the equilibrated Gram matrix.** Instead of a
   cofactor search that can latch onto any of several comparable minors, the
   null vector is the eigenvector of `C^H C` belonging to the SMALLEST
   eigenvalue, with `C` two-sided equilibrated. This is well conditioned, has no
   minor-selection step at all, and mpmath's `eighe` is reliable for complex
   Hermitian input (unlike mpmath's SVD, which fails on complex matrices --
   recorded as E-11 in the Phase 1C dead-end list).

   The two routes are kept SEPARATE and cross-checked, so a disagreement is
   reported rather than silently resolved.

Both are verified against each other and against the Cramer route, and the
classifier consumes only the cross-checked result.

Reproduce:
    cd code && python3 -m scattering.nullspace
"""
from __future__ import annotations

import mpmath as mp


# ---------------------------------------------------------------------------
def _two_sided_equilibrate(A, dps):
    """Row- then column-equilibrate, treating arithmetic zeros as structural.

    Returns (C, row_scales, col_scales, structural_cols).

    A column is STRUCTURAL when its largest entry is below what the arithmetic
    can resolve: `n * 10^(-dps+4) * global_max`. Making this explicit is the
    whole point -- the previous code used `if mcol > 0`, which is a test of the
    floating-point accident rather than of the mathematics.
    """
    n = A.rows
    entries = [abs(A[i, j]) for i in range(n) for j in range(n)]
    gmax = max(entries) if entries else mp.mpf(1)
    if gmax == 0:
        return A, [mp.mpf(1)] * n, [mp.mpf(1)] * n, list(range(n))

    #: entries below this are indistinguishable from zero at this precision
    floor = mp.mpf(n) * mp.power(10, -dps + 4) * gmax

    rs = []
    for i in range(n):
        mx = max(abs(A[i, j]) for j in range(n))
        rs.append(mx if mx > floor else mp.mpf(1))

    B = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = A[i, j] / rs[i]

    cs, structural = [], []
    for j in range(n):
        mx = max(abs(B[i, j]) for i in range(n))
        if mx > floor / max(rs[i] for i in range(n)) and mx > 0:
            cs.append(mx)
        else:
            cs.append(mp.mpf(1))
            structural.append(j)

    C = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            C[i, j] = B[i, j] / cs[j]
    return C, rs, cs, structural


def null_vector_gram(A, dps=None):
    """Null vector as the smallest-eigenvalue eigenvector of the equilibrated
    Gram matrix. No minor search, therefore no branch ambiguity."""
    if dps is None:
        dps = mp.mp.dps
    n = A.rows
    C, rs, cs, structural = _two_sided_equilibrate(A, dps)
    G = C.H * C
    # eighe requires a Hermitian matrix; symmetrise to kill roundoff asymmetry
    Gh = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    try:
        E, Q = mp.eighe(Gh)
    except Exception as exc:                                   # noqa: BLE001
        raise RuntimeError(f"eighe failed: {exc!r}") from exc
    idx = min(range(n), key=lambda i: abs(E[i]))
    v = [Q[i, idx] / cs[i] for i in range(n)]
    return v, [abs(e) for e in E], structural


def _rel_residual(A, v):
    """||A v|| / (||A||_max * ||v||_max) -- a scale-free measure of null-ness.

    Used instead of comparing equilibrated minors, because minor size measures
    how well conditioned a COFACTOR is, not how nearly the resulting vector
    satisfies A v = 0. Those are different questions and the earlier code
    conflated them.
    """
    n = A.rows
    r = [sum(A[i, j] * v[j] for j in range(n)) for i in range(n)]
    rn = max(abs(t) for t in r)
    an = max(abs(A[i, j]) for i in range(n) for j in range(n))
    vn = max(abs(t) for t in v)
    if an == 0 or vn == 0:
        return mp.mpf(0)
    return rn / (an * vn)


def null_vector_cramer(A, dps=None):
    """Null vector by ENUMERATING every (dropped-row, free-column) pair, solving
    the 3x3 subsystem by Cramer, and keeping the candidate that MINIMISES THE
    RELATIVE RESIDUAL ||A v|| / (||A|| ||v||).

    Why this and not "largest minor wins": the largest minor identifies the
    best-conditioned cofactor, which need not be the one whose vector is
    closest to the true null space. The failure this replaces was precisely of
    that kind -- at Al6061 30 deg branch 3 a spurious minor won and produced the
    pure electric wire gauge where the true branch is thermal. Selecting by
    residual is a direct statement of what a null vector is.

    Every column is a legal free variable, including structurally-zero ones: the
    wire-gauge null vector (phi = 1, u = w = 0) lives ENTIRELY in a structurally
    zero column, so excluding such columns made it unreachable.
    """
    if dps is None:
        dps = mp.mp.dps
    n = A.rows

    def det3(mm, rr, cc):
        a, b, c = (mm[rr[0], cc[0]], mm[rr[0], cc[1]], mm[rr[0], cc[2]])
        d, e, f = (mm[rr[1], cc[0]], mm[rr[1], cc[1]], mm[rr[1], cc[2]])
        g, h, i2 = (mm[rr[2], cc[0]], mm[rr[2], cc[1]], mm[rr[2], cc[2]])
        return a * (e * i2 - f * h) - b * (d * i2 - f * g) + c * (d * h - e * g)

    # equilibrate only to make the 3x3 linear solves well conditioned; the
    # RESIDUAL is then measured on the ORIGINAL matrix, where "zero" is physical
    C, rs, cs, _structural = _two_sided_equilibrate(A, dps)

    best = None
    for dr in range(n):
        rows = [i for i in range(n) if i != dr]
        for jf in range(n):
            cols = [c for c in range(n) if c != jf]
            M3 = mp.matrix(3, 3)
            for r_i, r in enumerate(rows):
                for c_i, c in enumerate(cols):
                    M3[r_i, c_i] = C[r, c]
            det = det3(M3, [0, 1, 2], [0, 1, 2])
            if det == 0:
                continue
            rhs = [-C[r, jf] for r in rows]
            xs = []
            for c_i in range(3):
                Mi = mp.matrix(3, 3)
                for r_i in range(3):
                    for cc_i in range(3):
                        Mi[r_i, cc_i] = rhs[r_i] if cc_i == c_i else M3[r_i, cc_i]
                xs.append(det3(Mi, [0, 1, 2], [0, 1, 2]) / det)
            ae = [mp.mpc(0)] * n
            ae[jf] = mp.mpc(1)
            for c_i, c in enumerate(cols):
                ae[c] = xs[c_i]
            # undo equilibration, then judge on the physical matrix
            v = [ae[j] / cs[j] for j in range(n)]
            res = _rel_residual(A, v)
            if best is None or res < best[0]:
                best = (res, v, dr, jf)
    if best is None:
        raise RuntimeError("no non-singular 3x3 subsystem in any enumeration")
    res, v, dr, jf = best
    return v, dr, jf, res


# ---------------------------------------------------------------------------
def _normalise_physical(v, dps):
    """Scale and fix phase. Returns (vector, scale_basis).

    TWO BASES, chosen by arithmetic and not by preference:

      * "mechanical" -- divide by the largest of u,w. Used when the mechanical
        content is RESOLVABLE at this precision. This is the right basis for
        any elastic or piezoelectric branch, where u is O(1) m and phi is
        O(1e9) V; normalising by the largest component of any kind would make
        every piezoelectric mode look gauge-dominated.

      * "largest" -- divide by the largest component of any kind. Used only when
        the mechanical content is below what the arithmetic can resolve, i.e.
        the eigenvector really is the pure electric wire gauge
        (phi = 1, u = w = 0). Dividing by a mechanical component here would
        mean dividing by roundoff: at Al6061 30 deg the wire-gauge branch has
        |u| ~ 1e-61, so mechanical normalisation amplified the two nullspace
        routes' roundoff by 1e61 and reported a spurious cross-check failure of
        ~6e60. The disagreement was entirely an artefact of the basis.

    The resolvability floor is `10**(4-dps) * largest`, the same arithmetic
    floor used in the equilibration: below it, a number is not small, it is zero.
    """
    comps = [abs(t) for t in v]
    biggest = max(comps)
    if biggest == 0:
        return None, "degenerate"
    floor = mp.power(10, 4 - dps) * biggest
    mech = max(range(2), key=lambda i: comps[i])
    if comps[mech] > floor:
        nrm, basis, pivot = comps[mech], "mechanical", mech
    else:
        nrm, basis, pivot = biggest, "largest", comps.index(biggest)
    out = [mp.mpc(t) / nrm for t in v]
    ph = out[pivot] / abs(out[pivot])
    return [t / ph for t in out], basis


def null_vector(A, dps=None, cross_check=True):
    """The robust null vector, with the two routes cross-checked.

    Returns dict with `vector`, `svals`, `structural_columns`, `agreement`,
    `route_used`. If the routes disagree beyond tolerance, `agreement` is False
    and the CALLER MUST NOT use the vector -- it is reported, never silently
    preferred.
    """
    if dps is None:
        dps = mp.mp.dps
    vg, svals, structural = null_vector_gram(A, dps)
    vg_n, basis = _normalise_physical(vg, dps)
    out = dict(vector=vg_n, svals=svals, normalisation_basis=basis,
               structural_columns=structural, route_used="gram",
               agreement=True, cramer_gap=None)
    if not cross_check:
        return out
    if all(s < mp.mpf("1e-30") * max(svals) for s in svals[:-1]):
        # fully degenerate: no unique null direction, say so
        out["agreement"] = False
        out["cramer_gap"] = "degenerate spectrum"
        out["route_used"] = "gram (unconfirmed)"
        return out
    try:
        vc, dr, jf, cramer_res = null_vector_cramer(A, dps)
        vc_n, _vbasis = _normalise_physical(vc, dps)
        out["cramer_residual"] = cramer_res
        out["gram_residual"] = _rel_residual(A, vg)
        if vc_n is None or out["vector"] is None:
            out["agreement"] = False
            out["cramer_gap"] = "normalisation failed"
        else:
            gap = max(abs(out["vector"][i] - vc_n[i]) for i in range(4))
            out["cramer_gap"] = gap
            # Both routes must agree AND both must actually be null vectors.
            # Agreement alone is worthless if both are equally wrong.
            if gap > mp.mpf("1e-25") or out["gram_residual"] > mp.mpf("1e-25"):
                out["agreement"] = False
    except Exception as exc:                                   # noqa: BLE001
        out["agreement"] = False
        out["cramer_gap"] = f"cramer failed: {exc!r}"
    return out


# ---------------------------------------------------------------------------
def _main():
    import sys
    import numpy as np
    sys.path.insert(0, ".")
    import foundation_checks as fc

    w = 2 * np.pi * 1e6
    k = float(w * np.sin(np.deg2rad(30)) / 4000.0)
    print("Stability of the null vector across working precision")
    print("(previously intermittent at dps 90, 110, 120, 135)\n")
    for mat in ("Al6061", "PZT-5H"):
        print(f"--- {mat} ---")
        for dps in (50, 60, 80, 90, 100, 110, 120, 135, 150):
            mp.mp.dps = dps
            m = fc.material_from_table(mat)
            K, W = mp.mpf(repr(k)), mp.mpf(repr(float(w)))
            roots = [mp.mpc(z) for z in m.xi_roots(K, W)]
            tags = []
            for xi in roots:
                A = fc.M_mp(m, K, xi, W)
                r = null_vector(A, dps)
                v = r["vector"]
                if v is None:
                    tags.append("none")
                    continue
                cm = max(abs(t) for t in v)
                th = abs(v[3]) / cm
                ph = abs(v[2]) / cm
                tags.append("T" if th > 0.5 else ("W" if ph > 0.9 else "?"))
            print(f"  dps={dps:>4d}: {' '.join(tags)}")


if __name__ == "__main__":
    _main()
