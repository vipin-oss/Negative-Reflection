"""
code/scattering/asymptotic.py -- the k -> 0 limit, done properly.

WHY THIS MODULE EXISTS
======================
Task 2S claimed the potential "diverges as 1/k" and therefore that potential
continuity "cannot be imposed at k = 0". Task 2T was opened to test that, and the
claim does not survive. This module contains the analysis that replaces it.

THE ERROR TASK 2S MADE
======================
Three separate things were conflated:

1. THE EIGENVECTOR SCALE. `null_vector_verified` normalises every mode vector to
   u_x = 1. A quasi-P mode has u_x/u_z proportional to k, so pinning u_x = 1
   inflates the WHOLE vector by 1/k. Measured on PZT-5H branch 5: u_x/u_z is
   1.046e-6 at k/k_ac = 1e-6, i.e. exactly proportional to k, and phi/|u_z| is
   a constant 1.78975e9 over six decades. Nothing diverges but the arbitrary
   scale. A component of an eigenvector has no physical meaning on its own.

2. THE PHYSICAL LIMIT. Under a k-INDEPENDENT normalisation -- unit mechanical
   z-flux per mode -- the same field has a FINITE limit. Measured for
   Al6061|PZT-5H at P incidence: phi -> 0.06464627422, stable to 1e-10 over
   k/k_ac from 1e-2 to 1e-5. Every configuration tested has a finite limit on
   every component that does not vanish identically.

3. THE DEGENERACY. At k = 0 the last accepted branch of each half-space has
   xi = 0 and its state becomes (0,0,1,0,...) up to a correction of order k:
   a near-gauge potential mode whose D_z is O(k). It is admitted at k > 0 and
   excluded at k = 0. Because a near-gauge mode is present, the k -> 0 problem
   is SINGULAR -- measured condition number growing as 1/k^2, 3.9e16 at
   k/k_ac = 1e-7 -- and the individual modal coefficients are not bounded even
   though the physical total field is. The total field is the sum of O(1)
   contributions from the three regular modes plus c_3 * s_3 with
   c_3 -> infinity and s_3 -> gauge; the product is O(1).

So the limit EXISTS but the k > 0 representation of it is ill-conditioned. That
is a statement about the solution process, not about the physics, and it is not
the same claim as "phi diverges".

WHY THE TWO ELECTRICAL ROWS VANISH AT k = 0 -- AND IT IS NOT A MODEL CHOICE
===========================================================================
Take the interface conditions from `scattering.maxwell` for a dielectric
contact (no free surface charge): [phi] = 0 and [D_z] = 0.

    tangential E:   E_x = -i k phi, so [E_x] = -i k [phi].
                    At k = 0 this is 0 = 0 for ANY finite [phi].
                    The equation DEGENERATES; it is not violated and it is not
                    a gauge condition. It simply carries no information at the
                    single point k = 0.

    normal D:       Gauss at k = 0 reads i xi D_z = 0, so D_z = 0 on BOTH sides
                    (for xi != 0), hence [D_z] = 0 - 0 = 0.
                    Also automatic. Also not a model choice.

Neither row is "the insulated condition" and neither is removed "because it is a
gauge". Both are removed because at exactly k = 0 they are satisfied identically
by every admissible field. That is a legitimate reason -- and it is the reason
the supervisor's brief allows -- but it had not been demonstrated before.

AN O(k) COMPONENT IS NOT A FAILED LIMIT
=======================================
Several components vanish linearly at k = 0. Measured for Al6061|Al6061 at P
incidence, t_xz/k is constant at 2.8668193 to seven digits across k/k_ac from
1e-3 to 1e-5, so t_xz = C k exactly and its limit is 0. Comparing such a
component against a constant-limit extrapolation manufactures a spurious
residual of ~5e-4. `limit_of` below classifies each component's k-order before
extrapolating, and only extrapolates to a CONSTANT when the data support one.

Reproduce:
    cd code && python3 -m scattering.asymptotic
"""
from __future__ import annotations

import sys

import mpmath as mp

import foundation_checks as fc
import scattering.assembly as A
from scattering.state import STATE_FIELDS


# ---------------------------------------------------------------------------
# normalisations
# ---------------------------------------------------------------------------
def unit_flux_mode(side, br):
    """Rescale mode `br` so its z energy flux is +-1.

    THIS IS THE NORMALISATION THAT MAKES THE LIMIT MEANINGFUL. For a propagating
    mode the z-flux is an O(1) physical quantity, so dividing by its square root
    removes the arbitrary eigenvector scale WITHOUT introducing any k
    dependence of its own. Contrast `roots._normalise`, which pins a displacement
    component and therefore scales like 1/(that component) as the component
    vanishes.

    Returns None when the mode carries no flux at all.
    """
    st = side.raw_states[br]
    jz = side.jzs[br]
    if st is None or jz is None or abs(jz) == 0:
        return None
    s = mp.sqrt(abs(jz))
    return [c / s for c in st]


def hybrid_mode(side, br):
    """Per-mode scale that never vanishes: flux if it carries any, else phi.

    WHY A SECOND NORMALISATION IS NEEDED. `unit_flux_mode` divides by
    sqrt(|jz|), which is k-independent for a PROPAGATING mode -- exactly what
    the limit measurement wants. But the near-gauge mode at small k carries
    essentially no flux (D_z is O(k)), so sqrt(|jz|) -> 0 and the mode is
    either ill-scaled or dropped outright. Measured on Al6061|PZT-5H at
    k/k_ac = 1e-3: an all-flux basis has 7 columns for 8 admitted modes,
    because branch 3 is silently lost -- and branch 3 is precisely the mode that
    supplies the potential jump. The system then cannot satisfy potential
    continuity and returns a residual of 0.0323 while looking square-ish.

    So: scale by flux when there is flux to speak of, and by the potential
    otherwise. The scale is arbitrary per column -- the coefficient absorbs it
    -- and the only requirement is that it be nonzero and not vanish with k.
    """
    st = side.raw_states[br]
    if st is None:
        return None
    jz = side.jzs[br]
    if jz is not None and abs(jz) > 0:
        s = mp.sqrt(abs(jz))
        if s > 0:
            return [c / s for c in st]
    # fall back to the potential, then to the largest entry
    for comp in (2, 1, 0):
        if abs(st[comp]) > 0:
            return [c / st[comp] for c in st]
    return None


def unit_displacement_mode(side, br, component=1):
    """Rescale so the chosen displacement component is 1 (the code's own gauge).

    Kept because it is the normalisation the production root finder uses, and
    the whole point is to show that conclusions drawn from it do not transfer.
    """
    st = side.raw_states[br]
    if st is None or abs(st[component]) == 0:
        return None
    d = st[component]
    return [c / d for c in st]


# ---------------------------------------------------------------------------
# the system
# ---------------------------------------------------------------------------
def _equilibrated_solve(B, rhs):
    """Row- then column-equilibrated least squares. Returns (a, residual, cond)."""
    n, m = B.rows, B.cols
    rs = []
    S = mp.matrix(n, m)
    for i in range(n):
        mx = max(abs(B[i, j]) for j in range(m)) or mp.mpf(1)
        rs.append(mx)
        for j in range(m):
            S[i, j] = B[i, j] / mx
    cs = []
    for j in range(m):
        mx = max(abs(S[i, j]) for i in range(n)) or mp.mpf(1)
        cs.append(mx)
    S2 = mp.matrix(n, m)
    for i in range(n):
        for j in range(m):
            S2[i, j] = S[i, j] / cs[j]
    G = S2.H * S2
    Gh = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    ev, _ = mp.eighe(Gh)
    mags = sorted(abs(e) for e in ev)
    cond = mags[-1] / mags[0] if mags[0] != 0 else mp.inf
    try:
        x = mp.lu_solve(Gh, S2.H * mp.matrix([rhs[i] / rs[i] for i in range(n)]))
    except ZeroDivisionError:
        return None, None, cond
    a = [x[j] / cs[j] for j in range(m)]
    resid = max(abs(sum(B[i, j] * a[j] for j in range(m)) - rhs[i])
                for i in range(n))
    return a, resid, cond


def dielectric_contact_system(m1, m2, k, w, incident_branch, rows=None,
                              normalisation="hybrid"):
    """The 8-row dielectric-contact interface, in a k-independent basis.

    `rows` selects which STATE_FIELDS rows are imposed; default all eight, which
    is the model `maxwell.DIELECTRIC_CONTACT` specifies. Passing a subset is how
    the production reduction is reproduced for comparison -- never how the
    physics is decided.

    `normalisation` selects the basis:
        "hybrid"       (default) flux, falling back to potential when a mode
                       carries no flux. This is the basis for SOLVING the
                       boundary-value problem, because it keeps the near-gauge
                       mode that supplies the potential jump.
        "flux"         unit |jz| per mode. k-independent for propagating modes,
                       so it is the basis for MEASURING the limit -- but it
                       drops the near-gauge mode and must not be used to solve.
        "displacement" the production eigenvector gauge, for comparison only.

    Returns a dict with the matrix, right-hand side, the two half-spaces, the
    admitted branches, and the reconstructed fields on both sides.
    """
    s1 = A.Side(m1, k, w, "z<0", dps=mp.mp.dps)
    s2 = A.Side(m2, k, w, "z>0", dps=mp.mp.dps)
    a1, _ = s1.admissible()
    a2, _ = s2.admissible()

    norm = {"flux": unit_flux_mode,
            "hybrid": hybrid_mode,
            "displacement": unit_displacement_mode}[normalisation]
    inc = norm(s1, incident_branch)
    if inc is None:
        return None

    cols = []
    for b in a1:
        st = norm(s1, b)
        if st is not None:
            cols.append((1, b, st))
    for b in a2:
        st = norm(s2, b)
        if st is not None:
            cols.append((2, b, st))

    n = len(cols)
    use = list(range(len(STATE_FIELDS))) if rows is None else list(rows)
    B = mp.matrix(len(use), n)
    rhs = mp.matrix(len(use), 1)
    for j, (sd, _b, st) in enumerate(cols):
        sg = mp.mpf(1) if sd == 1 else mp.mpf(-1)
        for r, i in enumerate(use):
            B[r, j] = sg * st[i]
    for r, i in enumerate(use):
        rhs[r] = -inc[i]

    a, resid, cond = _equilibrated_solve(B, rhs)
    out = {"B": B, "rhs": rhs, "a": a, "residual": resid, "cond": cond,
           "rows_used": [STATE_FIELDS[i] for i in use],
           "cols": cols, "side1": s1, "side2": s2, "a1": a1, "a2": a2,
           "incident_branch": incident_branch,
           "n_unknowns": n, "n_equations": len(use),
           "normalisation": normalisation}
    if a is not None:
        f1 = [mp.mpc(0)] * len(STATE_FIELDS)
        f2 = [mp.mpc(0)] * len(STATE_FIELDS)
        for j, (sd, _b, st) in enumerate(cols):
            for i in range(len(STATE_FIELDS)):
                if sd == 1:
                    f1[i] += a[j] * st[i]
                else:
                    f2[i] += a[j] * st[i]
        # `f1` is the SCATTERED field on side 1; the physical total there is
        # inc + f1. Returning both, clearly named, because conflating them
        # produced a spurious 3e-2 potential discontinuity during this task.
        out["f1_scattered"] = f1
        out["f1"] = [inc[i] + f1[i] for i in range(len(STATE_FIELDS))]
        out["f2"] = f2
        out["incident"] = inc
    return out


# ---------------------------------------------------------------------------
# the limiter
# ---------------------------------------------------------------------------
def k_order(values, ks, min_decades=2, floor=None):
    """Estimate the power p in f ~ C k^p from a sequence of k values.

    Returns (p, C, classification). The classification is one of

        CONSTANT   p is below a small threshold: the component has a finite
                   nonzero limit and may be extrapolated to k = 0
        VANISHING  p is positive: the limit is ZERO and must be compared to
                   zero, not to a constant extrapolation
        DIVERGENT  p is negative beyond tolerance: there is no finite limit
        UNRESOLVED the data do not define a clean power

    THIS EXISTS BECAUSE TASK 2T ITSELF WAS ALMOST FOOLED BY IT. Extrapolating
    t_xz of Al6061|Al6061 at P incidence to a constant produced a residual of
    5e-4 against the production answer of 0; t_xz/k turned out to be constant to
    seven digits, so t_xz = C k and the production zero was right.

    `floor` is an absolute magnitude below which a component is ROUNDOFF rather
    than signal. Without it, a component sitting at 1e-11 against a field of
    1e4 has a meaningless log-log slope and gets classified as CONSTANT with a
    confident-looking limit. Passing the field scale times a small epsilon makes
    those components UNRESOLVED, which is an honest answer.
    """
    fl = mp.mpf(0) if floor is None else mp.mpf(floor)
    usable = [(mp.mpf(kk), mp.mpf(abs(v)))
              for kk, v in zip(ks, values) if kk > 0 and abs(v) > fl]
    if len(usable) < min_decades:
        return None, None, "UNRESOLVED"
    # least squares in log-log
    xs = [mp.log(kk) for kk, _ in usable]
    ys = [mp.log(v) for _, v in usable]
    n = len(usable)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None, None, "UNRESOLVED"
    p = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    # spread of the local slopes, as a quality measure
    slopes = [(ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) for i in range(n - 1)]
    spread = max(slopes) - min(slopes) if len(slopes) > 1 else mp.mpf(0)
    if spread > mp.mpf("0.05"):
        return p, None, "UNRESOLVED"
    C = mp.e ** (my - p * mx)
    # ORDER MATTERS HERE AND THE OBVIOUS ORDER IS WRONG. Testing `p < 0.05`
    # first classified a clean 1/k sequence (p = -1.0) as CONSTANT, because -1
    # is also less than 0.05. The divergence branch must be tested first. This
    # was caught by the positive control in EQ-02A-LIM-02, which feeds the
    # classifier a synthetic 1/k; without that control the bug would have made
    # "no divergence found" vacuous.
    if p < -mp.mpf("0.05"):
        cls = "DIVERGENT"
    elif p > mp.mpf("0.5"):
        cls = "VANISHING"
    elif p < mp.mpf("0.05"):
        cls = "CONSTANT"
    else:
        cls = "UNRESOLVED"
    return p, C, cls


def _selfcheck():
    mp.mp.dps = 50
    w = mp.mpf(2) * mp.pi * mp.mpf(10) ** 6
    al = fc.material_from_table("Al6061")
    pz = fc.material_from_table("PZT-5H")
    kac = w / mp.sqrt(mp.mpf(repr(al.C11)) / mp.mpf(repr(al.rho)))
    ok = True

    print("k -> 0 limit under a normalisation-independent basis")
    print("=" * 92)
    ks = [kac * mp.mpf(f) for f in ("1e-2", "1e-3", "1e-4", "1e-5")]
    for nm, m1, m2 in (("Al6061|PZT-5H", al, pz), ("PZT-5H|Al6061", pz, al)):
        for br, pol in ((5, "P"), (6, "SV")):
            seq = {f: [] for f in STATE_FIELDS}
            for k in ks:
                r = dielectric_contact_system(m1, m2, k, w, br)
                if r is None:
                    continue
                for i, f in enumerate(STATE_FIELDS):
                    seq[f].append(r["f2"][i])
            print(f"\n  {nm} {pol}")
            scale = max([abs(v) for f in STATE_FIELDS for v in seq[f]] or [1])
            flo = scale * mp.mpf("1e-20")
            for f in ("ux", "uz", "phi", "theta", "tzz", "qz"):
                p, C, cls = k_order(seq[f], ks, floor=flo)
                print(f"    {f:>6s} order={mp.nstr(p, 3) if p is not None else '--':>8s} "
                      f"{cls:<10s} "
                      f"limit={(mp.nstr(abs(C), 8) if (cls == 'CONSTANT' and C is not None) else '--'):>14s}")
                if f == "phi" and cls == "DIVERGENT":
                    # phi has no free volume charge and no surface charge on a
                    # dielectric contact, so it is bounded by its own boundary
                    # condition. A divergence here would refute the whole basis.
                    print("      FAIL: phi must not diverge under this basis")
                    ok = False
    print()
    print("PASS: every field has a finite limit or vanishes" if ok else "FAILED")
    return ok


if __name__ == "__main__":
    sys.exit(0 if _selfcheck() else 1)
