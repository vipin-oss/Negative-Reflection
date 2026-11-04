"""
Independent isotropic elastic P-SV reference -- Phase 2A task 2, prerequisite 1.

WHY THIS FILE EXISTS, AND WHY IT IS NOT code/elastic_psv.py
===========================================================
`code/elastic_psv.py` was the Phase 1B reference. Task 2 audited it before
building anything on top, and it is WRONG IN ITS TRANSMITTED AMPLITUDES:

    welded Al/steel, normal incidence
      closed form  t = 2 Z1/(Z1+Z2) = +0.537019003885
      elastic_psv  tP              = -0.537019003885     <-- sign flipped

At oblique incidence the same flip appears in BOTH transmitted amplitudes at
every angle tested, while the reflected amplitudes are correct. Its energy
tests could not see this because R and T depend on |amplitude|^2 and therefore
on the sign only through a square. The Phase 1B validation was entirely
energy-based, so a sign error in half of the unknowns was invisible to it.

This is exactly the failure the task brief warns about -- "compare complex
amplitudes, not only total energy" -- and it is the reason the coupled solver
must be validated against amplitudes, not just against closure.

So the reference here is built from closed-form theory and verified against
closed-form expressions before being used to judge anything else:

  * normal incidence      Rayleigh:  r = (Z2-Z1)/(Z2+Z1), t = 2Z1/(Z2+Z1)
  * oblique incidence     Snell + the 4x4 welded system, assembled here
  * both are checked against the exact Fresnel/energy identities

CONVENTIONS (stated so they can be checked against any textbook)
===============================================================
    time dependence   e^{-i w t}
    spatial           exp(i(k x + xi z))
    displacement      u = (u_x, u_z), compressional potential along (k, xi)
    tractions         t_xz = sigma_xz, t_zz = sigma_zz
    cycle-averaged    <P_z> = -(1/2) Re[ sigma_xz v_x^* + sigma_zz v_z^* ]
    velocities        v = -i w u   (from the e^{-i w t} convention)

`<P_z>` takes a MINUS sign because with e^{-i w t} the physical field is
Re[A e^{i(...)}] and the time average of a product of two such fields picks up
Re[A B^*]/2 while the flux vector is -sigma . v.

Reproduce:
    cd code && python3 -m scattering.elastic_reference
"""
from __future__ import annotations

import mpmath as mp

#: below this the transmitted amplitude's imaginary part is treated as exactly
#: zero rather than as numerical noise. Not a physics threshold: it only decides
#: whether a printed value is expressed in real or complex form.
IMAG_SNAP = mp.mpf("1e-40")


class ISO:
    """Isotropic elastic solid, in mpmath throughout.

    Parameters are given as C11, C44 and rho. For an isotropic solid
    C33 = C11, C13 = C12 = C11 - 2*C44.
    """

    def __init__(self, name, rho, C11, C44):
        self.name = name
        self.rho = mp.mpf(repr(float(rho)))
        self.C11 = mp.mpf(repr(float(C11)))
        self.C44 = mp.mpf(repr(float(C44)))
        self.C12 = self.C11 - 2 * self.C44
        self.cL = mp.sqrt(self.C11 / self.rho)
        self.cT = mp.sqrt(self.C44 / self.rho)
        #: acoustic impedances for the closed-form normal-incidence limits
        self.ZL = self.rho * self.cL
        self.ZT = self.rho * self.cT

    def xi_P(self, k, w):
        """Axial wavenumber of the P branch, BEFORE branch selection.

        sqrt of a possibly negative real number returns the principal complex
        root, whose imaginary part is >= 0 by construction. That is the DECAYING
        (+z) branch; the -z branch is its negative. Selection between them is
        the caller's job and is done on flux, never on the sign of a square root.
        """
        return mp.sqrt((w / self.cL) ** 2 - k ** 2 + mp.mpc(0))

    def xi_S(self, k, w):
        return mp.sqrt((w / self.cT) ** 2 - k ** 2 + mp.mpc(0))

    def __repr__(self):
        return (f"ISO({self.name}, cL={mp.nstr(self.cL, 6)}, "
                f"cT={mp.nstr(self.cT, 6)}, ZL={mp.nstr(self.ZL, 6)})")


def mode_vector(k, xi, kind):
    """Unit-modulus displacement polarisation.

    P: along the wavevector (k, xi).   S: perpendicular, (xi, -k).

    Normalising to unit modulus is legitimate HERE and only here: the two
    components have the same dimension (metres), so a unit-modulus vector is a
    pure scale choice. In the coupled problem u, phi and theta have different
    units and no such normalisation exists.
    """
    if kind == "P":
        v = (mp.mpc(k), mp.mpc(xi))
    elif kind == "S":
        v = (mp.mpc(xi), mp.mpc(-k))
    else:
        raise ValueError(kind)
    n = mp.sqrt(abs(v[0]) ** 2 + abs(v[1]) ** 2)
    if n == 0:
        return None
    return (v[0] / n, v[1] / n)


def traction(m, k, xi, u):
    """(t_xz, t_zz) from the isotropic Hooke law, hand-expanded.

    sigma_zz = lambda e_xx + (lambda + 2 mu) e_zz
    sigma_xz = 2 mu e_xz
    with lambda = C12, mu = C44 and plane strain e_yy = 0.
    """
    I = mp.mpc(0, 1)
    lam, mu = m.C12, m.C44
    ux, uz = u
    exx = I * k * ux
    ezz = I * xi * uz
    exz = I * (xi * ux + k * uz) / 2
    tzz = lam * exx + (lam + 2 * mu) * ezz
    txz = 2 * mu * exz
    return txz, tzz


def flux_z(m, k, xi, u, w):
    """Cycle-averaged z power flux per unit displacement amplitude squared."""
    txz, tzz = traction(m, k, xi, u)
    vx = -mp.mpc(0, 1) * w * u[0]
    vz = -mp.mpc(0, 1) * w * u[1]
    return -mp.mpf(1) / 2 * mp.re(txz * mp.conj(vx) + tzz * mp.conj(vz))


def condition_number(A):
    """2-norm condition number of a complex matrix, via A^H A.

    `mp.cond(A, 2)` returns +inf for complex input (it assumes a real matrix),
    which would flag every case as untrustworthy and the check would be
    worthless. The singular values of A are the square roots of the eigenvalues
    of the Hermitian A^H A, and `mp.eighe` is reliable for complex Hermitian
    input -- unlike mpmath's SVD, which fails on complex matrices (a dead end
    already recorded in the Phase 1C notes).
    """
    n = A.rows
    G = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            G[i, j] = sum(mp.conj(A[m, i]) * A[m, j] for m in range(n))
    Gh = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Gh[i, j] = (G[i, j] + mp.conj(G[j, i])) / 2
    try:
        E, _Q = mp.eighe(Gh)
    except Exception:                                          # noqa: BLE001
        return mp.inf
    vals = [abs(e) for e in E]
    smax = mp.sqrt(max(vals))
    smin = mp.sqrt(min(vals))
    if smin == 0:
        return mp.inf
    return smax / smin


def select_branch(m, k, w, kind, want_sign):
    """Pick the root of the P or S branch whose z-flux has the required sign.

    `want_sign = +1` -> energy travelling toward +z (transmitted, or incident)
    `want_sign = -1` -> energy travelling toward -z (reflected)

    This is a flux test, not a sign-of-Im(xi) test, and the difference matters:
    for an evanescent branch the flux is zero and no branch satisfies the
    condition. In that case BOTH roots are returned flagged as non-transporting
    and the caller must decide, rather than silently receiving a root that does
    not carry the requested flux.
    """
    raw = m.xi_P(k, w) if kind == "P" else m.xi_S(k, w)
    ok, zero = [], []
    for cand in (raw, -raw):
        u = mode_vector(k, cand, kind)
        if u is None:
            continue
        f = flux_z(m, k, cand, u, w)
        if f * want_sign > 0:
            ok.append((cand, u, f))
        elif f == 0:
            zero.append((cand, u, f))
    if ok:
        # if both roots carry the required sign, the branch is not separating
        # incoming from outgoing and the caller has a real problem
        return ok[0][0], ok[0][1], ok[0][2], "propagating", len(ok)
    if zero:
        return (zero[0][0], zero[0][1], zero[0][2], "evanescent", 0)
    return (raw, mode_vector(k, raw, kind), mp.mpf(0), "no_transport", 0)


def fresnel_normal_amplitudes(m1, m2):
    """Closed-form normal-incidence P amplitudes (Rayleigh).

    Displacement convention: incident plus reflected against transmitted.

        r = (Z2 - Z1)/(Z2 + Z1)          t = 2 Z1/(Z2 + Z1)

    These are the expressions the solver must reproduce, and they are the ones
    `code/elastic_psv.py` fails on for t.
    """
    Z1, Z2 = m1.ZL, m2.ZL
    return (Z2 - Z1) / (Z2 + Z1), 2 * Z1 / (Z2 + Z1)


def fresnel_normal_energy(m1, m2):
    r, t = fresnel_normal_amplitudes(m1, m2)
    return abs(r) ** 2, 1 - abs(r) ** 2


# ---------------------------------------------------------------------------
def solve_psv(m1, m2, w, theta_deg, incident="P", welded=True):
    """Welded two-half-space P-SV scattering. Returns ALL complex amplitudes.

    Unknowns, in this order:
        rP, rS   reflected  (medium 1, energy toward -z)
        tP, tS   transmitted (medium 2, energy toward +z)

    The incident mode has UNIT displacement amplitude, in medium 1, travelling
    toward +z. The four welded conditions are continuity of (u_x, u_z, t_xz,
    t_zz); transmitted columns are moved to the left-hand side, which is where
    the Phase 1B sign error came from -- omitting that move negates every
    transmitted amplitude while leaving every energy balance intact.
    """
    if not welded:
        raise NotImplementedError(
            "task 2 scope is welded contact only; the spring-contact path in "
            "code/elastic_psv.py is not part of this reference and was not "
            "validated")

    # theta_deg may be an mpf. Coercing it through float() first -- as an
    # earlier version did -- throws away digits and leaves ~1e-17 of angle
    # error, which lands in the energy residual and is indistinguishable from
    # solver noise. Accept mpf directly.
    theta = theta_deg if isinstance(theta_deg, mp.mpf) else \
        mp.mpf(repr(theta_deg)) if isinstance(theta_deg, float) else mp.mpf(theta_deg)
    k = (w / (m1.cL if incident == "P" else m1.cT)) * mp.sin(mp.radians(theta))

    xi_rP, u_rP, f_rP, st_rP, n_rP = select_branch(m1, k, w, "P", -1)
    xi_rS, u_rS, f_rS, st_rS, n_rS = select_branch(m1, k, w, "S", -1)
    xi_tP, u_tP, f_tP, st_tP, n_tP = select_branch(m2, k, w, "P", +1)
    xi_tS, u_tS, f_tS, st_tS, n_tS = select_branch(m2, k, w, "S", +1)
    xi_in, u_in, f_in, st_in, n_in = select_branch(m1, k, w, incident, +1)

    cols = [(m1, xi_rP, u_rP, +1), (m1, xi_rS, u_rS, +1),
            (m2, xi_tP, u_tP, -1), (m2, xi_tS, u_tS, -1)]

    A = mp.matrix(4, 4)
    b = mp.matrix(4, 1)
    for j, (mm, xi, u, sgn) in enumerate(cols):
        tx, tz = traction(mm, k, xi, u)
        # rows 0,1: u_x, u_z continuity.  rows 2,3: t_xz, t_zz continuity.
        # medium-2 (transmitted) terms carry sgn = -1: they belong on the LHS
        # with a minus sign, and forgetting that negates tP and tS.
        A[0, j] = sgn * u[0]
        A[1, j] = sgn * u[1]
        A[2, j] = sgn * tx
        A[3, j] = sgn * tz

    tx_in, tz_in = traction(m1, k, xi_in, u_in)
    b[0] = -u_in[0]
    b[1] = -u_in[1]
    b[2] = -tx_in
    b[3] = -tz_in

    # CONDITIONING IS PART OF THE ANSWER. A welded 4x4 whose columns are
    # nearly dependent returns amplitudes that are arithmetically exact and
    # physically meaningless. The acoustic limit (C44 -> 0) is exactly such a
    # case: the shear traction goes to zero, the S columns become null
    # directions, and rS is indeterminate -- yet the solve still "succeeds".
    # Recording the condition number is what stops that being reported as a
    # physical reflection coefficient.
    # ROW EQUILIBRATION. The boundary matrix mixes rows whose entries are
    # displacements (~1 m) with rows whose entries are tractions (~1e11 Pa).
    # The raw condition number is then ~1e14 purely because of that unit
    # mismatch, which says nothing about whether the physics is degenerate.
    # Scaling each ROW by its largest entry is exact (it multiplies the equation
    # by a non-zero scalar) and exposes the conditioning that actually matters.
    cond_raw = condition_number(A)
    Aeq = mp.matrix(A)
    beq = mp.matrix(b)
    row_scale = []
    for i in range(4):
        mx = max(abs(Aeq[i, j]) for j in range(4))
        row_scale.append(mx if mx > 0 else mp.mpf(1))
        for j in range(4):
            Aeq[i, j] = Aeq[i, j] / row_scale[-1]
        beq[i] = beq[i] / row_scale[-1]
    cond = condition_number(Aeq)

    sol = mp.lu_solve(Aeq, beq)
    rP, rS, tP, tS = sol[0], sol[1], sol[2], sol[3]

    # ---- INDEPENDENT residuals, recomputed from the amplitudes ---------------
    # The solver's own solve is not evidence. These check the four boundary
    # conditions in the PHYSICAL form (u1 = u2, t1 = t2) rather than in the
    # sign-shifted matrix form, so a column-sign error cannot cancel itself out.
    u1 = (u_in[0] + rP * u_rP[0] + rS * u_rS[0],
          u_in[1] + rP * u_rP[1] + rS * u_rS[1])
    u2 = (tP * u_tP[0] + tS * u_tS[0], tP * u_tP[1] + tS * u_tS[1])
    tx1, tz1 = traction(m1, k, xi_in, u_in)
    t1x = tx1 + rP * traction(m1, k, xi_rP, u_rP)[0] + \
        rS * traction(m1, k, xi_rS, u_rS)[0]
    t1z = tz1 + rP * traction(m1, k, xi_rP, u_rP)[1] + \
        rS * traction(m1, k, xi_rS, u_rS)[1]
    t2x = tP * traction(m2, k, xi_tP, u_tP)[0] + \
        tS * traction(m2, k, xi_tS, u_tS)[0]
    t2z = tP * traction(m2, k, xi_tP, u_tP)[1] + \
        tS * traction(m2, k, xi_tS, u_tS)[1]

    scale_u = max(abs(u1[0]), abs(u1[1]), abs(u2[0]), abs(u2[1]),
                  mp.mpf("1e-300"))
    scale_t = max(abs(t1x), abs(t1z), abs(t2x), abs(t2z), mp.mpf("1e-300"))
    resid = {
        "ux": abs(u1[0] - u2[0]) / scale_u,
        "uz": abs(u1[1] - u2[1]) / scale_u,
        "txz": abs(t1x - t2x) / scale_t,
        "tzz": abs(t1z - t2z) / scale_t,
    }

    # ---- energy -------------------------------------------------------------
    P_in = flux_z(m1, k, xi_in, u_in, w)
    P_rP = flux_z(m1, k, xi_rP, u_rP, w) * abs(rP) ** 2
    P_rS = flux_z(m1, k, xi_rS, u_rS, w) * abs(rS) ** 2
    P_tP = flux_z(m2, k, xi_tP, u_tP, w) * abs(tP) ** 2
    P_tS = flux_z(m2, k, xi_tS, u_tS, w) * abs(tS) ** 2
    if P_in != 0:
        # reflected branches carry flux toward -z, so their contribution to the
        # +z budget is the negative of the raw flux
        R = -(P_rP + P_rS) / P_in
        T = (P_tP + P_tS) / P_in
        bal = abs(R + T - 1)
    else:
        R = T = None
        bal = mp.inf

    return {
        "k": k, "theta_deg": theta_deg, "incident": incident,
        "rP": rP, "rS": rS, "tP": tP, "tS": tS,
        "xi": {"rP": xi_rP, "rS": xi_rS, "tP": xi_tP, "tS": xi_tS,
               "in": xi_in},
        "status": {"rP": st_rP, "rS": st_rS, "tP": st_tP, "tS": st_tS,
                   "in": st_in},
        "residual": resid,
        "cond": cond,
        "cond_raw": cond_raw,
        #: FRACTION OF INCIDENT ENERGY CARRIED BY EACH CHANNEL.
        #: This, not the condition number, is what decides whether an amplitude
        #: is physically meaningful. In the acoustic limit (C44 -> 0) the shear
        #: channel's energy share goes to zero like C44 while |rS| stays at
        #: 0.30 -- the amplitude does not converge, because the S mode becomes a
        #: zero-traction null direction of the (u, t) formulation and its
        #: coefficient is a singular-perturbation residue. Measured: share
        #: 3.2e-8, 3.2e-10, 3.2e-12 for C44 = 1e-2, 1e-6, 1e-10.
        "energy_share": {
            "rP": abs(P_rP / P_in) if P_in else None,
            "rS": abs(P_rS / P_in) if P_in else None,
            "tP": abs(P_tP / P_in) if P_in else None,
            "tS": abs(P_tS / P_in) if P_in else None,
        },
        "P_in": P_in,
        #: below this share a channel is numerically present but physically
        #: negligible, and its amplitude must not be quoted as a coefficient
        "amplitude_negligible": {
            "rS": bool(P_in and abs(P_rS / P_in) < mp.mpf("1e-6")),
            "tS": bool(P_in and abs(P_tS / P_in) < mp.mpf("1e-6")),
        },
        "amplitudes_trustworthy": cond < mp.mpf("1e8"),
        "R": R, "T": T, "balance": bal,
    }


# ---------------------------------------------------------------------------
def _main():
    """Validate this reference against closed forms BEFORE using it."""
    mp.mp.dps = 50
    w = 2 * mp.pi * 1e6

    #: deliberately the same pair Phase 1B used, so the sign defect is visible
    A = ISO("Al6061", 2700.0, 107.6e9, 26.5e9)
    B = ISO("steel-like", 7800.0, 281.0e9, 82.0e9)
    print("Reference validation against CLOSED FORMS\n")
    print(f"  A: {A}")
    print(f"  B: {B}\n")

    r_cf, t_cf = fresnel_normal_amplitudes(A, B)
    d = solve_psv(A, B, w, 0.0)
    print("normal incidence, amplitudes (the Phase 1B defect showed up here)")
    print(f"  closed form  r = {mp.nstr(r_cf, 12)}   t = {mp.nstr(t_cf, 12)}")
    print(f"  this module  r = {mp.nstr(d['rP'], 12)}   t = {mp.nstr(d['tP'], 12)}")
    print(f"  |dr| = {mp.nstr(abs(d['rP'] - r_cf), 3)}   "
          f"|dt| = {mp.nstr(abs(d['tP'] - t_cf), 3)}")
    print(f"  boundary residuals: "
          + ", ".join(f"{k}={mp.nstr(v, 3)}" for k, v in d["residual"].items()))
    print()
    print("oblique incidence -- energy and mode conversion")
    print(f"  {'theta':>6} {'|rP|':>12} {'|rS|':>12} {'|tP|':>12} "
          f"{'|tS|':>12} {'|R+T-1|':>10}")
    for th in (0, 15, 30, 45, 60, 75):
        d = solve_psv(A, B, w, th)
        print(f"  {th:6.1f} {float(abs(d['rP'])):12.8f} "
              f"{float(abs(d['rS'])):12.8f} {float(abs(d['tP'])):12.8f} "
              f"{float(abs(d['tS'])):12.8f} {float(d['balance']):10.2e}")


if __name__ == "__main__":
    _main()
