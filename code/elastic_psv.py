"""
code/elastic_psv.py  --  Phase 1B, supervisor P1 item.

A REDUCED elastic P-SV scattering solver, written independently of the
piezothermoelastic machinery in foundation_checks.py.

Why this file exists
--------------------
Phase 1 claimed that the decoupled (e = beta = p = 0) piezothermoelastic problem
is singular and therefore "must be regularised by taking coupling -> 0".  The
supervisor's instruction is that this must not be assumed: the correct procedure
is to compare any such limit against a PROPERLY REDUCED elastic subsystem, in
which the electric and thermal variables are REMOVED from the formulation rather
than left in as zero rows.

That distinction matters.  In the full 8x8 coupled system, setting e = beta = p
= 0 exactly leaves phi and theta appearing in no other equation, so four of the
eight interface conditions become identically zero rows and the matrix is
singular *as a formulation artefact*.  In the reduced 4x4 elastic system the
same physics is non-singular and the classical Fresnel and P-SV results hold
exactly.  The reduced system is therefore the correct reference, and it is the
test the coupled solver must be compared against -- not the other way round.

Sign and convention choices (stated so they can be checked independently):
  * time dependence  e^{-i w t}
  * plane wave       exp(i(k x + xi z - w t))
  * traction         t_zz = sigma_zz, t_xz = sigma_xz
  * power flux       <P_z> = -1/2 Re(sigma_xz v_x^* + sigma_zz v_z^*)
  * dimensionless   -> all energies are per unit incident flux
"""

import numpy as np


class ElasticISO:
    """Isotropic elastic medium (the sagittally-isotropic reduction)."""

    def __init__(self, name, rho, C11, C44):
        self.name = name
        self.rho = float(rho)
        self.C11 = float(C11)
        self.C44 = float(C44)
        self.cL = np.sqrt(self.C11 / self.rho)
        self.cT = np.sqrt(self.C44 / self.rho)
        self.ZL = self.rho * self.cL          # P impedance
        self.ZT = self.rho * self.cT          # S impedance

    def xi_L(self, k, w):
        return np.sqrt((w / self.cL) ** 2 - k ** 2 + 0j)

    def xi_T(self, k, w):
        return np.sqrt((w / self.cT) ** 2 - k ** 2 + 0j)


def mode_vectors(k, xi, kind):
    """Displacement (u_x, u_z) of a plane wave with wavenumbers (k, xi).

    kind = 'P': polarisation along the wavevector (k, xi)
    kind = 'S': polarisation perpendicular to it, in the sagittal plane
    The returned vectors are normalised to unit magnitude, which is legitimate
    here because the elastic system has ONE physical dimension (length) and the
    two components share it -- unlike the coupled case, where u, phi and theta
    carry different units.
    """
    if kind == "P":
        v = np.array([k, xi], dtype=complex)
    elif kind == "S":
        v = np.array([xi, -k], dtype=complex)
    else:
        raise ValueError(kind)
    n = np.sqrt(np.sum(v * np.conj(v)))
    if n == 0:
        return None
    return v / n


def tractions(m, k, xi, u):
    """(t_xz, t_zz) for displacement amplitude u = (u_x, u_z) at (k, xi).

    sigma_zz = C13 e_xx + C33 e_zz with C33 = C11, C13 = C11 - 2 C44, and
    sigma_xz = C44 gamma_xz.  Plane strain: e_yy = 0.
    """
    ux, uz = u
    exx = 1j * k * ux
    ezz = 1j * xi * uz
    gxz = 1j * (xi * ux + k * uz)
    C12 = m.C11 - 2 * m.C44
    tzz = C12 * exx + m.C11 * ezz
    txz = m.C44 * gxz
    return txz, tzz


def flux_z(m, k, xi, u, w):
    """Cycle-averaged z-flux of one mode, per unit amplitude squared."""
    txz, tzz = tractions(m, k, xi, u)
    vx = -1j * w * u[0]
    vz = -1j * w * u[1]
    return -0.5 * np.real(txz * np.conj(vx) + tzz * np.conj(vz))


def solve_psv(m1, m2, w, theta_deg, incident="P", welded=True, K=0.0, C=0.0):
    """Welded or spring-contact P-SV scattering between two isotropic media.

    Returns a dict with the amplitudes and the energy account.  Unknowns:
    r_P, r_S (medium 1) and t_P, t_S (medium 2); the incident mode is P (or S)
    with unit displacement amplitude.

    Interface conditions (welded):
        u_x, u_z, t_xz, t_zz continuous   ->  4 equations, 4 unknowns
    Spring contact with scalar stiffness K and damping C (normal + shear share
    the same diagonal pair here, since the medium is isotropic):
        t continuous, and  t = (K + i w C) [[u]] with [[u]] = u2 - u1.
    """
    k = (w / (m1.cL if incident == "P" else m1.cT)) * np.sin(
        np.deg2rad(theta_deg))

    # medium 1: incident + reflected.  Outgoing (reflected) means carrying
    # energy in -z, i.e. Im(xi) >= 0 with the exp(+i xi z) convention AND a
    # negative real part for propagating waves.
    xiP1, xiS1 = m1.xi_L(k, w), m1.xi_T(k, w)
    # medium 2: transmitted, carrying energy in +z
    xiP2, xiS2 = m2.xi_L(k, w), m2.xi_T(k, w)

    def pick_out(xi, sign):
        """Choose the branch whose flux has the required sign."""
        for cand in (xi, -xi):
            u = mode_vectors(k, cand, "P")
            if u is None:
                continue
            f = flux_z(m1, k, cand, u, w)
            if sign * f > 0:
                return cand
        return xi if sign > 0 else -xi

    def pick_out_typed(mm, xi, kind, sign):
        for cand in (xi, -xi):
            u = mode_vectors(k, cand, kind)
            if u is None:
                continue
            f = flux_z(mm, k, cand, u, w)
            if sign * f > 0:
                return cand
        return xi

    xiP1 = pick_out_typed(m1, xiP1, "P", -1)     # reflected: flux in -z
    xiS1 = pick_out_typed(m1, xiS1, "S", -1)
    xiP2 = pick_out_typed(m2, xiP2, "P", +1)     # transmitted: flux in +z
    xiS2 = pick_out_typed(m2, xiS2, "S", +1)

    uP1 = mode_vectors(k, xiP1, "P")
    uS1 = mode_vectors(k, xiS1, "S")
    uP2 = mode_vectors(k, xiP2, "P")
    uS2 = mode_vectors(k, xiS2, "S")

    if incident == "P":
        u_in = mode_vectors(k, -xiP1 if False else
                            pick_out_typed(m1, m1.xi_L(k, w), "P", +1), "P")
        xi_in = pick_out_typed(m1, m1.xi_L(k, w), "P", +1)
    else:
        xi_in = pick_out_typed(m1, m1.xi_T(k, w), "S", +1)
        u_in = mode_vectors(k, xi_in, "S")

    # 8x8-free assembly: 4 equations for the 4 unknowns.
    # Column order: r_P, r_S, t_P, t_S
    cols_u = [uP1, uS1, uP2, uS2]
    cols_x = [xiP1, xiS1, xiP2, xiS2]

    A = np.zeros((4, 4), dtype=complex)
    b = np.zeros(4, dtype=complex)

    tx_in, tz_in = tractions(m1, k, xi_in, u_in)

    for j, (u, xi) in enumerate(zip(cols_u, cols_x)):
        mm = m1 if j < 2 else m2
        tx, tz = tractions(mm, k, xi, u)
        A[0, j] = u[0]
        A[1, j] = u[1]
        A[2, j] = tx
        A[3, j] = tz
    b[0] = -u_in[0]
    b[1] = -u_in[1]
    b[2] = -tx_in
    b[3] = -tz_in

    if not welded:
        # spring contact: replace the two traction-continuity rows by
        #   t1 - t2 = 0        (still continuity of traction, given t = K[[u]])
        # is NOT what a spring does.  Correct: t_c = K'[[u]] with t_c the common
        # traction, so  traction_1 = traction_2  AND  traction = K' (u2 - u1).
        # That is 4 equations again; build directly.
        kc = complex(K, w * C) if False else (K + 1j * w * C)
        A = np.zeros((4, 4), dtype=complex)
        # rows 0,1 : t2 - K'(u2-u1) = 0   (using medium-2 tractions)
        # rows 2,3 : t1 - t2 = 0
        for j, (u, xi) in enumerate(zip(cols_u, cols_x)):
            mm = m1 if j < 2 else m2
            tx, tz = tractions(mm, k, xi, u)
            if j >= 2:
                A[0, j] = tx - kc * u[0]
                A[1, j] = tz - kc * u[1]
            else:
                A[0, j] = -(-kc * u[0])
                A[1, j] = -(-kc * u[1])
            A[2, j] = tx * (1 if j < 2 else -1)
            A[3, j] = tz * (1 if j < 2 else -1)
        b[0] = -(tx_in - kc * u_in[0])
        b[1] = -(tz_in - kc * u_in[1])
        b[2] = -tx_in
        b[3] = -tz_in

    sol = np.linalg.solve(A, b)
    rP, rS, tP, tS = sol

    # ---- energy account, from the flux of each channel -------------------
    P_in = flux_z(m1, k, xi_in, u_in, w)

    def chan_flux(mm, xi, u, amp):
        return flux_z(mm, k, xi, u, w) * abs(amp) ** 2

    P_rP = chan_flux(m1, xiP1, uP1, rP)
    P_rS = chan_flux(m1, xiS1, uS1, rS)
    P_tP = chan_flux(m2, xiP2, uP2, tP)
    P_tS = chan_flux(m2, xiS2, uS2, tS)
    # reflected channels carry flux in -z, so R is MINUS their sum
    R = -(P_rP + P_rS) / P_in
    T = (P_tP + P_tS) / P_in
    return dict(R=float(np.real(R)), T=float(np.real(T)),
                rP=rP, rS=rS, tP=tP, tS=tS, k=k, theta=theta_deg,
                P_in=P_in, residual=float(abs(R + T - 1)))


def fresnel_normal(m1, m2):
    """Exact normal-incidence P reflection between two isotropic media."""
    Z1 = m1.ZL
    Z2 = m2.ZL
    r = (Z2 - Z1) / (Z2 + Z1)
    return abs(r) ** 2, 1 - abs(r) ** 2
