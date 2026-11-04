"""
plate_spectral.py  (v2)
=======================
Spectral (Chebyshev collocation) solver for guided waves in homogeneous plates,
isotropic elastic AND generally anisotropic piezoelectric.

Design notes / why this v2 exists
---------------------------------
v1 assembled the stiffness blocks by hand.  That is error-prone: the first run
produced an incorrect sign for the piezoelectric stiffening and an incorrect
Rayleigh-Lamb residual (tangent singularities).  v2 therefore assembles the
operator from the *energy bilinear form* directly, using the constitutive
tensors, so that the only inputs are the material tensors themselves:

    A(k) q = omega^2 M(k) q ,
    A = INT [ B_S^T C B_S - B_S^T e B_E - B_E^T e^T B_S - B_E^T eps B_E ] dx2
    M = INT [ rho (N^T N) ] dx2

with S = B_S q and E = B_E q the nodal strain / electric-field operators for a
sagittal-plane (x1, x2) problem, d/dx3 = 0, u3 = 0.

The massless electric potential is eliminated exactly by static condensation
(Guyan reduction), which is equivalent to solving Gauss's law without time
derivatives.

Conventions
-----------
* Time factor exp(-i omega t); propagation factor exp(+i k x1).
* Voigt order for the 6-vectors: (11, 22, 33, 23, 13, 12).
* Engineering shear strains S4 = 2*S23, S5 = 2*S13, S6 = 2*S12.
* SI units.
* For the isotropic case we still use the anisotropic machinery with a
  Lamé-constructed C, so that both paths are validated through one code path.

Boundary conditions
-------------------
* Traction-free surfaces: natural (nothing added).
* Electrically OPEN  (D2 = 0 on both faces): natural.
* Electrically SHORT (phi = 0 on both faces): essential, imposed by deleting
  the potential DOFs of the two boundary nodes.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla

VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


# ---------------------------------------------------------------------------
# Chebyshev
# ---------------------------------------------------------------------------
def cheb(N: int):
    """Chebyshev-Lobatto nodes on [-1,1], differentiation matrix, and
    Clenshaw-Curtis quadrature weights (which integrate polynomials exactly
    up to degree N, sufficient for a spectral method)."""
    if N == 0:
        return np.array([1.0]), np.zeros((1, 1)), np.array([2.0])
    th = np.pi * np.arange(N + 1) / N
    x = np.cos(th)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    # Clenshaw-Curtis weights
    n = np.arange(0, N + 1)
    w = np.zeros(N + 1)
    for j in range(N + 1):
        s = 0.0
        for m in range(1, N // 2 + 1):
            b = 1.0 if 2 * m == N else 2.0
            s += b / (4 * m * m - 1) * np.cos(2 * m * th[j])
        cj = 1.0 if j in (0, N) else 2.0
        w[j] = cj / N * (1.0 - s)
    return x, D, w


# ---------------------------------------------------------------------------
# Material tensors
# ---------------------------------------------------------------------------
def isotropic_tensors(E: float, nu: float, rho: float):
    """Return (C 6x6, e 6x3, eps 3x3, rho) for an isotropic elastic solid."""
    mu = E / (2 * (1 + nu))
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C[2, 2] = lam + 2 * mu
    C[0, 1] = C[0, 2] = C[1, 0] = C[1, 2] = C[2, 0] = C[2, 1] = lam
    C[3, 3] = C[4, 4] = C[5, 5] = mu
    e = np.zeros((6, 3))
    eps = np.eye(3) * 8.8541878128e-12
    return C, e, eps, rho


def crystal_6mm(c11, c12, c13, c33, c44, e15, e31, e33, eps11, eps33):
    """6mm tensors with the POLING AXIS along x3 (crystal frame)."""
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = c11
    C[2, 2] = c33
    C[0, 1] = C[1, 0] = c12
    C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = c13
    C[3, 3] = C[4, 4] = c44
    C[5, 5] = 0.5 * (c11 - c12)
    e = np.zeros((6, 3))
    e[0, 2] = e[1, 2] = e31      # T1, T2 <- E3
    e[2, 2] = e33                # T3 <- E3
    e[3, 1] = e[4, 0] = e15      # T4<-E2, T5<-E1
    eps = np.diag([eps11, eps11, eps33])
    return C, e, eps


def bond_rotate(C, e, eps, R):
    """Rotate stiffness / piezoelectric / dielectric tensors by rotation matrix R
    (columns = new basis vectors expressed in the old frame).

    Uses the standard Bond-matrix transformation, with the piezoelectric tensor
    rotated as a third-rank tensor via e_new = R e R^T in (6x3) reduced form.
    Verified in check_02 by invariance tests (isotropic C must be unchanged).
    """
    # Build full 4th-order C from Voigt
    C4 = np.zeros((3, 3, 3, 3))
    for I, (i, j) in enumerate(VOIGT):
        for J, (k, l) in enumerate(VOIGT):
            v = C[I, J]
            C4[i, j, k, l] = v
            C4[j, i, k, l] = v
            C4[i, j, l, k] = v
            C4[j, i, l, k] = v
    C4n = np.einsum("pi,qj,rk,sl,pqrs->ijkl", R, R, R, R, C4, optimize=True)

    # 3rd-rank piezoelectric e_{kij} with Voigt (k, (ij))
    e3 = np.zeros((3, 3, 3))
    for J, (i, j) in enumerate(VOIGT):
        for k in range(3):
            e3[k, i, j] = e[J, k]
            e3[k, j, i] = e[J, k]
    e3n = np.einsum("pi,qj,rk,pqr->ijk", R, R, R, e3, optimize=True)

    epsn = R @ eps @ R.T

    # back to Voigt
    Cn = np.zeros((6, 6))
    for I, (i, j) in enumerate(VOIGT):
        for J, (k, l) in enumerate(VOIGT):
            Cn[I, J] = C4n[i, j, k, l]
    en = np.zeros((6, 3))
    for J, (i, j) in enumerate(VOIGT):
        for k in range(3):
            en[J, k] = e3n[k, i, j]
    return Cn, en, epsn


# ---------------------------------------------------------------------------
# Plate layer
# ---------------------------------------------------------------------------
class Plate:
    """Homogeneous plate occupying -h <= x2 <= h, sagittal motion (u1,u2,phi)."""

    def __init__(self, h: float, C, e, eps, rho: float, N: int = 28,
                 electrical: str = "open"):
        assert electrical in ("open", "short")
        self.h, self.C, self.e, self.eps, self.rho = h, np.asarray(C), np.asarray(e), np.asarray(eps), rho
        self.N, self.electrical = N, electrical
        x, Dn, w = cheb(N)
        self.x2 = h * x
        self.D = Dn / h
        self.w = w * h
        n = N + 1
        self.n = n
        self.nu = 2 * n                    # mechanical dofs
        self.nphi = n                      # potential dofs (before BC)
        self.ndof = 3 * n

    # -- operators ---------------------------------------------------------
    def _BE(self, k: float):
        """Nodal strain / electric-field operators.
        Returns BS with shape (6, n, ndof) and BE with shape (3, n, ndof),
        so that S_I = BS[I] @ q and E_m = BE[m] @ q.
        (A (6, ndof) layout cannot represent derivative operators correctly --
        each strain component needs its own n x ndof operator.)"""
        n, D = self.n, self.D
        nd = self.ndof
        ik = 1j * k
        I = np.eye(n, dtype=complex)
        BS = np.zeros((6, n, nd), dtype=complex)
        BE = np.zeros((3, n, nd), dtype=complex)
        # S1 = ik u1
        BS[0][:, 0:n] = ik * I
        # S2 = du2/dx2
        BS[1][:, n:2 * n] = D
        # S3 = S4 = S5 = 0   (u3 = 0, d/dx3 = 0)
        # S6 = 2 S12 = du1/dx2 + ik u2
        BS[5][:, 0:n] = D
        BS[5][:, n:2 * n] = ik * I
        # E1 = -dphi/dx1 = -ik phi
        BE[0][:, 2 * n:3 * n] = -ik * I
        # E2 = -dphi/dx2
        BE[1][:, 2 * n:3 * n] = -D
        return BS, BE

    def assemble(self, k: float, condense: bool = True):
        """Return (A, M) for A q = omega^2 M q, with the massless potential
        statically condensed out."""
        BS, BE = self._BE(k)
        C, e, eps = self.C, self.e, self.eps
        w = self.w
        n, nd = self.n, self.ndof

        def wint(X, Y):
            """integral of X^T Y over the thickness (X, Y are n x nd)."""
            return X.conj().T @ (w[:, None] * Y)

        A = np.zeros((nd, nd), dtype=complex)
        # elastic:  sum_{I,J} BS_I^T W C_IJ BS_J
        for I in range(6):
            for J in range(6):
                if C[I, J] != 0.0:
                    A += C[I, J] * wint(BS[I], BS[J])
        # piezoelectric:  - sum_{I,m} e_Im ( BS_I^T W BE_m + BE_m^T W BS_I )
        for I in range(6):
            for m in range(3):
                if e[I, m] != 0.0:
                    A -= e[I, m] * wint(BS[I], BE[m])
                    A -= e[I, m] * wint(BE[m], BS[I])
        # dielectric:  - sum_{m,l} eps_ml BE_m^T W BE_l
        for m in range(3):
            for l in range(3):
                if eps[m, l] != 0.0:
                    A -= eps[m, l] * wint(BE[m], BE[l])
        A = 0.5 * (A + A.conj().T)          # enforce Hermiticity

        M = np.zeros((nd, nd))
        rn = np.diag(self.rho * w)
        M[0:n, 0:n] = rn
        M[n:2 * n, n:2 * n] = rn

        # --- electrical boundary condition ---
        if self.electrical == "short":
            kill = np.array([2 * n, 3 * n - 1])      # phi at both faces = 0
            keep = np.setdiff1d(np.arange(nd), kill)
            A, M = A[np.ix_(keep, keep)], M[np.ix_(keep, keep)]
        if condense:
            nn = 2 * n
            A_uu, A_up = A[:nn, :nn], A[:nn, nn:]
            A_pu, A_pp = A[nn:, :nn], A[nn:, nn:]
            if A_pp.size:
                A = A_uu - A_up @ np.linalg.pinv(A_pp) @ A_pu
            else:
                A = A_uu
            M = M[:nn, :nn]
        return A, M

    # -- spectrum ----------------------------------------------------------
    def spectrum(self, k: float):
        """Sorted real frequencies f (Hz) for real k. Negative/zero dropped."""
        A, M = self.assemble(k)
        vals = sla.eig(A, M, right=False)
        vals = np.asarray(vals, dtype=complex)
        vals = vals[np.isfinite(vals.real) & np.isfinite(vals.imag)]
        real = np.abs(vals.imag) < 1e-7 * (np.abs(vals.real) + 1.0)
        vals = vals[real].real
        vals = vals[vals > 1e-3]
        f = np.sqrt(vals) / (2 * np.pi)
        return np.sort(f)

    def eigen(self, k: float, f_target: float):
        """Eigenvector of the *condensed* problem closest to f_target, expanded
        back to (u1, u2) nodal values."""
        A, M = self.assemble(k)
        vals, vecs = sla.eig(A, M)
        vals = np.asarray(vals, dtype=complex)
        tgt = (2 * np.pi * f_target) ** 2
        idx = int(np.nanargmin(np.abs(vals - tgt)))
        q = vecs[:, idx]
        q = q / np.max(np.abs(q))
        return q, np.sqrt(vals[idx].real) / (2 * np.pi)


# ---------------------------------------------------------------------------
# Exact isotropic theory (for validation)
# ---------------------------------------------------------------------------
def rayleigh_lamb_entire(omega, k, cL, cT, h, sym=True):
    """Entire (regularised) Rayleigh-Lamb dispersion function for a plate of
    half-thickness h.  Multiplying the classical equation by cos(p h) cos(q h)
    removes the tangent singularities, so the function is smooth and its zeros
    are exactly the guided modes.  Convention exp(i(kx - omega t)).

    Symmetric:     Q (q^2 - k^2)^2 + P * 4 k^2 p q = 0
    Antisymmetric: P (q^2 - k^2)^2 + Q * 4 k^2 p q = 0
    with P = sin(p h) cos(q h),  Q = sin(q h) cos(p h).
    (The two cases swap P and Q -- an earlier version of this helper had the
    antisymmetric case wrong, which is why the validation is run both ways.)
    """
    omega = np.asarray(omega, dtype=complex)
    p = np.sqrt((omega / cL) ** 2 - k ** 2 + 0j)
    q = np.sqrt((omega / cT) ** 2 - k ** 2 + 0j)
    P = np.sin(p * h) * np.cos(q * h)
    Q = np.sin(q * h) * np.cos(p * h)
    A = (q ** 2 - k ** 2) ** 2
    B = 4 * k ** 2 * p * q
    return Q * A + P * B if sym else P * A + Q * B


def exact_cutoffs(cL, cT, h, nmax=8):
    """Exact cutoff frequencies (Hz) at k = 0 for a plate of half-thickness h:
    f_n = n c / (4h) for both the longitudinal and transverse families."""
    vals = set()
    for n in range(nmax):
        vals.add(n * cT / (4 * h))
        vals.add(n * cL / (4 * h))
    return np.array(sorted(v for v in vals if v > 0))


def rl_roots(k, cL, cT, h, f_lo, f_hi, n=40000, tol=1e-2):
    """All Rayleigh-Lamb roots in [f_lo, f_hi] at fixed k, refined by bisection.

    Branch points of the vertical slowness  p = sqrt((w/cL)^2 - k^2)  and
    q = sqrt((w/cT)^2 - k^2)  occur at  f = cL k / 2pi  and  f = cT k / 2pi .
    There the regularised dispersion function is continuous in magnitude but the
    sign change of the square root produces a spurious sign flip, so a naive
    sign-change root finder reports a false root.  We therefore discard roots
    lying within `tol` (relative) of either branch point.

    Returns (f_sym, f_anti) in Hz.
    """
    fs = np.linspace(f_lo, f_hi, n)
    W = 2 * np.pi * fs
    bp = np.array([cL * k / (2 * np.pi), cT * k / (2 * np.pi)])
    out = []
    for sym in (True, False):
        r = rayleigh_lamb_entire(W, k, cL, cT, h, sym=sym).real
        sgn = np.sign(r)
        idx = np.where(np.diff(sgn) != 0)[0]
        roots = []
        for i in idx:
            a, b = fs[i], fs[i + 1]
            fa = r[i]
            ok = True
            for _ in range(90):
                m = 0.5 * (a + b)
                fm = float(rayleigh_lamb_entire(2 * np.pi * m, k, cL, cT, h, sym=sym).real)
                if np.sign(fm) == np.sign(fa):
                    a, fa = m, fm
                else:
                    b = m
            f = 0.5 * (a + b)
            if np.any(np.abs(f - bp) / f < tol):
                ok = False           # branch-point artefact
            if ok:
                roots.append(f)
        out.append(np.array(sorted(roots)))
    return out[0], out[1]


def energy_velocity(pl, k: float, f: float):
    """Energy (group) velocity of a guided mode from the flux identity

        v_e = P / (2 <T>)  ,

    where P is the time-averaged power flux along x1 and <T> the time-averaged
    kinetic energy per unit length.  For a lossless, non-dissipative mode this
    equals d(omega)/dk, but it is obtained WITHOUT differentiating a dispersion
    curve, so it needs no branch tracking and is robust across branch crossings.

    Returns (v_e, P, T).
    """
    q, _ = pl.eigen(k, f)
    n, D, w = pl.n, pl.D, pl.w
    u1, u2 = q[0:n], q[n:2 * n]
    a = np.max(np.abs(q)) + 1e-300
    u1, u2 = u1 / a, u2 / a
    ik = 1j * k
    Sv = np.zeros((6, n), dtype=complex)
    Sv[0] = ik * u1
    Sv[1] = D @ u2
    Sv[5] = D @ u1 + ik * u2
    sig = pl.C @ Sv
    om = 2 * np.pi * f
    P = (om / 2.0) * np.imag(np.sum(w * (sig[0] * np.conj(u1) + sig[5] * np.conj(u2))))
    T = (pl.rho * om ** 2 / 4.0) * np.sum(w * (np.abs(u1) ** 2 + np.abs(u2) ** 2))
    ve = P / (2 * T) if T > 0 else np.nan
    return ve, P, T


def eigen_full(pl, k: float, f: float):
    """Return (u, phi) nodal arrays for the mode nearest f.

    The potential is recovered from Gauss's law with no time derivative,
        A_{phi phi} phi = - A_{phi u} u ,
    which is the same static-condensation relation used to build the reduced
    problem.  This is needed to inspect how uniform the potential is through the
    thickness -- the reason the open/short surface condition hardly matters.
    """
    BS, BE = pl._BE(k)
    n, nd = pl.n, pl.ndof
    C, e, eps = pl.C, pl.e, pl.eps

    def wint(X, Y):
        return X.conj().T @ (pl.w[:, None] * Y)

    # full matrices (no condensation)
    A = np.zeros((nd, nd), dtype=complex)
    for I in range(6):
        for J in range(6):
            if C[I, J] != 0.0:
                A += C[I, J] * wint(BS[I], BS[J])
    for I in range(6):
        for m in range(3):
            if e[I, m] != 0.0:
                A -= e[I, m] * wint(BS[I], BE[m])
                A -= e[I, m] * wint(BE[m], BS[I])
    for m in range(3):
        for l in range(3):
            if eps[m, l] != 0.0:
                A -= eps[m, l] * wint(BE[m], BE[l])
    A = 0.5 * (A + A.conj().T)
    if pl.electrical == "short":
        kill = np.array([2 * n, 3 * n - 1])
        keep = np.setdiff1d(np.arange(nd), kill)
        A = A[np.ix_(keep, keep)]

    Ared, Mred = pl.assemble(k)          # condensed problem
    vals, vecs = sla.eig(Ared, Mred)
    vals = np.asarray(vals, dtype=complex)
    tgt = (2 * np.pi * f) ** 2
    idx = int(np.nanargmin(np.abs(vals - tgt)))
    u = vecs[:, idx]

    # assemble the uncoupled potential blocks
    phis = slice(2 * n, 3 * n)
    App = np.zeros((n, n), complex)
    for m in range(3):
        for l in range(3):
            if eps[m, l] != 0.0:
                App -= eps[m, l] * wint(BE[m][:, phis], BE[l][:, phis])
    Aup = np.zeros((2 * n, n), complex)
    for I in range(6):
        for m in range(3):
            if e[I, m] != 0.0:
                Aup -= e[I, m] * wint(BS[I][:, :2 * n], BE[m][:, phis])
    phi = -np.linalg.pinv(App) @ (Aup.conj().T @ u)
    return u, phi
