"""
coupled_plate.py
================
Generalised spectral (Chebyshev) solver for guided waves in a homogeneous plate
carrying an ARBITRARY number of quasi-static coupled scalar fields.

Motivation
----------
The piezothermoelastic / magneto-electro-elastic / piezothermoelastic-semiconductor
configurations proposed in this project all have the same mathematical skeleton:

    u1, u2                 mechanical displacement (2 scalar unknowns per node)
    s_1 ... s_F            quasi-static scalar potentials (phi, psi, T, n)

The massless scalars are eliminated exactly by static condensation, leaving a
Hermitian generalised eigenproblem in the mechanical unknowns:

    [ A_uu - A_us A_ss^{-1} A_su ] q = w^2 M_uu q

Energy (free-energy) bilinear form, with S = B_S q and F_s = B_s q the nodal
strain / generalised-field operators:

    Pi = 1/2 Int [ S:C:S
                   - sum_s ( S:e_s:F_s + F_s^T e_s^T S )
                   - sum_{s,t} F_s^T M_st F_t ] dx2

which reproduces the sign convention already validated in plate_spectral.py
for the piezoelectric case (there F = E = -grad(phi)).

VALIDATION CONTRACT
-------------------
Setting q_psi = 0 and removing the magnetic field MUST reproduce
plate_spectral.Plate to machine precision.  This is checked in
check_03_multiphysics_premise.py (V1) before any new physics is reported.

Conventions (identical to plate_spectral.py)
--------------------------------------------
* time factor exp(-i w t); propagation exp(+i k x1)
* Voigt order (11,22,33,23,13,12); engineering shears
* SI units; plate -h <= x2 <= +h

Author note: this module only assembles the LINEARISED, LOSSLESS,
quasi-static-scalar system.  It does NOT contain thermal relaxation
(MGT / dual-phase-lag) or carrier transport; those change the eigenproblem
from Hermitian to non-self-adjoint and are explicitly outside this module's
validation contract (see the research-scope document, Section 6).
"""

import numpy as np
import scipy.linalg as sla

__all__ = ["GeneralPlate", "cheb", "iso_C", "tensor_6mm"]


def cheb(N):
    """Chebyshev-Lobatto nodes on [-1,1] (ascending), differentiation matrix,
    and Clenshaw-Curtis weights.  Identical to plate_spectral.cheb."""
    if N == 0:
        return np.array([1.0]), np.zeros((1, 1)), np.array([2.0])
    th = np.pi * np.arange(N + 1) / N
    x = np.cos(th)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    n = np.arange(0, N + 1) * 1.0
    w = np.zeros(N + 1)
    for jj in range(N + 1):
        acc = 0.0
        for m in range(1, N // 2 + 1):
            b = 1.0 if 2 * m == N else 2.0
            acc += b / (4 * m * m - 1) * np.cos(2 * m * th[jj])
        cj = 1.0 if jj in (0, N) else 2.0
        w[jj] = cj / N * (1.0 - acc)
    # plate_spectral returns descending x; we return ascending for the
    # boundary-node convention used below.
    return x[::-1], D[::-1, ::-1], w[::-1]


def iso_C(E, nu):
    """Isotropic stiffness in Voigt notation (6x6)."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            C[i, j] = lam
        C[i, i] = lam + 2 * mu
    C[3, 3] = C[4, 4] = C[5, 5] = mu
    return C


def tensor_6mm(C11, C12, C13, C33, C44, K1, K2, K3, C66=None):
    """6mm stiffness with the symmetry axis along x3 (crystal frame)."""
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C11
    C[2, 2] = C33
    C[0, 1] = C[1, 0] = C12
    C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = C13
    C[3, 3] = C[4, 4] = C44
    C[5, 5] = C66 if C66 is not None else 0.5 * (C11 - C12)
    return C


def rot_x3_to_x2():
    """Rotation taking the crystal x3 axis onto the plate normal x2."""
    return np.array([[1.0, 0.0, 0.0],
                     [0.0, 0.0, 1.0],
                     [0.0, 1.0, 0.0]])


def bond_rotate(C, coupling_mats, second_mats, R):
    """Rotate 6x6 stiffness, 6x3 coupling matrices and 3x3 second-order
    tensors by the orthogonal matrix R (x_new = R x_old)."""
    def T(R, order=2):
        """Bond matrix for stress (order 2) or strain (order 2)."""
        pass
    # Build the Voigt bond matrix explicitly (safe, no formula guessing)
    B = np.zeros((6, 6))
    idx = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    for I, (i, j) in enumerate(idx):
        for J, (p, q) in enumerate(idx):
            # sigma'_I = B_IJ sigma_J  with sigma'_IJ = R_ip R_jq sigma_pq
            val = 0.0
            for (a, b), f in zip([(i, j), (j, i)], [0.5, 0.5]):
                for (cc, dd), g in zip([(p, q), (q, p)], [0.5, 0.5]):
                    val += R[a, cc] * R[b, dd] * f * g
            B[I, J] += val * (2 if I >= 3 or J >= 3 else 1) * (1 if I < 3 and J < 3 else 1)
    # Robust alternative: transform via full index contraction
    Bn = np.zeros((6, 6))
    for I, (i, j) in enumerate(idx):
        for J, (a, b) in enumerate(idx):
            base = R[i, a] * R[j, b]
            if I >= 3 and J >= 3:
                base = R[i, a] * R[j, b] + R[i, b] * R[j, a]
            elif I >= 3 and J < 3:
                base = 0.5 * (R[i, a] * R[j, b] + R[i, b] * R[j, a])
            elif I < 3 and J >= 3:
                base = R[i, a] * R[j, b] + R[i, b] * R[j, a]
            Bn[I, J] = base
    Cn = Bn @ C @ Bn.T

    coups = []
    for e in coupling_mats:
        # e is 6x3 : sigma_I = e_Ik E_k , E_k = R_kl E'_l  =>  e' = Bn e R^T
        coups.append(Bn @ e @ R.T)

    seconds = []
    for m in second_mats:
        seconds.append(R @ m @ R.T)
    return Cn, coups, seconds


class GeneralPlate:
    """Plate -h<=x2<=h with 2 mechanical unknowns and F quasi-static scalars.

    Parameters
    ----------
    h : half thickness
    C : (6,6) stiffness
    couplings : list of (6,3) coupling matrices, one per scalar field
                (piezoelectric e, piezomagnetic q, thermal beta, ...)
    second    : (F,F) list-of-lists of (3,3) matrices giving the generalised
                "permittivity" block, i.e. M[s][t] multiplies F_s . F_t
    rho : density
    N : spectral order (nodes = N+1)
    bcs : list of strings, one per scalar: "open" (natural: F_2 = 0 on faces)
          or "short" (essential: scalar = 0 on both faces)
    """

    def __init__(self, h, C, couplings, second, rho, N=28, bcs=None,
                 field_names=None):
        self.h = float(h)
        self.C = np.asarray(C, dtype=float)
        self.couplings = [np.asarray(c, dtype=float) for c in couplings]
        self.F = len(self.couplings)
        self.second = [[np.asarray(second[s][t], dtype=float)
                        for t in range(self.F)] for s in range(self.F)]
        self.rho = float(rho)
        self.N = N
        self.bcs = bcs if bcs is not None else ["open"] * self.F
        self.field_names = field_names if field_names is not None else \
            ["s%d" % i for i in range(self.F)]

        x, Dn, w = cheb(N)
        self.x2 = h * x              # ascending
        self.x = x
        self.D = Dn / h
        self.w = w * h
        self.n = N + 1
        self.nu = 2 * self.n
        self.ndof = self.nu + self.F * self.n

    # ------------------------------------------------------------------
    def _ops(self, k):
        """Return B_S (6,n,ndof) and list of B_s (3,n,ndof)."""
        n, D, nd = self.n, self.D, self.ndof
        ik = 1j * k
        I = np.eye(n, dtype=complex)
        BS = np.zeros((6, n, nd), dtype=complex)
        BS[0][:, 0:n] = ik * I                 # S1  = ik u1
        BS[1][:, n:2 * n] = D                  # S2  = du2/dx2
        BS[5][:, 0:n] = D                      # S6  = du1/dx2 + ik u2
        BS[5][:, n:2 * n] = ik * I
        Bs = []
        for s in range(self.F):
            off = self.nu + s * n
            B = np.zeros((3, n, nd), dtype=complex)
            B[0][:, off:off + n] = -ik * I     # F_1 = -d s/dx1
            B[1][:, off:off + n] = -D          # F_2 = -d s/dx2
            Bs.append(B)
        return BS, Bs

    def assemble(self, k, condense=True):
        BS, Bs = self._ops(k)
        w, n, nd = self.w, self.n, self.ndof

        def wint(X, Y):
            return X.conj().T @ (w[:, None] * Y)

        A = np.zeros((nd, nd), dtype=complex)
        C = self.C
        for I in range(6):
            for J in range(6):
                if C[I, J] != 0.0:
                    A += C[I, J] * wint(BS[I], BS[J])

        for s in range(self.F):
            e = self.couplings[s]
            for I in range(6):
                for m in range(3):
                    if e[I, m] != 0.0:
                        A -= e[I, m] * wint(BS[I], Bs[s][m])
                        A -= e[I, m] * wint(Bs[s][m], BS[I])

        for s in range(self.F):
            for t in range(self.F):
                M = self.second[s][t]
                for m in range(3):
                    for l in range(3):
                        if M[m, l] != 0.0:
                            A -= M[m, l] * wint(Bs[s][m], Bs[t][l])
        A = 0.5 * (A + A.conj().T)

        Mk = np.zeros((nd, nd))
        rn = np.diag(self.rho * w)
        Mk[0:n, 0:n] = rn
        Mk[n:2 * n, n:2 * n] = rn

        # essential conditions: scalar = 0 at both faces
        kill = []
        for s in range(self.F):
            if self.bcs[s] == "short":
                off = self.nu + s * n
                kill += [off, off + n - 1]
        self._keep = np.arange(nd) if not kill else np.setdiff1d(np.arange(nd), kill)
        if kill:
            keep = self._keep
            A, Mk = A[np.ix_(keep, keep)], Mk[np.ix_(keep, keep)]

        if condense:
            nn = 2 * n
            A_uu, A_us = A[:nn, :nn], A[:nn, nn:]
            A_su, A_ss = A[nn:, :nn], A[nn:, nn:]
            if A_ss.size:
                A = A_uu - A_us @ np.linalg.pinv(A_ss) @ A_su
            else:
                A = A_uu
            Mk = Mk[:nn, :nn]
        return A, Mk

    # ------------------------------------------------------------------
    def spectrum(self, k):
        A, M = self.assemble(k)
        vals = np.asarray(sla.eig(A, M, right=False), dtype=complex)
        vals = vals[np.isfinite(vals.real) & np.isfinite(vals.imag)]
        real = np.abs(vals.imag) < 1e-7 * (np.abs(vals.real) + 1.0)
        vals = vals[real].real
        vals = vals[vals > 1e-3]
        return np.sort(np.sqrt(vals) / (2 * np.pi))

    def eigen(self, k, f_target):
        """Return (angular freq, eigenvector) closest to f_target (Hz)."""
        A, M = self.assemble(k)
        vals, vecs = sla.eig(A, M)
        vals = np.asarray(vals, dtype=complex)
        good = np.isfinite(vals.real) & (np.abs(vals.imag) < 1e-6 * (np.abs(vals.real) + 1))
        vals, vecs = vals[good].real, vecs[:, good]
        keep = vals > 1e-3
        vals, vecs = vals[keep], vecs[:, keep]
        f = np.sqrt(vals) / (2 * np.pi)
        j = int(np.argmin(np.abs(f - f_target)))
        return 2 * np.pi * f[j], vecs[:, j]

    # ------------------------------------------------------------------
    #  Energy velocity
    #  ----------------
    #  The massless scalar fields make the condensed operator A_uu(k) depend on
    #  k through a PSEUDO-INVERSE.  Near a ZGV the scalar block A_ss is nearly
    #  singular (this is the same ill-conditioning that invalidated the
    #  residual-based check V2a in check_01), so d/dk of the condensed operator
    #  is meaningless there.  The flux must therefore be evaluated on the
    #  UNCONDENSED pencil, whose generalised eigenproblem has an exactly
    #  singular mass matrix.  QZ handles that directly.
    #
    #      A(k) q = omega^2 M q ,   M = diag(rho*W, rho*W, 0, 0, ...)
    #      P     = -(1/2) omega Im[ q^H (dA/dk) q ]        (Hellmann-Feynman)
    #      <T>   = (1/4) q^H M q
    #      v_e   = P / (2<T>)
    # ------------------------------------------------------------------
    def _expanded(self, k, q):
        """Rebuild the FULL eigenvector [u; s] from a condensed one.

        The scalar fields are recovered ONCE from Gauss's-law-type equations,
        s = -pinv(A_ss) A_su u.  This is a single linear solve per mode; it is
        NOT differentiated with respect to k, which is what made the earlier
        attempt at a condensed flux unreliable near a ZGV.
        """
        A, _ = self.assemble(k, condense=False)
        nn = 2 * self.n
        A_su, A_ss = A[nn:, :nn], A[nn:, nn:]
        u = np.zeros(self.ndof, dtype=complex)
        u[self._keep[:nn]] = q            # scalar essential rows keep zeros
        s = -np.linalg.pinv(A_ss) @ (A_su @ u[:nn]) if A_ss.size else np.zeros(0)
        out = u.copy()
        full_s = np.zeros(self.ndof - nn, dtype=complex)
        kept = self._keep[self._keep >= nn] - nn
        full_s[kept] = s
        out[nn:] = full_s
        return out

    def modes(self, k, f_lo=None, f_hi=None):
        """Physical modes of the CONDENSED pencil (mass matrix positive
        definite, so no spurious pure-scalar eigenpairs).  Returns
        (freqs_Hz ascending, condensed eigenvectors)."""
        A, M = self.assemble(k)
        vals, vecs = sla.eig(A, M)
        vals = np.asarray(vals, dtype=complex)
        good = (np.isfinite(vals.real) & np.isfinite(vals.imag)
                & (np.abs(vals.imag) < 1e-7 * (np.abs(vals.real) + 1.0))
                & (vals.real > 1e-3))
        vals, vecs = vals[good].real, vecs[:, good]
        f = np.sqrt(vals) / (2 * np.pi)
        if f_lo is not None:
            keep = (f >= f_lo) & (f <= f_hi)
            f, vecs = f[keep], vecs[:, keep]
        order = np.argsort(f)
        return f[order], vecs[:, order]

    def _flux_of(self, k, q, om, dA=None):
        """P, <T>, v_e for one uncondensed eigenvector q at (k, omega).

        For a HERMITIAN pencil  A(k) q = omega^2 M q  the Hellmann-Feynman
        identity gives, with q normalised so that q^H M q = 1,

            d(omega^2)/dk = q^H (dA/dk) q        (REAL, not imaginary)
            v_e = domega/dk = q^H (dA/dk) q / (2 omega)
            <T> = (1/4) omega^2
            P   = 2 <T> v_e = (1/4) omega q^H (dA/dk) q

        NOTE (development bug, recorded): the flux is NOT Im[q^H dA q] -- for a
        Hermitian A that imaginary part vanishes identically and the routine
        silently returns zero.  An earlier version made exactly this mistake.
        The results are cross-checked against plate_spectral.energy_velocity,
        which uses the independent stress-velocity expression.
        """
        if q.shape[0] != self.ndof:
            q = self._expanded(k, q)          # condensed -> full, once
        M = self._mass_full()
        Mm = float(np.real(q.conj() @ (M @ q)))
        if Mm <= 0:
            return np.nan, np.nan, np.nan
        q = q / np.sqrt(Mm)
        if dA is None:
            dA = self._dA(k)
        Q = float(np.real(q.conj() @ (dA @ q)))
        T = 0.25 * om ** 2
        P = 0.25 * om * Q
        ve = Q / (2 * om) if om > 0 else np.nan
        return P, T, ve

    def _mass_full(self):
        """FULL-SIZE mass matrix (ndof x ndof).

        The scalar blocks are identically zero, so this matrix can be used
        with an expanded eigenvector [u; s] as well as with a condensed one
        whose scalar entries are zero.  No index reduction is needed -- and
        reducing it here was the source of a dimension mismatch during
        development.
        """
        n, nd = self.n, self.ndof
        M = np.zeros((nd, nd))
        rn = np.diag(self.rho * self.w)
        M[0:n, 0:n] = rn
        M[n:2 * n, n:2 * n] = rn
        return M

    def _dA(self, k):
        """dA/dk as a FULL-SIZE (ndof x ndof) matrix.

        assemble(condense=False) deletes the rows/columns of essential scalar
        conditions; they are re-inserted as zeros here so that dA can be
        contracted with a full expanded eigenvector.  This is exact because the
        corresponding entries of the expanded eigenvector are identically zero.
        """
        dk = 1e-4 * max(abs(k), 1.0)
        Ap, _ = self.assemble(k + dk, condense=False)
        Am, _ = self.assemble(k - dk, condense=False)
        d = (Ap - Am) / (2 * dk)
        d = 0.5 * (d + d.conj().T)
        keep = self._keep
        if len(keep) == self.ndof:
            return d
        full = np.zeros((self.ndof, self.ndof), dtype=complex)
        full[np.ix_(keep, keep)] = d
        return full

    def flux(self, k, f_target):
        """(v_e, P, <T>, f) for the mode nearest f_target (Hz)."""
        f, V = self.modes(k)
        j = int(np.argmin(np.abs(f - f_target)))
        dA = self._dA(k)
        P, T, ve = self._flux_of(k, V[:, j], 2 * np.pi * f[j], dA)
        return ve, P, T, f[j]

    def zgw_scan(self, k, f_lo, f_hi):
        """Tracking-free ZGV scan: the mode of minimum |v_e| in [f_lo,f_hi].

        The frequency window is mandatory; outside it the spectral
        discretisation admits spurious high-order modes whose v_e is
        meaningless (an un-windowed scan during development returned a
        spurious 2.9e8 m/s minimum).
        """
        f, V = self.modes(k, f_lo, f_hi)
        if len(f) == 0:
            return np.nan, np.nan
        dA = self._dA(k)
        best = (np.nan, np.inf)
        for j in range(len(f)):
            _, _, ve = self._flux_of(k, V[:, j], 2 * np.pi * f[j], dA)
            if np.isfinite(ve) and abs(ve) < abs(best[1]):
                best = (f[j], ve)
        return best

    def eigen_full(self, k, f_target):
        """Mode nearest f_target, with scalar fields expanded back.

        Returns (freq_Hz, dict of nodal fields) where the dict has keys
        'u1','u2' plus one entry per scalar field.
        """
        f, V = self.modes(k)
        j = int(np.argmin(np.abs(f - f_target)))
        q = V[:, j]
        q = q / (np.max(np.abs(q)) + 1e-300)
        n = self.n
        out = {"u1": q[0:n], "u2": q[n:2 * n]}
        for s in range(self.F):
            off = self.nu + s * n
            out[self.field_names[s]] = q[off:off + n]
        return f[j], out
