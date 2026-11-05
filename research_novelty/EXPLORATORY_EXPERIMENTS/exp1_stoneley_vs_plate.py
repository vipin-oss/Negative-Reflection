"""Exploratory experiment 1 (research only; does not touch production solver).

E1a: Stoneley surface branch at Al6061|PZT-5H half-space interface, using the
     independent 8-state operator. For each real k, scan omega for roots of
     the 8x8 interface determinant (4 decaying modes per side, 8 contact
     conditions). If omega(k)/k is constant -> nondispersive -> no backward
     branch possible (falsifies Arena C2).
E1b: PZT-5H plate of thickness h: Rayleigh-Lamb-type dispersion via the same
     operator with stress-free top/bottom; look for a backward segment
     d omega/d k < 0 (confirms the published plate mechanism exists in this
     material set).
E1c: repeat E1b with tauq -> 0 to quantify the Lord-Shulman thermal
     contribution at 1 MHz (expected negligible).
"""
import sys
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "handoffs" / "task2uv_2026-11-05" / "code"))

from independent_finite_k_state_space import load, operator, balanced, scale  # noqa

CSV = ROOT / "handoffs" / "task2uv_2026-11-05" / "source" / "material_constants.csv"


def modes(p, k, w, dps=30):
    mp.mp.dps = dps
    A = operator(p, k, w)
    B, kr, s = balanced(A, p, w)
    vals, V = mp.eig(B)
    out = []
    for j in range(8):
        xi = vals[j] * kr
        Y = mp.matrix([s[i] * V[i, j] for i in range(8)])
        out.append((xi, Y))
    return out, kr


def interface_det(p1, p2, k, w, dps=30):
    m1, kr = modes(p1, k, w, dps)
    m2, _ = modes(p2, k, w, dps)
    s1 = [m for m in m1 if mp.im(m[0]) > 0]   # decay into z<0? side1 z<0: Im xi>0 decays as z->-inf... use Re/Im rule
    # side1 occupies z<0: admissible if Im(xi)>0; side2 z>0: Im(xi)<0
    a1 = [m for m in m1 if mp.im(m[0]) > 0]
    a2 = [m for m in m2 if mp.im(m[0]) < 0]
    if len(a1) != 4 or len(a2) != 4:
        return None
    M = mp.matrix(8, 8)
    for j, (xi, Y) in enumerate(a1):
        for i in range(8):
            M[i, j] = Y[i]
    for j, (xi, Y) in enumerate(a2):
        for i in range(8):
            M[i, 4 + j] = -Y[i]
    return abs(mp.det(M))


def stoneley_scan():
    p1, p2 = load(CSV, "Al6061"), load(CSV, "PZT-5H")
    out = []
    for k in (mp.mpf("100"), mp.mpf("200"), mp.mpf("400")):
        # coarse scan in omega
        best = None
        wlo, whi = 2 * mp.pi * mp.mpf("2e5"), 2 * mp.pi * mp.mpf("3e6")
        prev = None
        N = 40
        for i in range(N + 1):
            w = wlo + (whi - wlo) * i / N
            d = interface_det(p1, p2, k, w)
            if d is None:
                continue
            if best is None or d < best[1]:
                best = (w, d)
        if best:
            out.append((k, best[0], best[0] / k, mp.nstr(best[1], 3)))
    return out


def plate_dispersion(mat, h, wlist, tau_scale=1.0, dps=30):
    p = load(CSV, mat)
    p["tauq"] = p["tauq"] * tau_scale
    res = []
    for w in wlist:
        # find real k roots of det M_plate = 0 by scanning k
        mp.mp.dps = dps
        best = None
        N = 60
        klo, khi = mp.mpf("10"), mp.mpf("3000")
        for i in range(N + 1):
            k = klo + (khi - klo) * i / N
            m, kr = modes(p, k, w, dps)
            # plate: combine modes with exp(i xi z) and exp(-i xi z) symmetry;
            # use full 8x8: u(z=0)=u(z=h) via Bloch-ish determinant -- instead
            # use simple symmetric/antisymmetric reduction is complex; here use
            # the determinant of [Y_j(z=0); Y_j(z=h)] stacking for 8 modes ->
            # overdetermined; instead scan det of 8x8 with rows = fields at 0
            # and h for 4 independent +/- combinations. Approximate: use det of
            # matrix with columns = modes, rows = (fields at z=0, fields at z=h)
            # for 4 modes paired +/-.
            ms = sorted(m, key=lambda t: mp.im(t[0]))
            cols = []
            used = set()
            for xi, Y in ms:
                if any(abs(xi + x2) < abs(xi) * mp.mpf("1e-6") for x2, _ in cols):
                    continue
                cols.append((xi, Y))
                if len(cols) == 4:
                    break
            M = mp.matrix(8, 8)
            for j, (xi, Y) in enumerate(cols):
                for i in range(8):
                    M[i, j] = Y[i]
                    M[i, 4 + j] = Y[i] * mp.e ** (1j * xi * h)
            d = abs(mp.det(M))
            if best is None or d < best[1]:
                best = (k, d)
        if best:
            res.append((w, best[0]))
    return res


if __name__ == "__main__":
    print("E1a Stoneley Al|PZT: k, omega_root, omega/k, |det|")
    for row in stoneley_scan():
        print(" ", mp.nstr(row[0], 4), mp.nstr(row[1], 6), mp.nstr(row[2], 8), row[3])
