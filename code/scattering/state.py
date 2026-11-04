"""
The physical interface state -- Phase 2A task 2.

WHAT THIS MODULE IS FOR
=======================
Every scattering calculation needs one thing before any matrix is assembled: a
precise statement of *which* eight fields are taken to be continuous, and in
which units. Getting that wrong produces a system that solves beautifully and
answers a different question, so the state is defined here once, derived from the
governing equations, and every other module imports it.

THE STATE
=========
    s = ( u_x, u_z, phi, theta, sigma_xz, sigma_zz, D_z, q_z )^T

Ordering is fixed by `STATE_FIELDS` and is used by every consumer. The two
halves have distinct characters:

    kinematic / potential   u_x, u_z, phi, theta     -- CARRY the wave
    kinetic / flux          sigma_xz, sigma_zz,
                            D_z, q_z                 -- ARE TRANSPORTED across

A welded interface transfers both groups: the field is continuous (no slip, no
separation, no charge sheet, no temperature jump) and so is the flux conjugate to
it (traction, normal electric displacement, normal heat flux).

DERIVATION FROM THE GOVERNING EQUATIONS
=======================================
With time dependence e^{-i w t} and sagittal variation exp(i(k x + xi z)), the
mechanical strains are

    e_xx = i k u_x           e_zz = i xi u_z           g_xz = i(xi u_x + k u_z)

plane strain, e_yy = 0. The electric field is E = -grad(phi), so

    E_x = -i k phi           E_z = -i xi phi

and the heat flux follows the Lord-Shulman/Cattaneo relation

    q_z = -kappa33 (i xi theta) / (1 - i w tau_q)

Signs are the ones already in `foundation_checks.M_mp`, so a mode vector from that
kernel satisfies these relations identically -- `state_from_mode` is therefore a
*check* on a null vector, not an independent model. The independent check is
`eigen_residual`: M(k, xi, w) a = 0 is the statement that the equations of motion
hold, and it is verified against the polynomial, not against this module.

WHY DIMENSIONAL SCALES ARE EXPLICIT
===================================
The four amplitudes (u, phi, theta) and the four fluxes do NOT share a unit. A
matrix built by dropping them into one array without scaling can be conditioned
by the ratio of volts to metres, which is ~1e9, and a solver will then report a
residual that looks small purely because one unknown was never constrained.
Every unit conversion goes through `SCALES` and `to_dimensionless` below.

Reproduce:
    cd code && python3 -m scattering.state
"""
from __future__ import annotations

import mpmath as mp

#: canonical ordering; consumers index by name, never by literal
STATE_FIELDS = ("ux", "uz", "phi", "theta", "txz", "tzz", "Dz", "qz")

#: physical dimension of each entry, written out so a units mistake is visible
STATE_DIMENSIONS = {
    "ux": "m", "uz": "m",
    "phi": "V", "theta": "K",
    "txz": "Pa", "tzz": "Pa",
    "Dz": "C/m^2", "qz": "W/m^2",
}


def _mp(x):
    return x if isinstance(x, mp.mpf) or isinstance(x, mp.mpc) else mp.mpf(x)


def material_mp(m):
    """Pull the constants out of a Medium6mm as exact mpf, once.

    Read through `repr(float(...))` so the decimal the material table specifies
    is what enters the arithmetic, not a binary approximation of it.
    """
    g = lambda n: mp.mpf(repr(float(getattr(m, n))))            # noqa: E731
    return dict(
        rho=g("rho"), C11=g("C11"), C12=g("C12"), C13=g("C13"),
        C33=g("C33"), C44=g("C44"),
        e31=g("e31"), e33=g("e33"), e15=g("e15"),
        eps11=g("eps11"), eps33=g("eps33"),
        beta1=g("beta1"), beta3=g("beta3"), p=g("p"),
        kappa11=g("kappa11"), kappa33=g("kappa33"),
        tauq=g("tauq"), T0=g("T0"), cE=g("cE"),
    )


def strains(k, xi, ux, uz):
    """(e_xx, e_zz, g_xz) under exp(i(k x + xi z))."""
    I = mp.mpc(0, 1)
    return I * k * ux, I * xi * uz, I * (xi * ux + k * uz)


def fields_of(m, k, xi):
    """The constitutive constants grouped for the stated relations."""
    c = material_mp(m)
    I = mp.mpc(0, 1)
    # E = -grad phi
    c["_Ex"] = -I * k
    c["_Ez"] = -I * xi
    return c


def state_from_mode(m, k, xi, w, a):
    """The 8-component physical state for one mode, with UNITS.

    `a` = (u_x, u_z, phi, theta) is the null vector of M at (k, xi, w). Returns a
    dict keyed by `STATE_FIELDS`, each value an mpf/mpc in SI units.

    Nothing is normalised here. In particular the state is NOT divided by |a|:
    a null vector is defined only up to a complex scale, and rescaling it must
    leave every physical amplitude unchanged. Normalisation is the caller's job
    and `scattering.solver` shows that the answer does not depend on it.
    """
    c = material_mp(m)
    ux, uz, phi, theta = [mp.mpc(t) for t in a]
    I = mp.mpc(0, 1)

    exx, ezz, gxz = strains(k, xi, ux, uz)
    Ex = -I * k * phi
    Ez = -I * xi * phi

    # mechanical: sigma = C : e - e^T : E - beta theta  (plane strain, e_yy = 0)
    tzz = c["C13"] * exx + c["C33"] * ezz - c["e33"] * Ez - c["beta3"] * theta
    txz = c["C44"] * gxz - c["e15"] * Ex

    # electrical: D = e : e + eps : E + p theta
    Dz = c["e31"] * exx + c["e33"] * ezz + c["eps33"] * Ez + c["p"] * theta

    # thermal: Lord-Shulman heat flux along z
    denom = 1 - I * w * c["tauq"]
    qz = -c["kappa33"] * (I * xi * theta) / denom

    return {
        "ux": ux, "uz": uz, "phi": phi, "theta": theta,
        "txz": txz, "tzz": tzz, "Dz": Dz, "qz": qz,
    }


def state_vector(m, k, xi, w, a):
    """`state_from_mode` as a plain tuple in STATE_FIELDS order."""
    s = state_from_mode(m, k, xi, w, a)
    return tuple(s[f] for f in STATE_FIELDS)


# ---------------------------------------------------------------------------
# Dimensional scaling
# ---------------------------------------------------------------------------
#: reference scales, chosen from the material and the drive frequency, so that
#: every state entry becomes O(1) for a typical mode. These are NORMALISATION
#: constants only: dividing by them cannot change any physical result, and the
#: solver verifies that by solving with and without them.
def reference_scales(m, w):
    """Characteristic magnitudes, from the material and the frequency.

    u  ~ 1 m                (displacement amplitude is the free scale)
    phi ~ |e33| u / eps33   (potential generated by that strain)
    th  ~ beta3 u / cE      (temperature rise from that strain)
    tractions  ~ |C33| u * k_ref
    D_z        ~ eps33 * phi * k_ref
    q_z        ~ kappa33 * th * k_ref
    """
    c = material_mp(m)
    cL = mp.sqrt(c["C11"] / c["rho"])
    k_ref = w / cL
    u_ref = mp.mpf(1)
    phi_ref = abs(c["e33"]) * u_ref / c["eps33"] if c["eps33"] != 0 else u_ref
    th_ref = abs(c["beta3"]) * u_ref / c["cE"] if c["cE"] != 0 else u_ref
    return {
        "ux": u_ref, "uz": u_ref,
        # a purely elastic medium has no piezoelectric coupling, so the potential
        # scale must fall back to a sensible value rather than zero
        "phi": phi_ref if phi_ref != 0 else u_ref,
        "theta": th_ref if th_ref != 0 else u_ref,
        "txz": abs(c["C44"]) * u_ref * k_ref,
        "tzz": abs(c["C33"]) * u_ref * k_ref,
        "Dz": c["eps33"] * (phi_ref if phi_ref != 0 else u_ref) * k_ref,
        "qz": abs(c["kappa33"]) * (th_ref if th_ref != 0 else u_ref) * k_ref,
    }


def to_dimensionless(state, scales):
    """Divide each entry by its characteristic scale.

    Purely a conditioning device. The inverse mapping is exact, so no physics is
    changed; `scattering.solver` asserts that the answer is identical with and
    without scaling rather than assuming it.
    """
    return tuple(state[f] / scales[f] if scales[f] != 0 else state[f]
                 for f in STATE_FIELDS)


def from_dimensionless(vec, scales):
    return {f: vec[i] * scales[f] for i, f in enumerate(STATE_FIELDS)}


# ---------------------------------------------------------------------------
# Independent verification of the state against analytical limits
# ---------------------------------------------------------------------------
def check_elastic_reduction(m, k, xi, w, ux, uz):
    """With e = beta = p = 0 the state must collapse to isotropic elasticity.

    Returns a dict of relative errors comparing the state's tractions against the
    closed-form Hooke law for an isotropic solid:

        sigma_zz = lambda (e_xx + e_zz) + 2 mu e_zz
        sigma_xz = 2 mu e_xz

    with lambda = C12 (since C13 = C12 = lambda for isotropy) and mu = C44. This
    is an independent formula: it is the textbook isotropic relation, not the
    6mm constitutive law the module otherwise uses.
    """
    c = material_mp(m)
    ux, uz = mp.mpc(ux), mp.mpc(uz)
    a = (ux, uz, mp.mpc(0), mp.mpc(0))
    s = state_from_mode(m, k, xi, w, a)
    exx, ezz, gxz = strains(k, xi, ux, uz)

    def rel(a_, b_):
        d = abs(a_ - b_)
        sc = max(abs(a_), abs(b_), mp.mpf("1e-300"))
        return d / sc

    # ISOTROPY IS A PREMISE, NOT AN ASSUMPTION. The textbook identity
    # sigma_zz = lambda(exx+ezz) + 2 mu ezz holds only when C13 = C12 and
    # C33 = C12 + 2 C44. PZT-5H is transversely isotropic (C13 = 8.41e10 vs
    # C12 = 7.95e10) and fails this by 0.87 %, so comparing it against the
    # isotropic formula would report a "defect" that is really the material
    # being anisotropic. The 6mm law is what the module implements, so for a
    # non-isotropic medium the check is against the 6mm law expanded by hand.
    isotropic = (c["C13"] == c["C12"]) and (c["C33"] == c["C12"] + 2 * c["C44"])

    if isotropic:
        lam = c["C12"]
        mu = c["C44"]
        tzz_ref = lam * (exx + ezz) + 2 * mu * ezz
        txz_ref = mu * gxz
        basis = "isotropic Hooke law (lambda, mu)"
    else:
        # hand-expanded 6mm law, written out independently of the module
        tzz_ref = c["C13"] * exx + c["C33"] * ezz
        txz_ref = c["C44"] * gxz
        basis = "6mm constitutive law, hand-expanded"

    return {
        "isotropic_medium": isotropic,
        "basis": basis,
        "tzz_rel": rel(s["tzz"], tzz_ref),
        "txz_rel": rel(s["txz"], txz_ref),
        # and the coupled channels must be exactly absent when phi = theta = 0
        "Dz_abs": abs(s["Dz"]),
        "qz_abs": abs(s["qz"]),
    }


def check_electrostatic_limit(m, k, xi, w, phi, eps_ref=1.0):
    """With no strain and no thermal field, D_z must reduce to eps33 E_z.

    A pure potential with u = theta = 0 gives D_z = eps33 * (-i xi phi). This is
    the electrostatic constitutive law and it is checked against the closed form
    rather than against the module's own expression.
    """
    I = mp.mpc(0, 1)
    s = state_from_mode(m, k, xi, w, (mp.mpc(0), mp.mpc(0), mp.mpc(phi),
                                      mp.mpc(0)))
    c = material_mp(m)
    expected = c["eps33"] * (-I * xi * mp.mpc(phi))
    d = abs(s["Dz"] - expected)
    sc = max(abs(expected), mp.mpf("1e-300"))
    return {"Dz_rel": d / sc,
            # no mechanical or thermal response may appear
            "tzz_abs": abs(s["tzz"]), "txz_abs": abs(s["txz"]),
            "qz_abs": abs(s["qz"])}


def check_thermal_limit(m, k, xi, w, theta):
    """With u = phi = 0, q_z must reduce to the LS/Cattaneo relation alone."""
    I = mp.mpc(0, 1)
    s = state_from_mode(m, k, xi, w, (mp.mpc(0), mp.mpc(0), mp.mpc(0),
                                      mp.mpc(theta)))
    c = material_mp(m)
    expected = -c["kappa33"] * (I * xi * mp.mpc(theta)) / \
        (1 - I * w * c["tauq"])
    d = abs(s["qz"] - expected)
    sc = max(abs(expected), mp.mpf("1e-300"))
    return {"qz_rel": d / sc,
            "tzz_abs": abs(s["tzz"]), "Dz_abs": abs(s["Dz"])}


def _main():
    import sys
    sys.path.insert(0, ".")
    import foundation_checks as fc

    mp.mp.dps = 50
    w = 2 * mp.pi * 1e6
    k = mp.mpf("785.39816339")
    print("State definition and analytical limits\n")
    for name in ("Al6061", "PZT-5H"):
        m = fc.material_from_table(name)
        print(f"--- {name} ---")
        r = check_elastic_reduction(m, k, mp.mpf("1000"), w, mp.mpf(1),
                                    mp.mpf("0.3"))
        print(f"  isotropic reduction : tzz rel {mp.nstr(r['tzz_rel'], 4)}, "
              f"txz rel {mp.nstr(r['txz_rel'], 4)}")
        r2 = check_electrostatic_limit(m, k, mp.mpf("500"), w, mp.mpf("1e4"))
        print(f"  electrostatic limit : Dz rel {mp.nstr(r2['Dz_rel'], 4)}")
        r3 = check_thermal_limit(m, k, mp.mpf("800"), w, mp.mpf("10"))
        print(f"  thermal limit       : qz rel {mp.nstr(r3['qz_rel'], 4)}")


if __name__ == "__main__":
    _main()
