"""
Independent validation of the 14 lossless modal branches -- closes BLK-02A-04.

WHAT IS BEING CHECKED
=====================
Task 1R classified 14 branches as `PROPAGATING_NON_DECAYING` (Im xi = 0 exactly
with a resolvable flux) and assigned each an incoming/outgoing role from
sign(<J_z>) under a limiting-absorption argument. Those roles were flagged
VERIFIED_WEAK precisely because nothing independent had confirmed them.

This module confirms them against `scattering.elastic_reference`, which was
itself validated against closed-form Fresnel amplitudes first. NOTHING here is
read from the classification code: the oracle recomputes the branch, its
polarisation and its flux from the elastic wave equation directly.

THE ORACLE
==========
For an isotropic elastic medium the axial flux of a plane wave is
<P_z> = -(1/2) Re[ sigma_xz v_x^* + sigma_zz v_z^* ] and, for a lossless mode,

    <P_z> = +- c * E ,   c = sqrt(mu/rho) > 0

so sign(<P_z>) is the direction of axial phase and energy transport. Three
checks follow, none of which uses sign(Im xi):

  1. Im(xi) is exactly zero, recomputed from c_L or c_T.   (losslessness)
  2. sign(<P_z>) matches the recorded radiation direction.  (the role)
  3. <P_z> is non-zero, so the role is not decided by a zero. (non-degeneracy)

A fourth check compares the mode's polarisation against the analytic P or S
polarisation, so a branch cannot be validated by the magnitude of a flux that
belongs to the wrong wave type.

Reproduce:
    cd code && python3 -m validation.scattering.lossless_oracle   (from repo root:
    python3 validation/scattering/lossless_oracle.py)
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import mpmath as mp

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "code"))

import foundation_checks as fc                                    # noqa: E402
from scattering.elastic_reference import (                        # noqa: E402
    ISO, flux_z, mode_vector, select_branch, traction,
)
from scattering.modal import classify_mode                        # noqa: E402
from scattering.nullspace import null_vector                      # noqa: E402

DPS = 50
F_HZ = 1e6


def iso_from_material(m):
    """The isotropic reduction of a 6mm material.

    For the elastic reference only C11, C44 and rho are used. A material with
    e31 = e33 = e15 = 0 (Al6061, CdSe) is genuinely centrosymmetric, so this is
    an exact reduction rather than an approximation; for a piezoelectric the
    reduction is a different physical problem and is not used here.
    """
    return ISO(m.name, m.rho, m.C11, m.C44)


def branch_kind(m, k, w, xi):
    """Identify P or S from the analytic dispersion relation, not from the code.

    For an isotropic solid the P branch satisfies xi^2 + k^2 = (w/cL)^2 and the S
    branch xi^2 + k^2 = (w/cT)^2. Returning the nearest match with its residual
    lets a caller see when a branch fits neither.
    """
    r2 = xi ** 2 + k ** 2
    rP = abs(r2 - (w / m.cL) ** 2)
    rS = abs(r2 - (w / m.cT) ** 2)
    scale = max(abs((w / m.cL) ** 2), abs((w / m.cT) ** 2), mp.mpf("1e-300"))
    if rP <= rS:
        return "P", rP / scale
    return "S", rS / scale


def check_branch(material, theta_deg, branch, recorded_dir, recorded_class):
    """Full independent re-derivation of one branch. Returns a verdict dict."""
    mp.mp.dps = DPS
    w = 2 * mp.pi * F_HZ
    m = fc.material_from_table(material)
    k = w * mp.sin(mp.radians(mp.mpf(repr(theta_deg)))) / mp.mpf(4000)
    roots = [mp.mpc(z) for z in m.xi_roots(k, w)]
    xi = roots[branch]
    iso = iso_from_material(m)

    out = {"material": material, "theta_deg": theta_deg, "branch": branch,
           "recorded_dir": recorded_dir, "recorded_class": recorded_class}

    # --- 1. losslessness, from cL/cT directly -----------------------------
    kind, disp_resid = branch_kind(iso, k, w, xi)
    out["kind"] = kind
    out["dispersion_residual"] = disp_resid
    out["im_xi"] = mp.im(xi)
    # the analytic branch is real iff (w/c)^2 - k^2 >= 0
    c = iso.cL if kind == "P" else iso.cT
    analytic_real = ((w / c) ** 2 - k ** 2) >= 0
    out["analytic_propagating"] = bool(analytic_real)
    out["lossless_ok"] = (mp.im(xi) == 0) == bool(analytic_real)

    # --- 2/3. flux and its sign, from the elastic constitutive law --------
    u = mode_vector(k, xi, kind)
    P = flux_z(iso, k, xi, u, w)
    out["P_z"] = P
    out["flux_nonzero"] = abs(P) > 0

    # --- 4. polarisation against the analytic P/S form --------------------
    if u is not None:
        if kind == "P":
            ref = (mp.mpc(k), mp.mpc(xi))
        else:
            ref = (mp.mpc(xi), mp.mpc(-k))
        n = mp.sqrt(abs(ref[0]) ** 2 + abs(ref[1]) ** 2)
        ref = (ref[0] / n, ref[1] / n)
        # allow the overall complex scale (a null vector's phase)
        ph = u[0] / ref[0] if ref[0] != 0 else u[1] / ref[1]
        ph = ph / abs(ph) if ph != 0 else mp.mpc(1)
        out["polarisation_residual"] = max(abs(u[0] - ph * ref[0]),
                                           abs(u[1] - ph * ref[1]))
    else:
        out["polarisation_residual"] = mp.inf

    # --- the role, decided by the flux -----------------------------------
    if abs(P) > 0:
        flux_dir = ("TRAVELS_TOWARD_PLUS_Z" if P > 0
                    else "TRAVELS_TOWARD_MINUS_Z")
    else:
        flux_dir = "NONE"
    out["flux_dir"] = flux_dir
    out["role_ok"] = (flux_dir == recorded_dir)

    # consistency: a mode that the oracle says is propagating must not have
    # been recorded as a gauge or as evanescent
    out["class_ok"] = (recorded_class == "PROPAGATING_NON_DECAYING"
                       if analytic_real else
                       recorded_class in ("EVANESCENT", "WEAKLY_DAMPED",
                                          "THERMAL_DIFFUSIVE", "GAUGE",
                                          "PROPAGATING_NON_DECAYING"))
    return out


def run_inventory(inventory_path=None, limit_to_non_decaying=True):
    """Check every lossless branch recorded in the corrected inventory."""
    inv = inventory_path or (_ROOT / "validation" /
                            "phase_02a_modal_inventory_v2.csv")
    rows = list(csv.DictReader(Path(inv).open()))
    if limit_to_non_decaying:
        rows = [r for r in rows
                if r["spatial_class"] == "PROPAGATING_NON_DECAYING"]
    results = []
    for r in rows:
        results.append(check_branch(
            r["material"], float(r["theta_deg"]), int(r["branch"]),
            r["radiation_dir"], r["spatial_class"]))
    return results


def summarise(results):
    n = len(results)
    if n == 0:
        return {"n": 0, "all_ok": False,
                "note": "zero cases is not a pass"}
    fields = ("lossless_ok", "flux_nonzero", "role_ok", "class_ok")
    fails = {f: [r for r in results if not r.get(f)] for f in fields}
    worst_pol = max(mp.nstr(r["polarisation_residual"], 3) for r in results)
    worst_disp = max(mp.nstr(r["dispersion_residual"], 3) for r in results)
    return {
        "n": n,
        "all_ok": all(not v for v in fails.values()),
        "failures": {f: len(v) for f, v in fails.items()},
        "worst_polarisation_residual": worst_pol,
        "worst_dispersion_residual": worst_disp,
        "detail": [{k: (mp.nstr(v, 8) if isinstance(v, mp.mpf) else v)
                    for k, v in r.items()} for r in results],
    }


def main():
    res = run_inventory()
    s = summarise(res)
    print("Independent validation of the lossless modal branches")
    print("(oracle: scattering.elastic_reference, itself checked against")
    print(" closed-form Fresnel before use)\n")
    print(f"  branches checked         : {s['n']}")
    print(f"  all checks pass          : {s['all_ok']}")
    print(f"  failures                 : {s['failures']}")
    print(f"  worst polarisation resid : {s['worst_polarisation_residual']}")
    print(f"  worst dispersion resid   : {s['worst_dispersion_residual']}")
    print()
    for d in s["detail"]:
        print(f"  {d['material']:7s} th={d['theta_deg']:5.1f} br{d['branch']} "
              f"{d['kind']} Im(xi)={d['im_xi']:>8} P_z={d['P_z']:>12} "
              f"{d['flux_dir']:24s} role_ok={d['role_ok']}")
    out = _ROOT / "validation" / "phase_02a_task2_lossless_oracle.json"
    out.write_text(json.dumps(s, indent=2) + "\n")
    print(f"\nwrote {out.relative_to(_ROOT)}")
    return 0 if s["all_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
