"""
Measure the divergence between the two modal-classification conventions.

WHY THIS FILE EXISTS
====================
Task 2 changed two classification rules in `scattering.modal` and
`scattering.radiation`. Both changes were forced by defects the coupled
assembly exposed, and both alter the labels of rows in the FROZEN Task 1R
inventory (`validation/phase_02a_modal_inventory_v2.csv`).

The brief forbids silently replacing historical evidence. So the artifact is
left exactly as published, `test_13_inventory_integrity.py` re-derives its
labels under the convention that generated it (`legacy_gauge=True`,
`strict_legacy=True`), and the divergence is measured here and reported rather
than assumed to be harmless.

THE TWO CONVENTIONS
===================
RULE 1 -- GAUGE vs ELECTROSTATIC
    legacy: any mechanically decoupled mode is GAUGE.
    current: mechanically decoupled AND D_z = -eps33*xi*phi = 0 is GAUGE;
             otherwise it is a decoupled Laplace ELECTROSTATIC mode.
    why: at k != 0 the pair at xi = -/+ ik carries D_z = -/+ i eps33 k phi,
    which is NOT zero, so it satisfies the D_z continuity equation like any
    other mode. Calling it a gauge let the assembly delete the electrostatic
    channel of an elastic half-space while keeping the equivalent channel of a
    piezoelectric one -- an exclusion of a physical channel. Measured
    example: Al6061 at 30 deg, D_z = -6.954e-9 = -eps33*k exactly.

RULE 2 -- LOCALISED vs PROPAGATING at small but non-zero Im(xi)
    legacy: any Im(xi) != 0 is localised; only exactly lossless modes transport.
    current: attenuation below ATTENUATION_PROPAGATING (1e-4 per wavelength)
             is propagation.
    why: thermomechanical coupling puts every acoustic mode just off
    Im(xi) = 0. PZT-5H's shear wave at 30 deg has attenuation 3e-13 per
    wavelength and a decay length of 9.6e8 m. The strict rule declared it
    non-transporting on the far side, which removed the incident channel from
    the coupled problem entirely. The threshold sits ten orders above that case
    and four below the weakest genuinely damped mode, so nothing is near it.

Both legacy behaviours remain reachable through explicit keyword arguments, so
the published artifact stays reproducible and the change is auditable rather
than invisible.

Run from the repository root:
    python3 validation/scattering/convention_divergence.py
"""
from __future__ import annotations

import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "code"))

import mpmath as mp                                     # noqa: E402

mp.mp.dps = 40

import foundation_checks as fc                          # noqa: E402
from scattering.modal import classify_mode               # noqa: E402
from scattering.nullspace import null_vector             # noqa: E402
from scattering.radiation import (                       # noqa: E402
    ATTENUATION_PROPAGATING, classify_role,
)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
INV = os.path.join(ROOT, "validation",
                   "phase_02a_modal_inventory_v2.csv")
OUT_JSON = os.path.join(ROOT, "validation",
                        "phase_02a_task2_convention_divergence.json")

W = 2 * mp.pi * 1e6


def measure():
    rows = list(csv.DictReader(open(INV)))
    class_diff = []
    role_diff = []
    band = []
    for r in rows:
        mat, th = r["material"], mp.mpf(r["theta_deg"])
        br = int(r["branch"])
        m = fc.material_from_table(mat)
        k = W * mp.sin(mp.radians(th)) / mp.mpf(4000)
        xi = m.xi_roots(k, W)[br]
        a = null_vector(fc.M_mp(m, k, xi, W), mp.mp.dps)["vector"]
        if a is None:
            continue
        jz = fc.flux_z_mp(m, k, xi, W, a)
        md_new = classify_mode(m, k, xi, W, a, jz=jz)
        md_old = classify_mode(m, k, xi, W, a, jz=jz, legacy_gauge=True)
        tag = f"{mat}@{r['theta_deg']} br{br}"

        rec = {"row": tag, "old": md_old["spatial_class"],
               "new": md_new["spatial_class"],
               "legacy_csv": r["spatial_class"]}
        if md_old["spatial_class"] != md_new["spatial_class"]:
            class_diff.append(rec)

        # attenuation scale, to show WHY rule 2 fires
        re_, im_ = abs(mp.re(xi)), abs(mp.im(xi))
        att = (im_ / re_) if re_ else mp.inf
        if 0 < att < 1e-4:
            band.append({"row": tag, "attenuation_per_wavelength":
                         mp.nstr(att, 4),
                         "decay_length_m": mp.nstr(1 / im_, 4) if im_ else "inf",
                         "wavelength_m": mp.nstr(2 * mp.pi / re_, 4)})

        # roles, under each convention, on both sides
        for side, pre in (("z>0", "z>0:"), ("z<0", "z<0:")):
            ro, _ = classify_role(md_old, side, jz, strict_legacy=True)
            rn, _ = classify_role(md_new, side, jz)
            if ro != rn:
                role_diff.append({"row": tag, "side": side,
                                  "old": ro, "new": rn})
    return class_diff, role_diff, band


def main():
    class_diff, role_diff, band = measure()
    print("Convention divergence -- task 2 classification changes\n")
    print(f"  inventory rows                        : 120")
    print(f"  spatial_class changed                 : {len(class_diff)}")
    print(f"  side roles changed                    : {len(role_diff)}")
    print(f"  rows inside the propagation band      : {len(band)}")
    print(f"  ATTENUATION_PROPAGATING               : "
          f"{float(ATTENUATION_PROPAGATING):g} per wavelength\n")

    if class_diff:
        print("  spatial_class changes (old -> new):")
        pairs = {}
        for d in class_diff:
            pairs[(d["old"], d["new"])] = pairs.get((d["old"], d["new"]), 0) + 1
        for (o, n), c in sorted(pairs.items()):
            print(f"    {o} -> {n}: {c}")
    if role_diff:
        pairs = {}
        for d in role_diff:
            pairs[(d["old"], d["new"])] = pairs.get((d["old"], d["new"]), 0) + 1
        print("\n  side-role changes (old -> new):")
        for (o, n), c in sorted(pairs.items()):
            print(f"    {o} -> {n}: {c}")
    if band:
        atts = sorted(mp.mpf(d["attenuation_per_wavelength"]) for d in band)
        print(f"\n  attenuation of the {len(band)} affected rows:")
        print(f"    min {mp.nstr(atts[0], 3)}   max {mp.nstr(atts[-1], 3)}"
              f"  per wavelength")
        print(f"    decay lengths {band[0]['decay_length_m']} .. "
              f"{band[-1]['decay_length_m']} m")

    doc = {
        "artifact": "phase_02a_task2_convention_divergence",
        "inventory": os.path.relpath(INV, ROOT),
        "inventory_unchanged": True,
        "note": "the published Task 1R inventory was NOT regenerated; the two "
                "legacy rules remain reachable via legacy_gauge / "
                "strict_legacy so the artifact stays reproducible",
        "n_rows": 120,
        "n_spatial_class_changed": len(class_diff),
        "n_role_changed": len(role_diff),
        "n_rows_in_propagation_band": len(band),
        "attenuation_propagating_threshold": float(ATTENUATION_PROPAGATING),
        "spatial_class_changes": class_diff,
        "role_changes": role_diff,
        "propagation_band": band,
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(doc, fh, indent=1)
    print(f"\nwrote {os.path.relpath(OUT_JSON, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
