"""
Export the corrected modal inventory for Phase 2A.

This does NOT touch `validation/phase_01d_independent_comparison.csv`. That file
is frozen historical evidence of what the Phase 1D classifier produced, including
its defects, and MAC00 in the test suite pins its checksum so it cannot drift.
The corrected results are written here under a NEW EXPERIMENT ID.

Every row carries the four independent classification fields, the diagnostics
that justify them, and the provenance needed to reproduce it.

Reproduce:
    cd code && python3 -m scattering.export_inventory
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, ".")

import foundation_checks as fc                                    # noqa: E402
from scattering.modal import classify_mode, screen_halfspaces      # noqa: E402
from scattering.nullspace import null_vector                       # noqa: E402

EXPERIMENT_ID = "PHASE-02A-MODAL-CORRECTED-1"
DPS = 60
F_HZ = 1e6
MATERIALS = ("PZT-5H", "Al6061", "CdSe")
ANGLES = (0.0, 15.0, 30.0, 45.0, 60.0)

ROOT = Path(__file__).resolve().parents[2]

COLUMNS = [
    "experiment_id", "material", "theta_deg", "branch",
    "Re_xi", "Im_xi", "abs_xi", "attenuation_ratio", "attenuation_basis",
    "decay_length_m",
    "decay_length_in_transverse_wavelengths",
    "spatial_class", "radiation_dir", "energy_flux_status", "confidence",
    "mech_content", "elec_content", "thermal_content", "thermal_dispersion_ratio",
    "J_z", "normalisation_basis", "nullspace_agreement", "cramer_gap",
    "structural_columns", "convention_verified", "admissible_side",
    "phase_1d_label", "label_changed", "reasons",
]


def _commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, timeout=20).stdout.strip() or "UNKNOWN"
    except Exception:                                              # noqa: BLE001
        return "UNKNOWN"


def _phase1d_labels():
    """Read the historical labels so the audit can show what changed."""
    p = ROOT / "validation" / "phase_01d_independent_comparison.csv"
    out = {}
    if not p.exists():
        return out
    for r in csv.DictReader(p.open()):
        out[(r["material"], float(r["theta_deg"]), int(r["branch"]))] = \
            r["mode_class"]
    return out


def _f(x, sig=10):
    if x is None:
        return ""
    try:
        v = float(x)
    except Exception:                                              # noqa: BLE001
        return str(x)
    if v == 0.0:
        return "0.0"
    return f"{v:.{sig}e}"


def build_rows():
    hist = _phase1d_labels()
    commit = _commit()
    rows, screens = [], {}
    mp.mp.dps = DPS
    w = 2 * mp.pi * F_HZ
    for mat in MATERIALS:
        m = fc.material_from_table(mat)
        xi_th = mp.sqrt(w * mp.mpf(repr(m.rho)) * mp.mpf(repr(m.cE))
                        / mp.mpf(repr(m.kappa33)))
        for theta in ANGLES:
            k = w * mp.sin(mp.radians(theta)) / 4000.0
            roots = [mp.mpc(z) for z in m.xi_roots(k, w)]
            modes = []
            for br, xi in enumerate(roots):
                r = null_vector(fc.M_mp(m, k, xi, w), DPS)
                if r["vector"] is None:
                    modes.append(None)
                    continue
                a = r["vector"]
                jz = fc.flux_z_mp(m, k, xi, w, a)
                md = classify_mode(m, k, xi, w, a, jz=jz)
                modes.append(md)
            hs = screen_halfspaces(modes)
            screens[f"{mat}@{theta}"] = hs
            for br, md in enumerate(modes):
                if md is None:
                    continue
                old = hist.get((mat, float(theta), br), "ABSENT")
                side = ("z>0" if br in hs["z>0_admissible"] else
                        "z<0" if br in hs["z<0_admissible"] else
                        "non_decaying" if br in hs["non_decaying"] else
                        "gauge_excluded")
                # the historical vocabulary had one field; map it to the new
                # spatial class for comparison, and record whether it moved
                changed = old not in ("ABSENT", md["spatial_class"]) and not (
                    old == "PROPAGATING" and md["spatial_class"] == "WEAKLY_DAMPED"
                    and md["attenuation_ratio"] is not None
                    and float(md["attenuation_ratio"]) < 1e-6)
                rows.append({
                    "experiment_id": EXPERIMENT_ID,
                    "material": mat, "theta_deg": f"{theta:g}", "branch": br,
                    "Re_xi": _f(mp.re(xi)),
                    "Im_xi": _f(mp.im(xi)),
                    "abs_xi": _f(abs(xi)),
                    "attenuation_ratio": _f(md["attenuation_ratio"]),
                    "attenuation_basis": md["attenuation_basis"],
                    "decay_length_m": _f(md["decay_length_m"]),
                    "decay_length_in_transverse_wavelengths":
                        _f(md["decay_length_in_transverse_wavelengths"]),
                    "spatial_class": md["spatial_class"],
                    "radiation_dir": md["radiation_dir"],
                    "energy_flux_status": md["energy_flux_status"],
                    "confidence": md["confidence"],
                    "mech_content": _f(md["mech_content"]),
                    "elec_content": _f(md["elec_content"]),
                    "thermal_content": _f(md["thermal_content"]),
                    "thermal_dispersion_ratio": _f(md["thermal_dispersion_ratio"]),
                    "J_z": _f(jz),
                    "normalisation_basis": r["normalisation_basis"],
                    "nullspace_agreement": r["agreement"],
                    "cramer_gap": _f(r["cramer_gap"]),
                    "structural_columns": "|".join(map(str, r["structural_columns"])),
                    "convention_verified": "",
                    "admissible_side": side,
                    "phase_1d_label": old,
                    "label_changed": changed,
                    "reasons": "; ".join(md["reasons"]),
                })
    return rows, screens, commit


def main():
    rows, screens, commit = build_rows()
    out_csv = ROOT / "validation" / "phase_02a_modal_inventory.csv"
    with out_csv.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLUMNS)
        wr.writeheader()
        wr.writerows(rows)

    # tally: what actually changed against the historical labels
    # NB: in memory this is a bool; only the CSV writer turns it into "True".
    # Comparing against the string here silently reported zero changes.
    changed = [r for r in rows if r["label_changed"] is True]
    forbidden = [r for r in rows if r["spatial_class"] == "PROPAGATING"]
    print(f"{EXPERIMENT_ID}  commit {commit[:12]}")
    print(f"  rows written          : {len(rows)}  -> {out_csv.relative_to(ROOT)}")
    print(f"  labels changed vs 1D  : {len(changed)}")
    print(f"  PROPAGATING labels    : {len(forbidden)} (vocabulary retired)")
    for cls in ("EVANESCENT", "WEAKLY_DAMPED", "THERMAL_DIFFUSIVE", "GAUGE"):
        n = sum(1 for r in rows if r["spatial_class"] == cls)
        print(f"    {cls:18s}: {n}")
    manifest = {
        "experiment_id": EXPERIMENT_ID,
        "commit": commit,
        "working_precision_dps": DPS,
        "frequency_hz": F_HZ,
        "materials": list(MATERIALS),
        "angles_deg": list(ANGLES),
        "rows": len(rows),
        "labels_changed_vs_phase_1d": len(changed),
        "output": str(out_csv.relative_to(ROOT)),
        "output_sha256": hashlib.sha256(out_csv.read_bytes()).hexdigest(),
        "phase_1d_evidence_untouched":
            "validation/phase_01d_independent_comparison.csv",
        "halfspace_screens": screens,
        "reproduce": "cd code && python3 -m scattering.export_inventory",
    }
    mp_path = ROOT / "validation" / "phase_02a_modal_inventory_manifest.json"
    mp_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"  manifest             : {mp_path.relative_to(ROOT)}")
    print(f"  sha256               : {manifest['output_sha256'][:16]}...")


if __name__ == "__main__":
    main()
