"""
Export the corrected modal inventory for Phase 2A (task 1R).

SUPERSEDES `PHASE-02A-MODAL-CORRECTED-1`, WHICH WAS INVALID
===========================================================
The first exporter had a scoping defect: it computed each branch's root and flux
into loop-local variables, stored only the CLASSIFICATION in a list, and then
re-read those same variables in a second loop after they had finished iterating.
Every row in a material/angle block therefore carried the LAST branch's
`Re_xi`, `Im_xi`, `abs_xi` and `J_z`. In the shipped CSV all eight PZT-5H rows at
30 degrees carried branch 7's `Re_xi = +2.569e6` -- including branch 0, whose
true `Re_xi` is -2.569e6, i.e. the wrong sign.

The repair is structural, not a patch: each branch is now ONE immutable record
that carries its index, root, eigenvector, flux, classification, nullspace
diagnostics, radiation direction and confidence together. Rows are generated
from those records, so there is no window in which a value can be read from a
different branch. `_Branch` is a frozen dataclass and `build_rows()` reads it by
field, never by position in a parallel list.

Nothing is recovered by parsing the prose `reasons` field -- that text is itself
generated FROM the record, and a test asserts the two cannot contradict.

Evidence preservation: the Phase 1D CSV and the invalid Task 1 CSV are both left
on disk untouched. This module writes a new immutable artifact under a new
experiment ID.

Reproduce:
    cd code && python3 -m scattering.export_inventory
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import mpmath as mp

sys.path.insert(0, ".")

import foundation_checks as fc                                    # noqa: E402
from scattering.modal import classify_mode, screen_halfspaces      # noqa: E402
from scattering.nullspace import null_vector                       # noqa: E402
from scattering.radiation import screen as radiation_screen        # noqa: E402

#: New experiment ID. The Task 1 artifact `PHASE-02A-MODAL-CORRECTED-1` is
#: marked INVALID/SUPERSEDED in the evidence registry and is not overwritten.
EXPERIMENT_ID = "PHASE-02A-MODAL-CORRECTED-2"

SUPERSEDES = "PHASE-02A-MODAL-CORRECTED-1"

DPS = 60
F_HZ = 1e6
MATERIALS = ("PZT-5H", "Al6061", "CdSe")
ANGLES = (0.0, 15.0, 30.0, 45.0, 60.0)

ROOT = Path(__file__).resolve().parents[2]

COLUMNS = [
    "experiment_id", "material", "theta_deg", "branch",
    "Re_xi", "Im_xi", "abs_xi", "attenuation_ratio", "attenuation_basis",
    "decay_length_m", "decay_length_in_transverse_wavelengths",
    "spatial_class", "radiation_dir", "radiation_source",
    "energy_flux_status", "confidence",
    "mech_content", "elec_content", "thermal_content", "thermal_dispersion_ratio",
    "J_z", "J_z_sign", "normalisation_basis", "nullspace_agreement", "cramer_gap",
    "gram_residual", "structural_columns", "admissible_side",
    "outgoing_status", "phase_1d_label", "label_changed", "reasons",
]


@dataclass(frozen=True)
class _Branch:
    """ONE branch of the dispersion relation, with everything measured about it.

    Frozen so a field cannot be rebound to another branch's value after
    construction, and so the exporter can only read it by name. Every quantity
    below is computed for this root and only this root.
    """
    material: str
    theta_deg: float
    branch: int
    xi: Any                       # the complex root itself
    eigenvector: Any              # the 4-component null vector
    nullspace: Any                # diagnostics dict from null_vector()
    jz: Any                       # cycle-averaged <J_z> for THIS mode
    classification: Any           # dict from classify_mode() for THIS mode
    admissibility: str            # computed from the whole branch set
    outgoing_status: str          # computed from the whole branch set
    extra: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
def _commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, timeout=20).stdout.strip() or "UNKNOWN"
    except Exception:                                              # noqa: BLE001
        return "UNKNOWN"


def _phase1d_labels():
    """Historical labels, for the audit column only. Never used as data."""
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


def build_branch_records(material, theta_deg):
    """Compute every branch for one (material, angle) and return records plus
    the half-space screen.

    The screen needs the whole set, so it is applied as a SECOND, explicit pass
    that rebuilds each record via `dataclasses.replace`. Nothing is read from a
    stale loop variable: each record already owns its root and flux.
    """
    from dataclasses import replace

    w = 2 * mp.pi * F_HZ
    k = w * mp.sin(mp.radians(theta_deg)) / 4000.0
    m = fc.material_from_table(material)
    roots = [mp.mpc(z) for z in m.xi_roots(k, w)]

    records = []
    for br, xi in enumerate(roots):
        r = null_vector(fc.M_mp(m, k, xi, w), DPS)
        a = r["vector"]
        if a is None:
            # keep the slot so later indices still refer to the same root
            records.append(None)
            continue
        jz = fc.flux_z_mp(m, k, xi, w, a)
        md = classify_mode(m, k, xi, w, a, jz=jz)
        records.append(_Branch(
            material=material, theta_deg=theta_deg, branch=br,
            xi=xi, eigenvector=a, nullspace=r, jz=jz, classification=md,
            admissibility="PENDING", outgoing_status="PENDING"))

    classifications = [rec.classification if rec is not None else None
                       for rec in records]
    jzs = [rec.jz if rec is not None else None for rec in records]
    screen = screen_halfspaces(classifications, jzs)
    roles = radiation_screen(classifications, jzs)

    for br, rec in enumerate(records):
        if rec is None:
            continue
        if br in screen["gauge_excluded"]:
            side = "gauge_excluded"
        elif br in screen["non_decaying"]:
            # bounded on both sides; the informative field is the role, not a
            # side. Saying "non_decaying" and stopping is the P0-3 defect.
            side = "both (non-decaying)"
        elif br in screen["z>0_admissible"]:
            side = "z>0"
        elif br in screen["z<0_admissible"]:
            side = "z<0"
        else:
            side = "excluded"
        # the outgoing status is per-side; record both so neither is implied
        up = roles["z>0"]
        dn = roles["z<0"]
        def role_of(side_roles, idx):
            if idx in side_roles["outgoing"]:
                return "OUTGOING"
            if idx in side_roles["incoming"]:
                return "INCOMING"
            return side_roles["excluded"].get(idx, ("UNRESOLVED", ""))[0]
        records[br] = replace(
            rec, admissibility=side,
            outgoing_status=(f"z>0:{role_of(up, br)}; z<0:{role_of(dn, br)}"))
    out = dict(screen)
    out["roles"] = roles
    return records, out


def build_rows():
    """Generate CSV rows directly from the branch records."""
    hist = _phase1d_labels()
    rows, screens, all_records = [], {}, {}
    mp.mp.dps = DPS
    for mat in MATERIALS:
        for theta in ANGLES:
            records, screen = build_branch_records(mat, theta)
            screens[f"{mat}@{theta}"] = screen
            all_records[f"{mat}@{theta}"] = records
            for rec in records:
                if rec is None:
                    continue
                md = rec.classification
                br = rec.branch
                old = hist.get((mat, float(theta), br), "ABSENT")
                changed = old not in ("ABSENT", md["spatial_class"]) and not (
                    old == "PROPAGATING"
                    and md["spatial_class"] == "WEAKLY_DAMPED"
                    and md["attenuation_ratio"] is not None
                    and float(md["attenuation_ratio"]) < 1e-6)
                rows.append({
                    "experiment_id": EXPERIMENT_ID,
                    "material": mat,
                    "theta_deg": f"{theta:g}",
                    "branch": br,
                    # taken from THIS record only
                    "Re_xi": _f(mp.re(rec.xi)),
                    "Im_xi": _f(mp.im(rec.xi)),
                    "abs_xi": _f(abs(rec.xi)),
                    "attenuation_ratio": _f(md["attenuation_ratio"]),
                    "attenuation_basis": md["attenuation_basis"],
                    "decay_length_m": _f(md["decay_length_m"]),
                    "decay_length_in_transverse_wavelengths":
                        _f(md["decay_length_in_transverse_wavelengths"]),
                    "spatial_class": md["spatial_class"],
                    "radiation_dir": md["radiation_dir"],
                    "radiation_source": md["radiation_source"],
                    "energy_flux_status": md["energy_flux_status"],
                    "confidence": md["confidence"],
                    "mech_content": _f(md["mech_content"]),
                    "elec_content": _f(md["elec_content"]),
                    "thermal_content": _f(md["thermal_content"]),
                    "thermal_dispersion_ratio": _f(md["thermal_dispersion_ratio"]),
                    "J_z": _f(rec.jz),
                    "J_z_sign": ("+" if rec.jz > 0 else
                                 "-" if rec.jz < 0 else "0"),
                    "normalisation_basis": rec.nullspace["normalisation_basis"],
                    "nullspace_agreement": rec.nullspace["agreement"],
                    "cramer_gap": _f(rec.nullspace["cramer_gap"]),
                    "gram_residual": _f(rec.nullspace.get("gram_residual")),
                    "structural_columns":
                        "|".join(map(str, rec.nullspace["structural_columns"])),
                    "admissible_side": rec.admissibility,
                    "outgoing_status": rec.outgoing_status,
                    "phase_1d_label": old,
                    "label_changed": changed,
                    "reasons": "; ".join(md["reasons"]),
                })
    return rows, screens, all_records


def main():
    rows, screens, _records = build_rows()
    out_csv = ROOT / "validation" / "phase_02a_modal_inventory_v2.csv"
    with out_csv.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLUMNS)
        wr.writeheader()
        wr.writerows(rows)

    changed = [r for r in rows if r["label_changed"] is True]
    print(f"{EXPERIMENT_ID}  (supersedes {SUPERSEDES})")
    print(f"  rows written          : {len(rows)}  -> {out_csv.relative_to(ROOT)}")
    print(f"  labels changed vs 1D  : {len(changed)}")
    print(f"  PROPAGATING labels    : "
          f"{sum(1 for r in rows if r['spatial_class'] == 'PROPAGATING')}")
    for cls in ("EVANESCENT", "WEAKLY_DAMPED", "THERMAL_DIFFUSIVE", "GAUGE"):
        print(f"    {cls:18s}: "
              f"{sum(1 for r in rows if r['spatial_class'] == cls)}")
    manifest = {
        "experiment_id": EXPERIMENT_ID,
        "supersedes": SUPERSEDES,
        "supersession_reason":
            "P0-1: exported Re_xi, Im_xi, abs_xi and J_z were read from a loop "
            "variable after iteration, so every row in a block carried the last "
            "branch's values. Artifact is INVALID as scientific evidence.",
        "commit": _commit(),
        "working_precision_dps": DPS,
        "frequency_hz": F_HZ,
        "materials": list(MATERIALS),
        "angles_deg": list(ANGLES),
        "rows": len(rows),
        "labels_changed_vs_phase_1d": len(changed),
        "output": str(out_csv.relative_to(ROOT)),
        "output_sha256": hashlib.sha256(out_csv.read_bytes()).hexdigest(),
        "columns": COLUMNS,
        "preserved_evidence": [
            "validation/phase_01d_independent_comparison.csv",
            "validation/phase_02a_modal_inventory.csv",
        ],
        "halfspace_screens": screens,
        "reproduce": "cd code && python3 -m scattering.export_inventory",
    }
    mp_path = ROOT / "validation" / "phase_02a_modal_inventory_v2_manifest.json"
    mp_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"  manifest             : {mp_path.relative_to(ROOT)}")
    print(f"  sha256               : {manifest['output_sha256'][:16]}...")


if __name__ == "__main__":
    main()
