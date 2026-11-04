#!/usr/bin/env python3
"""Record SHA-256 of every frozen scientific input.

Run deliberately when an input legitimately changes; the change must then be
logged in evidence/approved_baselines.json. Nothing else may rewrite this file.
"""
import hashlib, json, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PATTERNS = ["literature/material_constants.csv", "literature/verified_references.csv",
            "code/det_coeffs.json", "theory/characteristic_matrix.tex",
            "theory/energy_flux_derivation.tex", "theory/governing_equations.tex",
            "theory/interface_conditions.tex", "validation/schemas/material_constants.schema.json"]
OUT = REPO / "evidence" / "input_checksums.json"

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    try:
        commit = subprocess.run(["git","rev-parse","HEAD"],cwd=REPO,capture_output=True,
                                text=True).stdout.strip()
    except Exception:
        commit = "UNKNOWN"
    files = {}
    for rel in PATTERNS:
        p = REPO / rel
        if not p.exists():
            print(f"MISSING {rel}", file=sys.stderr); return 1
        files[rel] = sha(p)
    doc = dict(schema_version="1.0",
               generated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
               commit=commit,
               policy="Changing any of these files requires regenerating this register "
                      "AND recording the reason in evidence/approved_baselines.json.",
               files=files)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(REPO)} with {len(files)} entries")
    return 0

if __name__ == "__main__":
    sys.exit(main())
