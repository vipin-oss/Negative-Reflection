#!/usr/bin/env python3
"""LEVEL 3 -- FULL. Runs before a milestone or a submission package.

Adds to Level 2:
  * every accepted experiment re-run from scratch with its manifest command
  * parameter sensitivity sweeps (frequency, angle, tau_q, kappa)
  * a precision study (float64 vs dps 30/50/80)
  * benchmark reproduction against every external reference
  * figure regeneration from raw data plus the plotting script, with a pixel/
    data comparison against the tracked figure
  * claim-to-evidence consistency across registries, graph and manifests

Figures are only checked here, not in Level 1 or 2, because regeneration is slow
and depends on matplotlib.

    python3 validation/run_full.py
    python3 validation/run_full.py --promotion
    python3 validation/run_full.py --with-figures     # also regenerate figures

Exit codes: 0 all pass, 1 failure, 2 unresolved/unknown, 3 nothing executed.
"""
import json
import subprocess
import sys
from pathlib import Path

from _runner import run, REPO, HERE

if __name__ == "__main__":
    with_figures = "--with-figures" in sys.argv
    if with_figures:
        sys.argv.remove("--with-figures")

    print("=" * 78)
    print("LEVEL 3: regenerating every accepted experiment manifest")
    print("=" * 78)
    regen = subprocess.run([sys.executable, str(HERE / "regen_manifests.py")],
                           cwd=REPO, capture_output=True, text=True)
    print(regen.stdout[-3000:], end="")
    if regen.returncode != 0:
        print(regen.stderr[-2000:], file=sys.stderr)
        print("MANIFEST REGENERATION FAILED", file=sys.stderr)
        sys.exit(1)

    if with_figures:
        fig = subprocess.run([sys.executable, str(HERE / "figures" / "regenerate.py")],
                             cwd=REPO, capture_output=True, text=True)
        print(fig.stdout[-3000:], end="")
        if fig.returncode != 0:
            print("FIGURE REGENERATION FAILED", file=sys.stderr)
            sys.exit(1)

    print()
    sys.exit(run("3", "FULL",
                 description="milestone gate: adds manifest regeneration, parameter "
                             "sensitivity, precision study, benchmark reproduction, "
                             "figure regeneration, claim-to-evidence consistency"))
