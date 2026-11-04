#!/usr/bin/env python3
"""LEVEL 2 -- STANDARD. Runs before a phase submission.

Adds to Level 1:
  * the full Phase 1B foundation suite (16 checks, A1-J2) via subprocess
  * expanded symbolic work: null space, block structure, reduction chain
  * root cases and precision convergence
  * interface rank measurement on three material pairs
  * energy identity and passivity with a negative control
  * the mutation suite (proves the tests can detect defects)
  * manifest integrity for every accepted experiment

Typical wall clock 10-20 minutes on one core (the 4x4 symbolic determinants
dominate).

    python3 validation/run_standard.py
    python3 validation/run_standard.py --promotion

Exit codes: 0 all pass, 1 failure, 2 unresolved/unknown, 3 nothing executed.
"""
import subprocess
import sys
from pathlib import Path

from _runner import run, REPO, HERE

if __name__ == "__main__":
    print("=" * 78)
    print("LEVEL 2 also runs the mutation suite in an isolated temporary tree.")
    print("=" * 78)
    mut = subprocess.run([sys.executable, str(HERE / "mutation_tests" / "run_mutations.py")],
                         cwd=REPO, capture_output=True, text=True)
    p = HERE / "mutation_tests" / "last_mutation_run.json"
    if p.exists():
        import json
        j = json.loads(p.read_text())
        print(f"mutations: {j['detected']}/{j['total']} detected, "
              f"survivors: {j['survivors'] or 'none'}")
    else:
        print("mutations: no report produced", file=sys.stderr)
    print()
    code = run("2", "STANDARD",
               description="pre-submission gate: adds the full foundation suite, "
                           "null-space and block structure, convergence, interface "
                           "rank, energy with negative control, mutation suite")
    if not p.exists():
        return_code = max(code, 1)
    else:
        return_code = code if j["survivors"] == [] else max(code, 1)
    sys.exit(return_code)
