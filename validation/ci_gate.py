#!/usr/bin/env python3
"""CI gate: fail on any contradiction or any UNREGISTERED unknown.

Why this exists. `run_fast.py` deliberately exits 2 ("something is unknown") which
is neither success nor failure. A CI step treats every non-zero exit as failure,
so without this gate the workflow would be permanently red for a reason that is
honest but not actionable - and a permanently red build is one nobody reads.

The wrong fix is to make the suite exit 0 when it finds unknowns. The right fix is
an enumerated register (validation/known_open.json) that:

  * names every tolerated case, its reason and the blocker that owns it;
  * still reports each one as UNRESOLVED in the machine-readable report;
  * fails the build on any FAIL;
  * fails the build on any unknown NOT in the register;
  * fails the build on a register entry whose blocker has been closed but which
    was not removed.

So an unknown can be tolerated by decision, never by accident.

    python3 validation/ci_gate.py validation/evidence/last_run_1.json

Exit: 0 gate satisfied, 1 real failure, 2 unregistered unknown (also a failure).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
REGISTER = HERE / "known_open.json"
BLOCKERS = REPO / "reports" / "PHASE_01C_OPEN_BLOCKERS.json"


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    rep_path = Path(argv[0]) if argv else HERE / "evidence" / "last_run_1.json"
    if not rep_path.exists():
        print(f"NO REPORT at {rep_path}: the suite did not run. This is a failure.")
        return 1
    rep = json.loads(rep_path.read_text())
    reg = json.loads(REGISTER.read_text())
    allowed = {e["test_id"]: e for e in reg["entries"]}
    open_blockers = {b["id"] for b in json.loads(BLOCKERS.read_text())["blockers"]}

    cases = rep["cases"]
    fails = [c for c in cases if c["status"] == "FAIL"]
    others = [c for c in cases if c["status"] in ("UNRESOLVED", "WITHDRAWN")]

    print("=" * 78)
    print("CI GATE")
    print("=" * 78)
    print(f"report: {rep_path}")
    print(f"commit: {rep['commit']}")
    s = rep["summary"]
    print(f"{s['PASS']} PASS, {s['FAIL']} FAIL, {s['UNRESOLVED']} UNRESOLVED, "
          f"{s['WITHDRAWN']} WITHDRAWN, {s['NOT_APPLICABLE']} NOT_APPLICABLE, "
          f"{s['total']} total; {s['executed_cases']} executed cases")

    problems = []

    # 1. a contradiction is never acceptable
    for c in fails:
        problems.append(f"FAIL {c['test_id']}: {c['purpose'][:90]}")

    # 2. unknown is acceptable only if named in advance
    for c in others:
        e = allowed.get(c["test_id"])
        if e is None:
            problems.append(f"UNREGISTERED {c['status']} {c['test_id']}: "
                            f"{c['purpose'][:90]}")
        elif e["status"] != c["status"]:
            problems.append(f"{c['test_id']} reported {c['status']} but the register "
                            f"expects {e['status']}")

    # 3. the register must not outlive the blockers that justify it
    reported = {c["test_id"]: c["status"] for c in cases}
    for tid, e in allowed.items():
        if e["blocker"] not in open_blockers:
            problems.append(f"STALE REGISTER ENTRY {tid}: its blocker {e['blocker']} "
                            f"is closed, so it must be removed from known_open.json")
        if tid not in reported:
            problems.append(f"REGISTERED TEST NOT RUN {tid}: a tolerated unknown that "
                            f"did not execute is not tolerated, it is invisible")

    if allowed:
        print("\nTolerated unknowns (NOT passes; reported as such everywhere):")
        for tid, e in allowed.items():
            print(f"  {tid}: {e['status']} - {e['reason'][:80]}")
            print(f"        owned by {e['blocker']}, closes when: {e['closes_when'][:70]}")

    print()
    if problems:
        print("GATE NOT SATISFIED:")
        for p in problems:
            print("  - " + p)
        if fails:
            return 1
        return 2
    print("GATE SATISFIED: no contradictions, and every unknown is a registered, "
          "owned, open item.")
    print("This is NOT a statement that the research claims are accepted. "
          "Promotion still requires all blockers closed and supervisor authorisation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
