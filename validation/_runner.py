"""Shared logic for validation/run_{fast,standard,full}.py.

Exit codes are deliberately split. A run that could not be executed and a run
that executed and found a contradiction are different scientific events and must
not share an exit code:

    0  every case PASS
    1  at least one case FAIL         -- something is WRONG
    2  no failures, but at least one UNRESOLVED / WITHDRAWN / NOT_APPLICABLE
       -- something is UNKNOWN; this must never be reported as success
    3  infrastructure problem (pytest itself crashed, no cases recorded)

`--promotion` adds a fourth outcome, because the brief requires promotion to be
blocked while a mandatory blocker is open even when no test fails:

    exits 1 if promotion_allowed is False
"""
from __future__ import annotations

import argparse
import shutil
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "code"))

import harness  # noqa: E402


def run(level: str, level_name: str, markers: str | None = None,
        description: str = "") -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--promotion", action="store_true",
                    help="exit 1 unless promotion is allowed (blockers closed)")
    ap.add_argument("--json-out", default=None)
    a = ap.parse_args()

    # pytest must run in a SUBPROCESS (a crashed test must not take the runner
    # down), so the records cannot be read from this process's memory: the cases
    # live in the child's `harness.RECORDS`. The conftest hook serialises them to
    # PHASE1C_REPORT, and the parent loads that file. Getting this wrong made the
    # first version of this runner report "no cases recorded" for a suite that had
    # in fact run - a silent false negative in the reporting layer itself.
    import os
    import tempfile
    tmpdir = Path(tempfile.mkdtemp(prefix="phase1c_run_"))
    raw_report = tmpdir / "cases.json"
    env = dict(os.environ, PHASE1C_REPORT=str(raw_report))

    cmd = [sys.executable, "-m", "pytest", str(HERE / "tests"), "-q",
           "--no-header", "-p", "no:cacheprovider"]
    if markers:
        cmd += ["-m", markers]
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True, env=env)
    dt = time.time() - t0

    cases = []
    if raw_report.exists():
        import json as _json
        blob = _json.loads(raw_report.read_text())
        for d in blob.get("cases", []):
            c = harness.Case(**{k: d[k] for k in (
                "test_id", "category", "purpose", "oracle", "oracle_independent",
                "fixture", "expected", "tolerance", "tol_justification", "level",
                "status", "n_exec", "actual", "detail", "duration_s", "commit",
                "remote_verified")})
            cases.append(c)
        shutil_rm = __import__("shutil")
        shutil_rm.rmtree(tmpdir, ignore_errors=True)

    if not cases:
        # pytest crashed before any case recorded its evidence
        sys.stderr.write(proc.stdout[-4000:] + "\n" + proc.stderr[-4000:])
        print("NO CASES RECORDED: the suite did not execute. Treating as "
              "infrastructure failure, not as success.", file=sys.stderr)
        return 3

    out = Path(a.json_out) if a.json_out else HERE / "evidence" / f"last_run_{level}.json"
    rep = harness.write_report(cases, out, level=level)
    rep["level_name"] = level_name
    rep["description"] = description
    rep["wall_clock_s"] = round(dt, 2)
    rep["pytest_exit"] = proc.returncode
    with open(out, "w") as fh:
        json.dump(rep, fh, indent=2)

    harness.print_summary(rep)

    s = rep["summary"]
    if s["FAIL"]:
        code = 1
    elif s["UNRESOLVED"] or s["WITHDRAWN"] or s["NOT_APPLICABLE"]:
        code = 2
    else:
        code = 0

    if a.promotion:
        if not rep["promotion_allowed"]:
            print("\nPROMOTION GATE: REFUSED")
            for r in rep["promotion_blocked_reasons"]:
                print(f"  - {r}")
            return 1
        print("\nPROMOTION GATE: satisfied by automation. Supervisor "
              "authorisation is still required and is NOT implied by this exit code.")
        return 0

    print(f"\nreport: {out.relative_to(REPO)}")
    print(f"exit {code}  (0=all pass, 1=failure, 2=unresolved/unknown, 3=no cases)")
    return code
