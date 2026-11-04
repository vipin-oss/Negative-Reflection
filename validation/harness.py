"""
validation/harness.py -- Phase 1C test harness.

Every test in validation/tests/ is built with `make_case`, which enforces the
Phase 1C evidence contract (task D). A test that does not declare all of:

    test_id         stable identifier
    purpose         scientific purpose, one sentence
    oracle          the exact reference used, and whether it is INDEPENDENT
    fixture         the input fixture description
    expected        the expected result
    tolerance       numeric tolerance or None
    tol_justification   why that tolerance, in words

...cannot be constructed. This is deliberate: the Phase 1B lesson was that
claims stronger than their evidence are the main scientific risk, so the harness
makes an un-evidenced test impossible to write by accident.

Statuses
--------
PASS  FAIL  UNRESOLVED  WITHDRAWN  NOT_APPLICABLE

Rules (all enforced in code, not by convention):

* a case with n_exec == 0 can never be PASS or NOT_APPLICABLE-as-success;
* FAIL / UNRESOLVED / WITHDRAWN are never collapsed into PASS;
* the report carries the source commit for every case;
* an exception inside a case is a FAIL, never a skip.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass, field, asdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "code"))

PASS = "PASS"
FAIL = "FAIL"
UNRESOLVED = "UNRESOLVED"
WITHDRAWN = "WITHDRAWN"
NOT_APPLICABLE = "NOT_APPLICABLE"

VALID = (PASS, FAIL, UNRESOLVED, WITHDRAWN, NOT_APPLICABLE)
#: statuses that mean "we have no evidence this is right"
NON_EVIDENCE = (FAIL, UNRESOLVED, WITHDRAWN, NOT_APPLICABLE)


# ---------------------------------------------------------------------------
def source_commit() -> str:
    """HEAD SHA, or 'UNKNOWN' if git is unavailable (must not silently pass)."""
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                             capture_output=True, text=True, timeout=20)
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:                                          # noqa: BLE001
        pass
    return "UNKNOWN"


ABSENT = "ABSENT"


@dataclass
class Case:
    """One executed test case with its full evidence record."""
    test_id: str
    category: str
    purpose: str
    oracle: str
    oracle_independent: bool          # True / False / ABSENT
    fixture: str
    expected: str
    tolerance: float | None
    tol_justification: str
    level: str                        # "1" | "2" | "3"
    status: str = UNRESOLVED
    n_exec: int = 0
    actual: str = ""
    detail: str = ""
    duration_s: float = 0.0
    commit: str = ""
    remote_verified: bool = False

    def as_dict(self) -> dict:
        return asdict(self)


def make_case(**kw) -> Case:
    """Construct a Case, refusing incomplete declarations."""
    required = ("test_id", "category", "purpose", "oracle",
                "oracle_independent", "fixture", "expected",
                "tolerance", "tol_justification", "level")
    missing = [r for r in required if r not in kw]
    if missing:
        raise ValueError(
            f"test case is missing mandatory evidence fields: {missing}")
    if kw["oracle_independent"] is True and not kw.get("oracle"):
        raise ValueError("oracle_independent=True requires a stated oracle")
    if kw["tolerance"] is None and not kw.get("tol_justification"):
        raise ValueError(
            "a case with tolerance=None must justify why (tol_justification)")
    return Case(**kw)


#: every executed case, in order; the conftest hook serialises this
RECORDS: list[Case] = []


def finish(case: Case, fn) -> None:
    """Execute a case, record it, and turn its status into a pytest verdict.

    pytest must agree with the harness: a case that is not PASS makes the test
    fail, so the two can never diverge. The full record still goes to the JSON
    report, which is what carries the evidence.
    """
    execute(case, fn)
    RECORDS.append(case)
    if case.status != PASS:
        raise AssertionError(
            f"[{case.status}] {case.test_id} ({case.category}): {case.purpose}\n"
            f"  oracle:   {case.oracle}\n"
            f"  fixture:  {case.fixture}\n"
            f"  expected: {case.expected}\n"
            f"  actual:   {case.actual}\n"
            f"  detail:   {case.detail}")


def execute(case: Case, fn) -> Case:
    """Run fn(); fn returns (status, n_exec, actual, detail)."""
    t0 = time.time()
    case.commit = source_commit()
    try:
        got = fn()
        if len(got) == 4:
            status, n_exec, actual, detail = got
        elif len(got) == 3:
            status, n_exec, actual = got
            detail = ""
        else:
            raise ValueError(
                f"{case.test_id}: test body must return 3 or 4 values, got {got!r}")
        case.status, case.n_exec, case.actual, case.detail = (
            status, int(n_exec), str(actual), str(detail))
    except Exception as exc:                                   # noqa: BLE001
        case.status = FAIL
        case.n_exec = 0
        case.actual = f"exception {exc!r}"
        case.detail = "".join(traceback.format_exc().splitlines(True)[-4:])
    case.duration_s = time.time() - t0

    # ---- enforcement -----------------------------------------------------
    if case.status not in VALID:
        case.status = FAIL
        case.detail = f"invalid status {case.status!r}\n" + case.detail
    if case.n_exec == 0 and case.status in (PASS, NOT_APPLICABLE):
        case.status = UNRESOLVED
        case.detail = ("FORCED UNRESOLVED: 0 cases executed, so a PASS or "
                       "NOT_APPLICABLE would be unevidenced.\n" + case.detail)
    return case


# ---------------------------------------------------------------------------
BLOCKERS = REPO / "reports" / "PHASE_01C_OPEN_BLOCKERS.json"


def open_blockers() -> list[dict]:
    """Mandatory P0 blockers, read from the tracked blocker file."""
    if not BLOCKERS.exists():
        return [dict(id="BLOCKER-FILE-MISSING", severity="P0",
                     title="blocker register missing", promotion_blocked=True)]
    with open(BLOCKERS) as fh:
        return json.load(fh).get("blockers", [])


def promotion_allowed(cases: list[Case]) -> tuple[bool, list[str]]:
    """A run may be promoted only if no case is FAIL/UNRESOLVED/WITHDRAWN
    *and* no mandatory P0 blocker remains open."""
    reasons = []
    bad = [c for c in cases if c.status in NON_EVIDENCE]
    if bad:
        reasons.append(
            f"{len(bad)} case(s) not passing: "
            + ", ".join(f"{c.test_id}={c.status}" for c in bad))
    for b in open_blockers():
        if b.get("promotion_blocked"):
            reasons.append(f"blocker {b['id']} open ({b.get('severity')})")
    return (not reasons), reasons


def summarise(cases: list[Case]) -> dict:
    from collections import Counter
    cnt = Counter(c.status for c in cases)
    return dict(
        total=len(cases),
        PASS=cnt.get(PASS, 0),
        FAIL=cnt.get(FAIL, 0),
        UNRESOLVED=cnt.get(UNRESOLVED, 0),
        WITHDRAWN=cnt.get(WITHDRAWN, 0),
        NOT_APPLICABLE=cnt.get(NOT_APPLICABLE, 0),
        executed_cases=sum(c.n_exec for c in cases),
        zero_case_tests=[c.test_id for c in cases if c.n_exec == 0],
        commit=source_commit(),
        by_category=_by_category(cases),
    )


def _by_category(cases):
    out = {}
    for c in cases:
        d = out.setdefault(c.category, {PASS: 0, FAIL: 0, UNRESOLVED: 0,
                                        WITHDRAWN: 0, NOT_APPLICABLE: 0})
        d[c.status] += 1
    return out


def write_report(cases: list[Case], path: Path, level: str,
                 extra: dict | None = None) -> dict:
    ok, reasons = promotion_allowed(cases)
    rep = dict(
        schema_version="1.0",
        level=level,
        generated_unix=int(time.time()),
        commit=source_commit(),
        python=sys.version.split()[0],
        summary=summarise(cases),
        promotion_allowed=ok,
        promotion_blocked_reasons=reasons,
        blockers=open_blockers(),
        cases=[c.as_dict() for c in cases],
    )
    if extra:
        rep.update(extra)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=2)
    return rep


def print_summary(rep: dict, out=sys.stdout) -> None:
    s = rep["summary"]
    print("=" * 78, file=out)
    print(f"LEVEL {rep['level']}  commit {rep['commit'][:12]}  "
          f"python {rep['python']}", file=out)
    print("=" * 78, file=out)
    print(f"{s['PASS']} PASS, {s['FAIL']} FAIL, {s['UNRESOLVED']} UNRESOLVED, "
          f"{s['WITHDRAWN']} WITHDRAWN, {s['NOT_APPLICABLE']} NOT_APPLICABLE, "
          f"{s['total']} total", file=out)
    print(f"executed cases: {s['executed_cases']}", file=out)
    if s["zero_case_tests"]:
        print(f"ZERO-CASE TESTS: {s['zero_case_tests']}", file=out)
    for cat, d in sorted(s["by_category"].items()):
        tot = sum(d.values())
        print(f"  {cat:42s} {d[PASS]}/{tot} PASS", file=out)
    for c in rep["cases"]:
        if c["status"] != PASS:
            print(f"  {c['status']:15s} {c['test_id']:10s} {c['purpose'][:60]}",
                  file=out)
    print("-" * 78, file=out)
    if rep["promotion_allowed"]:
        print("PROMOTION ALLOWED (necessary, NOT sufficient - supervisor "
              "authorisation still required)", file=out)
    else:
        print("PROMOTION BLOCKED:", file=out)
        for r in rep["promotion_blocked_reasons"]:
            print(f"  - {r}", file=out)
    print("=" * 78, file=out)
