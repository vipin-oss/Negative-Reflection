import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent.parent / "code"))

import harness  # noqa: E402


def pytest_sessionfinish(session, exitstatus):
    """Write the machine-readable evidence report for whatever ran."""
    if not harness.RECORDS:
        return
    out = Path(__import__("os").environ.get(
        "PHASE1C_REPORT", HERE.parent / "evidence" / "last_run.json"))
    harness.write_report(harness.RECORDS, out, level="pytest")
