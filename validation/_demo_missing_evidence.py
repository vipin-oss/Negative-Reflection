#!/usr/bin/env python3
"""Demo 6: deleting a mandatory evidence artifact must prevent promotion.

Runs in an isolated copy. The real repository is never modified.
"""
import shutil, subprocess, sys, tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
IGN = shutil.ignore_patterns(".git", "__pycache__", "*.pyc", ".pytest_cache")
tmp = Path(tempfile.mkdtemp(prefix="phase1c_demo6_"))
dst = tmp / "repo"
try:
    shutil.copytree(REPO, dst, ignore=IGN)
    victim = dst / "evidence" / "equation_registry.csv"
    before = victim.stat().st_size
    victim.unlink()
    print(f"deleted evidence/equation_registry.csv ({before} bytes) in {dst}")
    r = subprocess.run([sys.executable, "validation/run_fast.py", "--promotion"],
                       cwd=dst, capture_output=True, text=True, timeout=1800)
    print(r.stdout[-3000:])
    print("exit:", r.returncode)
    if r.returncode == 0:
        print("PROMOTION WAS ALLOWED DESPITE MISSING EVIDENCE -- this is a defect")
    sys.exit(r.returncode)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
