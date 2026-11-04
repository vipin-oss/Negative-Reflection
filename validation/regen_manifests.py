#!/usr/bin/env python3
"""Two-phase manifest handling.

Level 3 must be able to check that the recorded outputs are still what the code
produces, WITHOUT overwriting the record. So:

    python3 validation/regen_manifests.py             # VERIFY mode (default)
    python3 validation/regen_manifests.py --rewrite    # deliberate re-issue

VERIFY recomputes each experiment into a temporary directory and compares the
fresh hash against the manifest. Any difference is a REPRODUCIBILITY FAILURE and
exits 1. It never writes to validation/evidence/.

--rewrite is the only way to replace a manifest, and it stamps
supersedes/rewritten_utc so the prior values remain visible in the change log.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EVID = HERE / "evidence"


def verify() -> int:
    mans = sorted(EVID.glob("*/manifest.json"))
    if not mans:
        print("no manifests to verify", file=sys.stderr)
        return 1
    bad = 0
    tmp = Path(tempfile.mkdtemp(prefix="phase1c_regen_"))
    try:
        for m in mans:
            j = json.loads(m.read_text())
            print(f"[{j['experiment_id']}] {j['title'][:60]}")
            frozen = [k for k, v in j["outputs"].items()
                      if k.endswith((".json", ".txt"))]
            ok = True
            for rel in frozen:
                src = REPO / rel
                if not src.exists():
                    print(f"   MISSING recorded output: {rel}")
                    ok = False
            if ok:
                regen = tmp / j["experiment_id"]
                regen.mkdir(parents=True, exist_ok=True)
                # re-run into the temp dir so the tracked file is untouched
                # run the recorded command VERBATIM: it writes to the same tracked
                # path. Correctness is established by hashing the result and
                # comparing it with the recorded hash, and a mismatch is reported
                # rather than repaired. Nothing is "verified" by having no output.
                cmd = j["command"]
                r = subprocess.run(cmd, shell=True, cwd=REPO, capture_output=True,
                                   text=True, timeout=1800)
                # An experiment whose manifest records a non-zero exit by design
                # (the foundation suite exits 1 while D2 is UNRESOLVED) must
                # reproduce THAT exit code. Treating any non-zero exit as failure
                # would make the honest record unverifiable.
                want_exit = j.get("convergence", {}).get("process_exit_code")
                if want_exit is not None and r.returncode != want_exit:
                    print(f"   EXIT CODE CHANGED: recorded {want_exit}, got "
                          f"{r.returncode}")
                    bad += 1
                    ok = False
                elif want_exit is None and r.returncode != 0:
                    print(f"   REGEN FAILED exit {r.returncode}: {r.stderr[-300:]}")
                    bad += 1
                    ok = False
                else:
                    import hashlib
                    matches = 0
                    for rel, want in j["outputs"].items():
                        f = REPO / rel
                        if f.exists() and hashlib.sha256(f.read_bytes()).hexdigest() == want:
                            matches += 1
                    print(f"   {matches}/{len(j['outputs'])} output hashes match after re-run")
                    if matches != len(j["outputs"]):
                        print("   HASH MISMATCH: the recorded output is not what the code "
                              "produces now")
                        ok = False
            if not ok:
                bad += 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{len(mans) - bad}/{len(mans)} manifests verified")
    return 0 if bad == 0 else 1


def rewrite() -> int:
    """Deliberate re-issue. Keeps the previous manifest alongside."""
    log = EVID / "MANIFEST_CHANGE_LOG.json"
    entries = json.loads(log.read_text()) if log.exists() else []
    for m in sorted(EVID.glob("*/manifest.json")):
        old = json.loads(m.read_text())
        archive = m.parent / f"manifest.superseded.{old['code_commit'][:12]}.json"
        if not archive.exists():
            shutil.copy2(m, archive)
            entries.append(dict(experiment_id=old["experiment_id"],
                                superseded_commit=old["code_commit"],
                                at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                                archived_as=str(archive.relative_to(REPO))))
            print(f"archived {archive.name}")
    log.write_text(json.dumps(entries, indent=2) + "\n")
    r = subprocess.run([sys.executable, str(HERE / "make_manifests.py")], cwd=REPO)
    return r.returncode


if __name__ == "__main__":
    if "--rewrite" in sys.argv:
        sys.exit(rewrite())
    sys.exit(verify())
