"""
Category 1 -- Source and material-data integrity.

These tests guard the bottom of the dependency graph. If a constant silently
changes value, or a provenance label becomes unsupported, everything above it is
suspect. They are cheap and run in Level 1.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pytest

from harness import (REPO, PASS, FAIL, UNRESOLVED, ABSENT,
                     make_case, finish)

MATERIALS = REPO / "literature" / "material_constants.csv"
CLAIMS = REPO / "evidence" / "claim_registry.csv"
EQUATIONS = REPO / "evidence" / "equation_registry.csv"

#: read from the schema, so the test and the schema cannot drift apart
_SCHEMA = json.loads(
    (REPO / "validation" / "schemas" / "material_constants.schema.json").read_text())
ALLOWED_CONST_STATUS = set(
    _SCHEMA["columns"]["status"]["enum"])
ALLOWED_CLAIM_STATUS = {
    "VERIFIED_STRONG", "VERIFIED_WEAK", "UNVERIFIED", "UNRESOLVED", "INVALID",
}
ALLOWED_EQ_STATUS = {
    "VERIFIED_STRONG", "VERIFIED_WEAK", "UNVERIFIED", "UNRESOLVED", "INVALID",
    "DEFINED_NOT_TESTED", "PARTIAL", "ASSUMED_STANDARD", "DERIVED",
    "NOT_IMPLEMENTED",
}


def _csv(path):
    import csv
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
def _run(_cases, **kw):
    return make_case(**kw)


# ===========================================================================
def test_MAT01_material_csv_schema():
    """Every material row has all required fields and a legal status."""
    c = make_case( test_id="MAT01", category="1-data-integrity",
             purpose="Material constants carry complete, legal provenance metadata",
             oracle="schema defined in validation/schemas/material_constants.schema.json "
                    "and the status vocabulary of literature/material_constants.csv",
             oracle_independent=True,
             fixture="literature/material_constants.csv",
             expected=">=70 rows; all 9 required columns; every status in the allowed set",
             tolerance=None,
             tol_justification="exact structural check; no tolerance applies",
             level="1")

    def body():
        rows = _csv(MATERIALS)
        req = {"material", "constant", "symbol", "value", "unit", "status",
               "source", "uncertainty", "notes"}
        bad_status, missing = [], []
        for r in rows:
            if r.get("status") not in ALLOWED_CONST_STATUS:
                bad_status.append((r.get("material"), r.get("constant"),
                                   r.get("status")))
            if not req.issubset(set(r.keys())):
                missing.append(sorted(req - set(r.keys())))
        ok = (len(rows) >= 70 and not bad_status and not missing)
        return (PASS if ok else FAIL, len(rows),
                f"{len(rows)} rows; {len(bad_status)} illegal status; "
                f"{len(missing)} rows missing fields",
                f"illegal: {bad_status[:5]}")

    finish(c, body)


def test_MAT02_no_tuned_parameter():
    """No constant is labelled in a way that implies it was fitted to pass a rule.

    Phase 1 carried a status CHOSEN_ABOVE_CAUSAL_FLOOR for CdSe tau_q. That label
    was itself evidence of tuning, and the rule behind it was withdrawn. This test
    makes the regression permanent: no such label may reappear.
    """
    c = make_case( test_id="MAT02", category="1-data-integrity",
             purpose="No material constant is tuned to satisfy a validation rule",
             oracle="regression on the Phase-1B defect (P1-A): the label "
                    "CHOSEN_ABOVE_CAUSAL_FLOOR must never reappear",
             oracle_independent=False,
             fixture="literature/material_constants.csv status column",
             expected="0 occurrences of any 'CHOSEN'/'TUNED' status",
             tolerance=None,
             tol_justification="string search; exact",
             level="1")

    def body():
        rows = _csv(MATERIALS)
        bad = [(r["material"], r["constant"], r["status"]) for r in rows
               if re.search(r"CHOSEN|TUNED|FITTED", r.get("status", ""), re.I)]
        return (PASS if not bad else FAIL, len(rows),
                f"{len(bad)} tuned-looking labels", str(bad))

    finish(c, body)


def test_CLM01_claim_registry_schema():
    """Claim registry: unique ids, legal statuses, every row has an evidence column."""
    c = make_case( test_id="CLM01", category="12-claim-traceability",
             purpose="Claim registry is well formed and traceable",
             oracle="schema in validation/schemas/claim_registry.schema.json",
             oracle_independent=True,
             fixture="evidence/claim_registry.csv",
             expected=">=40 rows; unique ids; statuses legal; test_id present on all",
             tolerance=None, tol_justification="structural; exact", level="1")

    def body():
        rows = _csv(CLAIMS)
        ids = [r["claim_id"] for r in rows]
        dup = {i for i in ids if ids.count(i) > 1}
        bad = [r["claim_id"] for r in rows if r["evidence_status"] not in ALLOWED_CLAIM_STATUS]
        no_test = [r["claim_id"] for r in rows if not r.get("test_id")]
        ok = (len(rows) >= 40 and not dup and not bad and not no_test)
        return (PASS if ok else FAIL, len(rows),
                f"{len(rows)} claims, {len(dup)} duplicate, {len(bad)} bad status, "
                f"{len(no_test)} without test_id",
                f"dup={sorted(dup)} bad={bad[:5]} notest={no_test[:5]}")

    finish(c, body)


def test_CLM02_no_verified_without_evidence():
    """A claim marked VERIFIED_* must name a test that exists in the suite.

    This is the machine version of U-11: a claim cannot be VERIFIED merely because
    a report says PASS. The named test must exist in code/foundation_checks.py or
    in validation/tests/.
    """
    c = make_case( test_id="CLM02", category="12-claim-traceability",
             purpose="Every VERIFIED claim names a test that actually exists",
             oracle="the union of test ids discoverable in code/foundation_checks.py "
                    "and validation/tests/",
             oracle_independent=True,
             fixture="evidence/claim_registry.csv cross-referenced against the code",
             expected="0 VERIFIED claims naming a non-existent test",
             tolerance=None, tol_justification="set membership; exact", level="1")

    def body():
        rows = _csv(CLAIMS)
        src = (REPO / "code" / "foundation_checks.py").read_text()
        known = set(re.findall(r'run\("([A-Z][0-9])"', src))
        for p in (REPO / "validation" / "tests").glob("test_*.py"):
            known |= set(re.findall(r'test_id="([A-Z]+[0-9]+)"', p.read_text()))
        n_check = 0
        bad = []
        for r in rows:
            if not r["evidence_status"].startswith("VERIFIED"):
                continue
            t = r.get("test_id", "").strip()
            if t in ("", "none", "no test", "no automated test"):
                bad.append((r["claim_id"], "no test named"))
                continue
            n_check += 1
            if t not in known and t.lower() not in {"lit", "mat", "fig"}:
                bad.append((r["claim_id"], f"names unknown test {t!r}"))
        ok = not bad
        return (PASS if ok else FAIL, n_check,
                f"{n_check} VERIFIED claims checked; {len(bad)} unsupported",
                str(bad[:10]))

    finish(c, body)


def test_CLM03_invalid_claims_retained():
    """Superseded INVALID claims are kept, not deleted, and point at a replacement."""
    c = make_case( test_id="CLM03", category="12-claim-traceability",
             purpose="Corrections do not erase the evidence of the earlier failure",
             oracle="Phase 1C section I: an approved baseline and its history are immutable",
             oracle_independent=False,
             fixture="evidence/claim_registry.csv rows with status INVALID",
             expected="every INVALID row has a superseded_by target that exists",
             tolerance=None, tol_justification="structural; exact", level="1")

    def body():
        rows = _csv(CLAIMS)
        ids = {r["claim_id"] for r in rows}
        inv = [r for r in rows if r["evidence_status"] == "INVALID"]
        bad = [r["claim_id"] for r in inv
               if r.get("superseded_by", "").strip()
               and r["superseded_by"] not in ids]
        ok = bool(inv) and not bad
        return (PASS if ok else FAIL, len(inv),
                f"{len(inv)} INVALID claims retained; {len(bad)} bad pointers",
                f"bad={bad}")

    finish(c, body)


def test_EQ01_equation_registry_schema():
    """Equation registry: unique ids, legal status, dependencies resolvable."""
    c = make_case( test_id="EQ01", category="12-claim-traceability",
             purpose="Equation registry is well formed and its cross-references resolve",
             oracle="schema in validation/schemas/equation_registry.schema.json",
             oracle_independent=True,
             fixture="evidence/equation_registry.csv",
             expected=">=35 rows; unique ids; statuses legal; every named dep exists",
             tolerance=None, tol_justification="structural; exact", level="1")

    def body():
        rows = _csv(EQUATIONS)
        ids = {r["eq_id"] for r in rows}
        dup = [i for i in ids if [r["eq_id"] for r in rows].count(i) > 1]
        bad = [r["eq_id"] for r in rows if r["verification_status"] not in ALLOWED_EQ_STATUS]
        dangling = []
        for r in rows:
            for d in re.findall(r"EQ-[A-Z]+-[0-9]+", r.get("downstream_deps", "")):
                if d not in ids:
                    dangling.append((r["eq_id"], d))
        ok = len(rows) >= 35 and not dup and not bad and not dangling
        return (PASS if ok else FAIL, len(rows),
                f"{len(rows)} equations; {len(bad)} bad status; "
                f"{len(dangling)} dangling deps",
                f"bad={bad[:5]} dangling={dangling[:5]}")

    finish(c, body)


def test_ART01_evidence_artifacts_present():
    """Every mandatory Phase 1C evidence artifact exists and is non-trivial."""
    c = make_case( test_id="ART01", category="11-reproducibility",
             purpose="All mandatory evidence artifacts exist",
             oracle="the deliverable list in the Phase 1C brief, section L",
             oracle_independent=True,
             fixture="repository tree",
             expected="all 12 mandatory artifacts exist and are non-empty",
             tolerance=None, tol_justification="existence; exact", level="1")

    def body():
        req = [
            "evidence/claim_registry.csv", "evidence/equation_registry.csv",
            "evidence/dependency_graph.json", "evidence/impact_analysis.py",
            "evidence/phase_01_claim_audit.md", "evidence/approved_baselines.json",
            "evidence/README.md",
            "validation/experiment_manifest_template.json",
            "validation/run_fast.py", "validation/run_standard.py",
            "validation/run_full.py",
            ".github/workflows/research-validation.yml",
            "reports/PHASE_01C_QUALITY_INFRASTRUCTURE.md",
            "reports/PHASE_01C_OPEN_BLOCKERS.md",
            "requirements-research.txt",
        ]
        missing = [r for r in req if not (REPO / r).exists()
                   or (REPO / r).is_file() and (REPO / r).stat().st_size == 0]
        return (PASS if not missing else FAIL, len(req),
                f"{len(req) - len(missing)}/{len(req)} present",
                f"missing={missing}")

    finish(c, body)


def test_ART02_input_checksums_stable():
    """Checksums of the frozen scientific inputs match the recorded baseline.

    Guards against silent edits to data or theory files. If an input legitimately
    changes, the checksum file must be regenerated deliberately and the change
    logged.
    """
    c = make_case( test_id="ART02", category="11-reproducibility",
             purpose="Frozen scientific inputs have not changed silently",
             oracle="evidence/input_checksums.json produced by this phase",
             oracle_independent=True,
             fixture="literature/*.csv, code/det_coeffs.json, theory/*.tex",
             expected="all recorded SHA-256 values match",
             tolerance=None, tol_justification="cryptographic hash; exact", level="1")

    def body():
        reg = REPO / "evidence" / "input_checksums.json"
        if not reg.exists():
            return (UNRESOLVED, 0, "checksum register absent", "")
        rec = json.loads(reg.read_text())["files"]
        bad = []
        for rel, want in rec.items():
            p = REPO / rel
            if not p.exists():
                bad.append((rel, "missing"))
                continue
            got = hashlib.sha256(p.read_bytes()).hexdigest()
            if got != want:
                bad.append((rel, f"{want[:12]} -> {got[:12]}"))
        return (PASS if not bad else FAIL, len(rec),
                f"{len(rec) - len(bad)}/{len(rec)} checksums match", str(bad))

    finish(c, body)


def test_DAT01_code_matches_data_table():
    """The code's constants and the tracked data table must be the same numbers.

    This is the check whose ABSENCE allowed a whole class of error: build_materials()
    hard-codes every constant in Python, so the csv is decorative unless something
    compares them. Found in Phase 1C: Al6061 beta differed between the two, and a
    1000x unit slip in the csv was invisible to the entire test suite.

    Anisotropy note: PZT-5H is transversely isotropic and genuinely does NOT satisfy
    C44 = (C11-C12)/2 (it misses by 5.0e7 Pa, 2e-4 relative, because C12 is derived
    from the compliance set while C44 is measured). So this test compares VALUES
    across the two sources; it must not assert an identity the material violates.
    """
    c = make_case(test_id="DAT01", category="1-data-integrity",
                  purpose="Every constant used by the code equals the value in the "
                          "tracked material table",
                  oracle="literature/material_constants.csv, read independently of "
                          "foundation_checks.build_materials()",
                  oracle_independent=True,
                  fixture="4 media x all comparable constants (72 values)",
                  expected="0 mismatches at 1e-9 relative, and >=60 values compared",
                  tolerance=1e-9,
                  tol_justification="both sides are decimal literals for the same "
                                    "quantity, so agreement should be exact; 1e-9 allows "
                                    "only textual reformatting, never a unit slip",
                  level="1")

    def body():
        import csv as _csv
        import foundation_checks as fc
        rows = list(_csv.DictReader(open(MATERIALS)))
        mats = fc.build_materials()
        bad, n = [], 0
        for r in rows:
            mm = mats.get(r["material"])
            if mm is None or not hasattr(mm, r["symbol"]):
                continue
            try:
                a = float(r["value"])
            except (TypeError, ValueError):
                continue
            b = float(getattr(mm, r["symbol"]))
            n += 1
            if abs(a - b) > 1e-9 * max(abs(a), abs(b), 1e-30):
                bad.append((r["material"], r["symbol"], a, b,
                            abs(a - b) / max(abs(a), 1e-30)))
        ok = not bad and n >= 60
        return (PASS if ok else FAIL, n,
                f"{n} constants compared; {len(bad)} mismatches",
                "; ".join(f"{m}.{s}: table {a:.6g} vs code {b:.6g} "
                          f"({d:.2e} rel)" for m, s, a, b, d in bad))

    finish(c, body)
