"""
Category 12 -- research-claim traceability, enforced against the CODE.

test_01 checks that a claim names a test that exists. That is not enough: U-11
happened because the test kept its name while the calculation it was supposed to
perform was deleted. This module closes that hole by requiring, for every
VERIFIED claim, that the named function still contains the tokens that carry the
asserted quantity (evidence/claim_evidence_map.json).
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

from harness import REPO, PASS, FAIL, UNRESOLVED, make_case, finish

MAP = REPO / "evidence" / "claim_evidence_map.json"
CLAIMS = REPO / "evidence" / "claim_registry.csv"


def _claims():
    with open(CLAIMS, newline="") as fh:
        return [r for r in csv.DictReader(fh) if r.get("claim_id")]


def _fn_source(fn: str, rel: str = "code/foundation_checks.py") -> str | None:
    p = REPO / rel
    if not p.exists():
        return None
    src = p.read_text()
    m = re.search(rf"^def {re.escape(fn)}\(.*?(?=^def |\Z)", src, re.S | re.M)
    return m.group(0) if m else None


# ===========================================================================
def test_BIND01_every_verified_claim_is_bound():
    """Every VERIFIED claim names a test that still executes its assertion."""
    c = make_case(test_id="BIND01", category="12-claim-traceability",
                  purpose="No claim may stay VERIFIED when the code that asserts it has "
                          "been replaced, renamed or gutted",
                  oracle="evidence/claim_evidence_map.json: for each VERIFIED claim, the "
                         "named function must exist and contain every required token",
                  oracle_independent=True,
                  fixture="44 registry rows cross-referenced against the source of "
                          "code/foundation_checks.py and validation/tests/",
                  expected="0 problems; every VERIFIED claim resolves to BOUND",
                  tolerance=None, tol_justification="exact token presence", level="1")

    def body():
        bmap = json.loads(MAP.read_text())
        bound = {b["claim_id"]: b for b in bmap["bindings"]}
        rows = _claims()
        problems, n_verified, n_bound_ok = [], 0, 0
        for r in rows:
            cid, st = r["claim_id"], r["evidence_status"]
            if not st.startswith("VERIFIED"):
                continue
            n_verified += 1
            b = bound.get(cid)
            if b is None:
                # a VERIFIED claim with no binding is only legal if the *audit*
                # deliberately records it as VERIFIED_WEAK with a stated defect
                problems.append((cid, "VERIFIED but NOT BOUND in claim_evidence_map"))
                continue
            # PHASE 1D GAP CLOSED. Until Phase 1D this check trusted the map
            # entirely and never compared the map's `check` field against the
            # registry's own `test_id`. Mutation M18 exploited exactly that:
            # rewriting CLM-D2-01's test_id to a test that does not compute the
            # balance was INVISIBLE, because the map still pointed at the real
            # function. A registry row that names a different test from its
            # binding is now a hard failure.
            declared = (r.get("test_id") or "").strip()
            bound_check = (b.get("check") or "").strip()
            if declared and bound_check and declared != bound_check:
                problems.append((cid,
                                 f"TEST_ID_DESYNC registry says {declared!r} but the "
                                 f"binding says {bound_check!r}"))
                continue
            rel = b.get("file", "code/foundation_checks.py")
            if b.get("binding_kind") == "data_file":
                # a claim verified by the CONTENT of a data file rather than by code
                src = (REPO / rel).read_text() if (REPO / rel).exists() else None
            else:
                src = _fn_source(b["fn"], rel)
            if src is None:
                problems.append((cid, f"FUNCTION/COLUMN MISSING {b.get('fn')} in {rel}"))
                continue
            miss = [t for t in b.get("required_tokens", []) if t not in src]
            if miss:
                problems.append((cid, f"TOKEN_MISSING {miss} in {b['fn']}"))
                continue
            n_bound_ok += 1
        ok = not problems
        return (PASS if ok else FAIL, n_verified,
                f"{n_verified} VERIFIED claims; {n_bound_ok} BOUND; "
                f"{len(problems)} problems", str(problems[:8]))

    finish(c, body)


def test_BIND02_known_defects_still_detectable():
    """The two documented U-11 defects must still be detectable from the source.

    CLM-C1-04 asserts a symbolic null space. If someone quietly re-adds a
    nullspace() call the claim would look bound; if someone removes the *claim*,
    the record of the failure disappears. This test pins the current, defective
    state so that any change to it is deliberate.
    """
    c = make_case(test_id="BIND02", category="12-claim-traceability",
                  purpose="A documented evidence defect cannot be silently repaired or "
                          "silently forgotten",
                  oracle="the recorded defect state in claim_evidence_map.json",
                  oracle_independent=False,
                  fixture="check_C1_structure, check_I2_interface_rank, check_F1_passivity",
                  expected="C1 still lacks a symbolic nullspace; I2 still hardcodes PASS; "
                           "F1 still samples a PSD matrix",
                  tolerance=None, tol_justification="exact source inspection", level="1")

    def body():
        c1 = _fn_source("check_C1_structure") or ""
        i2 = _fn_source("check_I2_interface_rank") or ""
        f1 = _fn_source("check_F1_passivity") or ""
        state = dict(
            c1_has_nullspace="nullspace" in c1,
            i2_hardcoded_pass='return "PASS"' in i2,
            f1_samples_psd=("A @ A.T" in f1 or "A@A.T" in f1),
        )
        # these are the DEFECTS as audited; if any flips, the registry must be updated
        expected = dict(c1_has_nullspace=False, i2_hardcoded_pass=True,
                        f1_samples_psd=True)
        changed = {k: (state[k], expected[k]) for k in state if state[k] != expected[k]}
        ok = not changed
        return (PASS if ok else FAIL, 3,
                f"defect state {state}; {'as recorded' if ok else 'CHANGED: ' + str(changed)}",
                "if this fails because a defect was FIXED, update "
                "claim_evidence_map.json and the claim registry in the same commit")

    finish(c, body)


def test_BIND03_no_orphan_bindings():
    """The binding file and the registry must not disagree about which claims exist."""
    c = make_case(test_id="BIND03", category="12-claim-traceability",
                  purpose="Binding entries and registry rows are in one-to-one "
                          "correspondence; neither can drift",
                  oracle="set equality between bound claim ids and registry ids, minus "
                         "the explicitly declared negative bindings",
                  oracle_independent=True,
                  fixture="claim_evidence_map.json + claim_registry.csv",
                  expected="0 unknown ids; 0 VERIFIED claims missing a binding; every "
                           "negative binding names a registry row",
                  tolerance=None, tol_justification="set equality", level="1")

    def body():
        bmap = json.loads(MAP.read_text())
        rows = _claims()
        ids = {r["claim_id"] for r in rows}
        status = {r["claim_id"]: r["evidence_status"] for r in rows}
        unknown = [b["claim_id"] for b in bmap["bindings"] if b["claim_id"] not in ids]
        neg_unknown = [b["claim_id"] for b in bmap.get("negative_bindings", [])
                       if b["claim_id"] not in ids]
        neg_ids = {b["claim_id"] for b in bmap.get("negative_bindings", [])}
        bound_ids = {b["claim_id"] for b in bmap["bindings"]}
        # a non-VERIFIED claim must not sit in the positive binding list
        wrong = [cid for cid in bound_ids
                 if cid in status and not status[cid].startswith("VERIFIED")]
        ok = not unknown and not neg_unknown and not wrong
        return (PASS if ok else FAIL, len(bound_ids) + len(neg_ids),
                f"{len(bound_ids)} positive, {len(neg_ids)} negative bindings; "
                f"{len(unknown) + len(neg_unknown)} unknown ids; {len(wrong)} "
                f"non-VERIFIED claims bound as if verified",
                f"unknown={unknown + neg_unknown} wrong={wrong}")

    finish(c, body)


def test_BIND04_status_tallies_consistent():
    """The registry tally must match what the audit report claims."""
    c = make_case(test_id="BIND04", category="12-claim-traceability",
                  purpose="The published status tally equals the registry contents",
                  oracle="recount from the csv, compared against the numbers stated in "
                         "reports/PHASE_01C_VALIDATION_SUMMARY.json",
                  oracle_independent=False,
                  fixture="evidence/claim_registry.csv",
                  expected="tallies agree exactly",
                  tolerance=None, tol_justification="integer counts", level="1")

    def body():
        cnt = Counter(r["evidence_status"] for r in _claims())
        summ = REPO / "reports" / "PHASE_01C_VALIDATION_SUMMARY.json"
        if not summ.exists():
            return (UNRESOLVED, len(cnt), "summary json not yet written", str(dict(cnt)))
        j = json.loads(summ.read_text())
        claimed = j.get("claims", {}).get("status_counts")
        if claimed is None:
            return (UNRESOLVED, len(cnt), "summary has no claim tallies", "")
        agree = {k: (cnt.get(k, 0), claimed.get(k)) for k in set(cnt) | set(claimed)
                 if cnt.get(k, 0) != claimed.get(k)}
        ok = not agree
        return (PASS if ok else FAIL, len(cnt),
                f"registry {dict(cnt)}", f"disagreements {agree}")

    finish(c, body)
