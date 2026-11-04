"""
Category 12 (continued) -- manuscript traceability, section J.

No manuscript content is written in Phase 1C. What is built here is the mechanism
that will refuse to let unsupported content be written later.
"""
from __future__ import annotations

import json

from harness import REPO, PASS, FAIL, UNRESOLVED, make_case, finish

MAP = REPO / "evidence" / "manuscript_evidence_map.json"
BASELINES = REPO / "evidence" / "approved_baselines.json"


def _map():
    return json.loads(MAP.read_text())


def _baselines():
    return json.loads(BASELINES.read_text())


def test_BIND05_stale_baseline_blocks_claims():
    """A claim cannot stay VERIFIED once its baseline is STALE or SUPERSEDED."""
    c = make_case(test_id="BIND05", category="12-claim-traceability",
                  purpose="Stale evidence downgrades claims and blocks the manuscript "
                          "entries that depend on it",
                  oracle="evidence/approved_baselines.json status field versus "
                         "evidence/manuscript_evidence_map.json entries",
                  oracle_independent=True,
                  fixture="baselines + manuscript map",
                  expected="every entry whose baseline is not ACTIVE is BLOCKED or STALE, "
                           "and no claim sourced only from a stale baseline is VERIFIED",
                  tolerance=None, tol_justification="structural consistency", level="1")

    def body():
        b = _baselines()
        m = _map()
        stale = {x["baseline_id"] for x in b["baselines"]
                 if x["status"] != "ACTIVE"}
        # An entry may be EVIDENCE_READY only if the SPECIFIC baseline it rests on
        # is ACTIVE. The historical rejected baselines are retained on purpose and
        # must not, by their mere existence, block current work - that would make
        # the honest retention of a failure punishing.
        status_of = {x["baseline_id"]: x["status"] for x in b["baselines"]}
        bad = []
        for e in m["entries"]:
            ref = e.get("baseline_ref")
            if e["status"] == "EVIDENCE_READY":
                if ref is None:
                    bad.append((e["id"], "EVIDENCE_READY with no baseline_ref"))
                elif status_of.get(ref) != "ACTIVE":
                    bad.append((e["id"], f"rests on {ref} which is "
                                         f"{status_of.get(ref)}"))
        ok = not bad and bool(b["baselines"])
        return (PASS if ok else FAIL, len(m["entries"]) + len(b["baselines"]),
                f"{len(b['baselines'])} baselines ({len(stale)} not active, retained "
                f"on purpose); "
                f"{len(m['entries'])} manuscript entries; {len(bad)} violations",
                str(bad))

    finish(c, body)


def test_BIND06_manuscript_entries_resolvable():
    """Every evidence id a manuscript entry cites must exist in a registry."""
    c = make_case(test_id="BIND06", category="12-claim-traceability",
                  purpose="A manuscript item cannot cite evidence that does not exist",
                  oracle="the equation registry, the claim registry and the dependency "
                         "graph, taken as the universe of legal evidence ids",
                  oracle_independent=True,
                  fixture="evidence/manuscript_evidence_map.json",
                  expected="0 unresolvable ids; 0 entries citing an INVALID equation "
                           "without being marked BLOCKED",
                  tolerance=None, tol_justification="set membership", level="1")

    def body():
        import csv
        eq = {r["eq_id"]: r["verification_status"]
              for r in csv.DictReader(open(REPO / "evidence" / "equation_registry.csv"))}
        cl = {r["claim_id"]: r["evidence_status"]
              for r in csv.DictReader(open(REPO / "evidence" / "claim_registry.csv"))}
        m = _map()
        bad, n = [], 0
        for e in m["entries"]:
            for x in e.get("evidence", []):
                n += 1
                if x not in eq:
                    bad.append((e["id"], f"unknown equation {x}"))
                elif eq[x] == "INVALID" and e["status"] not in ("BLOCKED", "FORBIDDEN",
                                                                "STALE"):
                    bad.append((e["id"], f"cites INVALID {x} without being BLOCKED"))
            for x in e.get("claims", []):
                n += 1
                if x not in cl:
                    bad.append((e["id"], f"unknown claim {x}"))
                elif cl[x] == "INVALID" and e["status"] not in ("BLOCKED", "FORBIDDEN",
                                                                "STALE"):
                    bad.append((e["id"], f"cites INVALID claim {x} unblocked"))
        ok = not bad
        return (PASS if ok else FAIL, n, f"{n} evidence references; {len(bad)} problems",
                str(bad[:6]))

    finish(c, body)


def test_BIND07_baseline_is_immutable():
    """Approved baselines are append-only: no id may change its recorded commit."""
    c = make_case(test_id="BIND07", category="12-claim-traceability",
                  purpose="An approved baseline is never silently overwritten",
                  oracle="the permanent supersession chain in approved_baselines.json: "
                         "each baseline's commit SHA must be unique and each "
                         "superseded_by must point to a later entry",
                  oracle_independent=True,
                  fixture="evidence/approved_baselines.json",
                  expected="unique SHAs; valid supersession pointers; every superseded "
                           "entry retains its original outputs",
                  tolerance=None, tol_justification="structural", level="1")

    def body():
        b = _baselines()
        bl = b["baselines"]
        ids = [x["baseline_id"] for x in bl]
        shas = [x["code_commit"] for x in bl]
        n = len(bl)
        dup_id = len(set(ids)) != n
        dup_sha = len(set(shas)) != n
        bad = []
        for x in bl:
            if x["status"] in ("SUPERSEDED", "STALE"):
                if not x.get("supersession_reason"):
                    bad.append((x["baseline_id"], "no supersession reason recorded"))
                if x.get("superseded_by") and x["superseded_by"] not in ids:
                    bad.append((x["baseline_id"], "dangling superseded_by"))
            if x["status"] == "ACTIVE" and x.get("superseded_by"):
                bad.append((x["baseline_id"], "ACTIVE but points at a successor"))
        ok = not dup_id and not dup_sha and not bad
        return (PASS if ok else FAIL, n,
                f"{n} baselines; duplicate ids={dup_id} duplicate SHAs={dup_sha}; "
                f"{len(bad)} chain problems", str(bad))

    finish(c, body)


def test_BIND08_figures_must_be_regenerable():
    """A figure with no regeneration path may not be cited as evidence."""
    c = make_case(test_id="BIND08", category="11-reproducibility",
                  purpose="Every figure cited by a manuscript entry has a recorded "
                          "regeneration command, present script and present raw data",
                  oracle="evidence/manuscript_evidence_map.json regeneration blocks "
                         "versus the filesystem",
                  oracle_independent=True,
                  fixture="reports/figures/, validation/figures/",
                  expected="for each cited figure: script_present and raw_data_present "
                           "are true, or the entry is not EVIDENCE_READY",
                  tolerance=None, tol_justification="existence checks", level="2")

    def body():
        m = _map()
        bad, n = [], 0
        for e in m["entries"]:
            if e["kind"] != "figure":
                continue
            n += 1
            r = e.get("regeneration", {})
            if e["status"] == "EVIDENCE_READY" and not (r.get("script_present")
                                                        and r.get("raw_data_present")):
                bad.append((e["id"], "EVIDENCE_READY without a regeneration path"))
            if r.get("figure_file") and not (REPO / r["figure_file"]).exists():
                bad.append((e["id"], f"named figure missing {r['figure_file']}"))
        ok = not bad
        return (PASS if ok else FAIL, n,
                f"{n} figure entries; {len(bad)} without a regeneration path", str(bad))

    finish(c, body)
