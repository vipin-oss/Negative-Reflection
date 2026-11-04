"""
Latex <-> implementation and active-theory consistency.

Task 9 of the Phase 1D brief requires that "corrected LaTeX matches
implementation" and that "no scientific blocker was silently closed". Those are
two different failure modes and both are checked here.

  THY01  the flux formula stated in the theory is the one the code computes
  THY02  the withdrawn tau_q causality floor does not appear as an ACTIVE claim
         in any current theory file (it must survive only in the INVALID
         registry record)
  THY03  the theory files are structurally sound and their labels resolve
  THY04  no blocker was closed without evidence, and G4/G5 are not PASS

Reproduce:
    cd validation && python3 -m pytest tests/test_10_theory_consistency.py -v
"""
from __future__ import annotations

import json
import re

from harness import REPO, PASS, FAIL, make_case, finish

THY = REPO / "theory"
ENERGY = THY / "energy_flux_derivation.tex"
GOVERN = THY / "governing_equations.tex"


# ===========================================================================
def test_THY01_flux_formula_matches_code():
    """The thermal term in the theory must be the bilinear one the code uses.

    A stale theory note is a scientific defect even when the code is right: the
    note is what a reader cites. The specific regression guarded here is the
    Phase-1 error where the note and the code both carried `1/2 Re[q_z]`.
    """
    c = make_case(
        test_id="THY01", category="3-thermodynamics",
        purpose="The cycle-averaged flux stated in the theory note is the bilinear "
                "form the production code evaluates",
        oracle="pattern comparison between theory/energy_flux_derivation.tex and the "
               "return expression of flux_z_mp in code/foundation_checks.py",
        oracle_independent=True,
        fixture="theory/energy_flux_derivation.tex and code/foundation_checks.py",
        expected="the note contains the bilinear q*Theta-conjugate term and the code "
                 "contains the matching mp.conj(Th) expression; neither contains the "
                 "withdrawn linear form as an active statement",
        tolerance=None, tol_justification="exact token presence/absence",
        level="1")

    def body():
        theory = ENERGY.read_text()
        code = (REPO / "code" / "foundation_checks.py").read_text()
        # the bilinear thermal term, written defensively against formatting
        bilinear_theory = bool(re.search(r"q_z\s*\\?Theta\^?\*", theory)
                               or "Re[q_z\\Theta^*]" in theory)
        bilinear_code = "mp.conj(Th)" in code and "qz" in code
        # the withdrawn LINEAR form must not be presented as correct. It may
        # appear only where it is explicitly marked wrong.
        linear_active = False
        for ln in theory.splitlines():
            if "+1/2 Re[q_z]" in ln or "+\\tfrac12\\Re[q_z]" in ln:
                low = ln.lower()
                if not any(k in low for k in ("was", "wrong", "withdrawn", "phase-1",
                                              "phase 1", "linear", "cannot", "supersede")):
                    linear_active = True
        return bilinear_theory, bilinear_code, linear_active

    def judge(res):
        bt, bc, lin = res
        detail = (f"theory states the bilinear thermal term: {bt}\n"
                  f"code evaluates a conjugate thermal term:    {bc}\n"
                  f"withdrawn linear form shown as active:       {lin}")
        if lin:
            return FAIL, "the withdrawn linear thermal term is presented as active", detail
        if not (bt and bc):
            return FAIL, "theory and code do not both carry the bilinear thermal term", detail
        return PASS, "theory and implementation agree on the thermal flux term", detail

    def run():
        status, actual, detail = judge(body())
        return status, 1, actual, detail

    finish(c, run)


# ===========================================================================
def test_THY02_withdrawn_causality_floor_not_active():
    """The withdrawn tau_q causality floor must be absent from ACTIVE theory.

    The review requires the withdrawn statement to be removed from the active
    formulation while its record is retained in the INVALID claim registry. Both
    halves are checked: absence from the active statement, presence in the
    registry. Removing it from the registry as well would destroy the audit
    trail and is a failure here too.
    """
    c = make_case(
        test_id="THY02", category="1-data-integrity",
        purpose="The withdrawn tau_q causality floor is gone from active theory but "
                "still recorded in the INVALID claim registry",
        oracle="absence from an ACTIVE assertion in the theory files, presence in "
               "evidence/claim_registry.csv as INVALID; the two are checked "
               "independently so that neither can be 'fixed' by deleting history",
        oracle_independent=True,
        fixture="theory/governing_equations.tex, theory/energy_flux_derivation.tex, "
                "evidence/claim_registry.csv",
        expected="no ACTIVE causality-floor assertion remains; CLM-H1-02 still exists "
                 "with evidence_status INVALID",
        tolerance=None, tol_justification="presence/absence with an explicit "
                                          "active-vs-withdrawn distinction",
        level="1")

    def body():
        registry_rows = []
        import csv
        with (REPO / "evidence" / "claim_registry.csv").open(newline="") as fh:
            registry_rows = list(csv.DictReader(fh))
        h102 = [r for r in registry_rows if r["claim_id"] == "CLM-H1-02"]
        registry_ok = bool(h102) and h102[0]["evidence_status"] == "INVALID"

        # an ACTIVE assertion looks like a requirement, without a withdrawal mark
        # The test must distinguish an ACTIVE assertion from a RETRACTION that
        # necessarily quotes the withdrawn rule. Keyword matching on one line is
        # too crude: a retraction reads "previously asserted a causality floor",
        # which contains no negation. We therefore require the offending line to
        # also be free of retraction vocabulary AND not sit inside a paragraph
        # that is explicitly marked withdrawn.
        WITHDRAWAL = ("withdraw", "retract", "invalid", "no such", "not a",
                      "must not", "was wrong", "earlier revision", "previously",
                      "unsupported", "superseded", "no longer", "is not supportable",
                      "for reference only", "never")
        active_hits = []
        for f in (GOVERN, ENERGY):
            paras = f.read_text().split("\n\n")
            for para in paras:
                low = para.lower()
                if not ("causal" in low and ("floor" in low or "require" in low)):
                    continue
                # if ANY line in this paragraph carries withdrawal vocabulary, the
                # paragraph is a retraction or a reference note, not an assertion
                if any(k in low for k in WITHDRAWAL):
                    continue
                first = para.strip().splitlines()[0][:110]
                active_hits.append(f"{f.name}: {first}")
        return registry_ok, active_hits

    def judge(res):
        registry_ok, active_hits = res
        detail = (f"CLM-H1-02 still INVALID in the registry: {registry_ok}\n"
                  f"active causality-floor assertions found: {len(active_hits)}")
        if active_hits:
            detail += "\n" + "\n".join(active_hits[:5])
        if not registry_ok:
            return FAIL, ("the withdrawn claim is no longer recorded as INVALID -- "
                          "history has been erased"), detail
        if active_hits:
            return FAIL, f"{len(active_hits)} active causality-floor assertion(s) remain", detail
        return PASS, "withdrawn from active theory, preserved in the INVALID registry", detail

    def run():
        status, actual, detail = judge(body())
        return status, 2, actual, detail

    finish(c, run)


# ===========================================================================
def test_THY03_latex_structurally_sound():
    """Balanced environments, balanced braces, and no dangling cross-references."""
    c = make_case(
        test_id="THY03", category="11-reproducibility",
        purpose="The theory sources are structurally valid and every internal reference "
                "resolves, so a reader can build them",
        oracle="independent structural parse of the .tex sources: environment counts, "
               "brace depth accounting with escapes handled, label/reference set "
               "comparison",
        oracle_independent=True,
        fixture="theory/governing_equations.tex, theory/energy_flux_derivation.tex",
        expected="one document environment each, every begin matched by an end, final "
                 "brace depth 0, and every internal ref resolving to a label",
        tolerance=None, tol_justification="exact structural equality",
        level="2")

    def body():
        problems = []
        for f in (GOVERN, ENERGY):
            s = f.read_text()
            if s.count("\\begin{document}") != 1 or s.count("\\end{document}") != 1:
                problems.append(f"{f.name}: document environment count wrong")
            for env in ("equation", "align", "itemize", "enumerate"):
                nb = len(re.findall(r"\\begin\{" + env + r"\}", s))
                ne = len(re.findall(r"\\end\{" + env + r"\}", s))
                if nb != ne:
                    problems.append(f"{f.name}: {env} {nb} begin vs {ne} end")
            depth, i = 0, 0
            while i < len(s):
                ch = s[i]
                if ch == "\\":
                    i += 2
                    continue
                if ch == "{":
                    depth += 1
                elif ch == "}":
                    depth -= 1
                    if depth < 0:
                        problems.append(f"{f.name}: unmatched closing brace")
                        depth = 0
                i += 1
            if depth != 0:
                problems.append(f"{f.name}: final brace depth {depth}")
            labels = set(re.findall(r"\\label\{([^}]+)\}", s))
            refs = set(re.findall(r"\\(?:eq)?ref\{([^}]+)\}", s))
            missing = refs - labels
            if missing:
                problems.append(f"{f.name}: unresolved refs {sorted(missing)}")
        return problems

    def judge(problems):
        detail = "\n".join(problems) if problems else "no structural problems found"
        if problems:
            return FAIL, f"{len(problems)} structural problem(s)", detail
        return PASS, "2 documents structurally sound, all refs resolve", detail

    def run():
        status, actual, detail = judge(body())
        return status, 2, actual, detail

    finish(c, run)


# ===========================================================================
def test_THY04_no_blocker_silently_closed():
    """No blocker may be CLOSED without recorded closure evidence, and G4/G5
    may not be reported PASS.

    This is the task-9 requirement stated as a test. It is deliberately
    independent of the blocker report's own prose: it reads the JSON and demands
    that a CLOSED blocker carries a non-empty closure_evidence object, and that
    the two energy gates are still failing.
    """
    c = make_case(
        test_id="THY04", category="12-claim-traceability",
        purpose="Every CLOSED blocker carries explicit closure evidence, and the energy "
                "gates G4 and G5 are not silently promoted",
        oracle="reports/PHASE_01C_OPEN_BLOCKERS.json read directly; the G2/G3/G4/G5 "
               "gate states in evidence/claim_registry.csv cross-checked against it",
        oracle_independent=True,
        fixture="blocker register and claim registry",
        expected="every CLOSED blocker has non-empty closure_evidence; BLK-G4 and "
                 "BLK-G5 are not marked PASS or VERIFIED anywhere",
        tolerance=None, tol_justification="exact structural requirement",
        level="1")

    def body():
        blk = json.loads((REPO / "reports" / "PHASE_01C_OPEN_BLOCKERS.json").read_text())
        problems, n_closed = [], 0
        for b in blk["blockers"]:
            st = str(b.get("status", "OPEN")).upper()
            if "CLOSED" in st and "PART" not in st:
                n_closed += 1
                ev = b.get("closure_evidence")
                if not ev or not isinstance(ev, dict) or not ev:
                    problems.append(f"{b['id']}: CLOSED with no closure_evidence")
        # the energy gates must remain unpassed
        gates = {b["id"]: b for b in blk["blockers"]}
        for gid in ("BLK-G4", "BLK-G5"):
            g = gates.get(gid)
            if g is None:
                problems.append(f"{gid} missing from the register")
                continue
            if str(g.get("status", "")).upper() in ("PASS", "CLOSED", "VERIFIED"):
                problems.append(f"{gid} has been marked {g.get('status')!r}; "
                                f"a special elastic limit may not promote it")
        # and G0/G1, which are also P0
        for gid in ("BLK-G0", "BLK-G1"):
            g = gates.get(gid)
            if g and str(g.get("status", "")).upper() in ("PASS", "CLOSED"):
                problems.append(f"{gid} closed without a scattering solution")
        return problems, n_closed, len(blk["blockers"])

    def judge(res):
        problems, n_closed, total = res
        detail = (f"{total} blockers registered, {n_closed} CLOSED\n"
                  + ("\n".join(problems) if problems else
                     "every closed blocker carries evidence; G0/G1/G4/G5 remain unpassed"))
        if problems:
            return FAIL, f"{len(problems)} blocker-integrity problem(s)", detail
        return PASS, (f"{n_closed} closed blocker(s) all evidenced; "
                      f"G0/G1/G4/G5 correctly unpassed"), detail

    def run():
        status, actual, detail = judge(body())
        return status, 1, actual, detail

    finish(c, run)
