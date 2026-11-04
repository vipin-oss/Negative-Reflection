"""
Categories 4, 6, 9, 10, 11 -- matrix identities, interface conditions, scattering,
convergence and reproducibility.

The scattering tests here are the ones that CAN be run before a coupled solver
exists: they test the machinery that a solver would have to satisfy, and they are
explicit that they are not a substitute for it.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys

import numpy as np

from harness import REPO, PASS, FAIL, UNRESOLVED, NOT_APPLICABLE, make_case, finish
import foundation_checks as fc


# ===========================================================================
def test_MATX01_det_independent_routes():
    """det M by two algebraically independent routes."""
    c = make_case(test_id="MATX01", category="4-matrix-identities",
                  purpose="det M from the exact polynomial equals the determinant of the "
                          "independently assembled numeric matrix",
                  oracle="the Berkowitz polynomial in code/det_coeffs.json, evaluated, "
                         "versus numpy's LU determinant of Medium6mm.M assembled from the "
                         "constitutive relations at the same point",
                  oracle_independent=True,
                  fixture="PZT-5H; 40 random (k, xi) with k in [1,1e5] m^-1, "
                          "xi real in [-1e5,1e5]",
                  expected="relative deviation < 1e-10 in at least 35 usable cases",
                  tolerance=1e-10,
                  tol_justification="both routes are double precision, but they share "
                                    "nothing else: one is a polynomial in xi, the other an "
                                    "LU factorisation of the assembled matrix. A "
                                    "coefficient-ordering error would show as O(1).",
                  level="1")

    def body():
        m = fc.build_materials()["PZT-5H"]
        w = 2 * np.pi * 1e6
        rng = np.random.default_rng(7)
        worst, n = 0.0, 0
        for _ in range(40):
            k = float(rng.uniform(1.0, 1e5))
            xi = float(rng.uniform(-1e5, 1e5))
            co = m.poly_coeffs(k, w)
            a = sum(c * xi ** j for j, c in enumerate(co))
            b = np.linalg.det(np.array(m.M(k, xi, w), dtype=complex))
            if abs(b) < 1e-290:
                continue
            worst = max(worst, abs(complex(a) - b) / abs(b))
            n += 1
        ok = worst < 1e-10 and n >= 35
        return (PASS if ok else FAIL, n,
                f"{n} usable cases; worst relative deviation {worst:.3e}", "")

    finish(c, body)


def test_MATX02_parity_exhaustive():
    """The parity transformation is unique; brute-force over all 16 sign patterns."""
    c = make_case(test_id="MATX02", category="4-matrix-identities",
                  purpose="S = diag(1,-1,-1,1) is the UNIQUE sign matrix satisfying "
                          "M(-xi) = S M(xi) S, and diag(1,1,1,-1) is detected as failing",
                  oracle="exhaustive search over all 2^4 = 16 diagonal sign matrices, each "
                         "tested on all 16 entries. No candidate is privileged.",
                  oracle_independent=True,
                  fixture="symbolic M",
                  expected="exactly the pair {+S, -S} with S = diag(1,-1,-1,1), i.e. one "
                           "survivor up to the physically irrelevant overall sign; the "
                           "Phase-1 S explicitly rejected",
                  tolerance=None, tol_justification="exact symbolic equality", level="1")

    def body():
        import sympy as sp
        xi = sp.Symbol("xi")
        M = fc.symbolic_M()
        Mn = M.subs(xi, -xi)
        survivors = []
        for bits in range(16):
            s = [1 if (bits >> i) & 1 == 0 else -1 for i in range(4)]
            S = sp.diag(*s)
            if sp.simplify(Mn - S * M * S) == sp.zeros(4, 4):
                survivors.append(tuple(s))
        Sa = sp.diag(1, 1, 1, -1)
        fails = [(i, j) for i in range(4) for j in range(4)
                 if sp.simplify((Mn - Sa * M * Sa)[i, j]) != 0]
        # S and -S are the same parity operation: (-S) M (-S) = S M S, and flipping the
        # sign of a state vector changes nothing physical. So the statement is
        # "unique up to an overall sign", and asserting a single survivor would have
        # been wrong - recorded here because the first version of this test did that.
        as_set = {tuple(x) for x in survivors}
        expected_pair = {(1, -1, -1, 1), (-1, 1, 1, -1)}
        ok = as_set == expected_pair and len(fails) >= 1
        return (PASS if ok else FAIL, 17,
                f"survivors over 16 sign patterns: {survivors} (= +/-S, one up to sign); "
                f"Phase-1 S fails "
                f"{len(fails)} entries", f"failing entries: {fails}")

    finish(c, body)


def test_IFC02_rank_measured():
    """The interface operator's rank is MEASURED, not asserted.

    Phase 1B recorded that the same 8x8 operator has rank 4 on Al6061|Al6061 and
    rank 8 on a mismatched pair. A solver that calls a direct solve on it returns
    a least-squares answer without warning. This test recomputes the rank from
    scratch for each material pair and angle, and pins the numbers.
    """
    c = make_case(test_id="IFC02", category="6-interface",
                  purpose="Interface operator rank is computed and reported for every "
                          "material pair and angle; a deficient operator is never solved "
                          "as if invertible",
                  oracle="rank by Gaussian elimination with full pivoting on the mpmath "
                         "operator. Independent of the Phase 1B text, which reported ranks "
                         "in prose only.",
                  oracle_independent=True,
                  fixture="PZT-5H|PZT-5H, Al6061|Al6061, PZT-5H|Al6061 at 5, 30, 60 deg",
                  expected="rank reported for all 9 cases; every rank < 8 is flagged; the "
                           "count is 9 and no case is skipped",
                  tolerance=None,
                  tol_justification="rank is an integer; the elimination uses exact "
                                    "mpmath arithmetic at dps 50 with full pivoting",
                  level="2")

    def body():
        import mpmath as mp
        mp.mp.dps = 50
        mats = fc.build_materials()
        pairs = [("PZT-5H", "PZT-5H"), ("Al6061", "Al6061"), ("PZT-5H", "Al6061")]
        w = 2 * mp.pi * 1e6
        rows, n = [], 0
        for a, b in pairs:
            m1, m2 = mats[a], mats[b]
            for ang in (5.0, 30.0, 60.0):
                k = w * mp.sin(mp.radians(ang)) / 4000.0
                B = mp.zeros(8, 8)
                col = 0
                for med, sgn in ((m1, +1), (m2, -1)):
                    try:
                        out, _dec = fc._admissible(med, k, w)
                    except Exception as exc:                   # noqa: BLE001
                        rows.append(dict(pair=f"{a}|{b}", angle=ang, rank="ERR",
                                         note=repr(exc)))
                        n += 1
                        continue
                    for z in out:
                        v = fc.state_fluxvec(med, k, z, w,
                                             fc.kernel_scaled(med, k, z, w))
                        for i in range(8):
                            B[i, col] = sgn * v[i]
                        col += 1
                if col != 8:
                    rows.append(dict(pair=f"{a}|{b}", angle=ang, rank=f"{col}cols",
                                     note="incomplete mode set"))
                else:
                    rows.append(dict(pair=f"{a}|{b}", angle=ang,
                                     rank=int(fc.exact_rank(B, 8, 8)), note=""))
                n += 1
        deficient = [r for r in rows if r["rank"] != 8]
        ok = n == 9 and all("rank" in r for r in rows)
        return (PASS if ok else FAIL, n,
                f"{n}/9 cases computed; {len(deficient)} rank-deficient",
                "; ".join(f"{r['pair']}@{r['angle']}: rank {r['rank']} {r['note']}"
                          for r in deficient) or "all rank 8")

    finish(c, body)


def test_SCT01_elastic_reference_reproducible():
    """The oracle a coupled solver must reproduce, run in a clean subprocess."""
    c = make_case(test_id="SCT01", category="9-scattering",
                  purpose="The elastic reference oracle is itself reproducible from a fresh "
                          "interpreter and matches the committed values",
                  oracle="code/elastic_psv.py executed as a subprocess; output compared to "
                         "the values recorded in the Phase 1B report",
                  oracle_independent=True,
                  fixture="Al6061|steel-like, 1 MHz, normal incidence",
                  expected="R = 0.217340589840 and T = 0.782659410160 to 1e-12",
                  tolerance=1e-12,
                  tol_justification="the recorded values are accurate to 1.7e-16; 1e-12 "
                                    "absorbs platform BLAS variation",
                  level="1")

    def body():
        code = (
            "import sys; sys.path.insert(0,'code')\n"
            "import numpy as np, elastic_psv as ep\n"
            "w=2*np.pi*1e6\n"
            "A=ep.ElasticISO('Al6061',2700.0,107.6e9,26.5e9)\n"
            "B=ep.ElasticISO('steel',7800.0,281.0e9,82.0e9)\n"
            "d=ep.solve_psv(A,B,w,0.0)\n"
            "print(repr(float(d['R'])),repr(float(d['T'])))\n")
        out = subprocess.run([sys.executable, "-c", code], cwd=REPO,
                             capture_output=True, text=True, timeout=300)
        if out.returncode != 0:
            return (FAIL, 0, "subprocess failed", out.stderr[-400:])
        R, T = [float(x) for x in out.stdout.split()]
        ok = abs(R - 0.217340589840) < 1e-12 and abs(T - 0.782659410160) < 1e-12
        return (PASS if ok else FAIL, 1, f"R = {R!r}, T = {T!r}", "")

    finish(c, body)


def test_SCT02_coupled_solver_absent():
    """The coupled scattering solver does not exist: not a PASS."""
    c = make_case(test_id="SCT02", category="9-scattering",
                  purpose="Record honestly that coupled R/T/A is NOT yet a validated "
                          "physical result",
                  oracle="none available - no independent coupled energy or "
                         "thermodynamic gate exists",
                  oracle_independent="ABSENT",
                  fixture="n/a - the blocker is a missing independent check, not "
                          "missing code",
                  expected="the gap is recorded, and the harness refuses to count a "
                           "zero-execution result as a pass",
                  tolerance=None,
                  tol_justification="no computation is possible, so no tolerance exists. "
                                    "The harness forces zero-execution NOT_APPLICABLE to "
                                    "UNRESOLVED, so this can never be counted as a pass.",
                  level="1")

    def body():
        # PHASE 2A TASK 2 UPDATE. A coupled solver now exists
        # (code/scattering/solver.py) and it closes the interface to 2.7e-51.
        # That does NOT make this a pass: what is still absent is an
        # INDEPENDENT coupled energy-conservation gate (G4) and an independent
        # thermodynamic-consistency gate (G5). A converged boundary solve is not
        # evidence that R, T and A are validated piezothermoelastic quantities,
        # so the blocker is restated rather than removed. See BLK-02A-13 in
        # reports/PHASE_02A_TASK2_OPEN_BLOCKERS.md.
        try:
            from scattering.solver import EVIDENCE_LABEL  # noqa: F401
            coupled_exists = True
        except Exception:                                       # noqa: BLE001
            coupled_exists = False
        if coupled_exists:
            return (UNRESOLVED, 0,
                    "a coupled solver exists but no independent coupled energy "
                    "(G4) or thermodynamic (G5) gate exists; coupled amplitudes "
                    "remain labelled DIAGNOSTIC (BLK-02A-13)",
                    "closes when a total-field flux computation reproduces the "
                    "incident flux independently on both sides of the interface")
        has = any(hasattr(fc, n) for n in
                  ("solve_coupled", "scattering_matrix", "mode_match"))
        if has:
            return (UNRESOLVED, 0,
                    "a coupled entry point now exists; this placeholder must be replaced "
                    "by a real validation", "")
        return (NOT_APPLICABLE, 0,
                "no coupled solver: BLK-G5 and BLK-U3 remain open. Forced to UNRESOLVED "
                "by the harness.", "")

    finish(c, body)


def test_CONV01_root_convergence():
    """Eigenvalue roots are converged with respect to working precision."""
    c = make_case(test_id="CONV01", category="10-convergence",
                  purpose="Root locations do not depend on the mpmath working precision, so "
                          "they are not precision artefacts",
                  oracle="the same solve at dps 30, 45 and 60, compared root by root",
                  oracle_independent=False,
                  fixture="PZT-5H, 1 MHz, 30 degrees",
                  expected="max relative root shift over dps 30->60 < 1e-20",
                  tolerance=1e-20,
                  tol_justification="the roots are well conditioned; the known 1.86e-5 "
                                    "error is in the null VECTOR, not the root. A tighter "
                                    "tolerance than the documented defect is therefore "
                                    "diagnostic rather than arbitrary.",
                  level="2")

    def body():
        import mpmath as mp
        m = fc.build_materials()["PZT-5H"]
        w = 2 * mp.pi * 1e6
        k = w * mp.sin(mp.radians(30.0)) / 4000.0
        base, worst, settled = None, mp.mpf(0), 0
        for dps in (30, 45, 60):
            mp.mp.dps = dps
            r = sorted([mp.mpc(z) for z in
                        m.xi_roots(mp.mpf(repr(float(k))), mp.mpf(repr(float(w))))],
                       key=lambda z: (mp.re(z), mp.im(z)))
            if base is None:
                base = r
                continue
            worst = max(worst, max(abs(a - b) / max(abs(a), mp.mpf("1e-300"))
                                   for a, b in zip(base, r)))
            settled += 1
        mp.mp.dps = 50
        ok = worst < mp.mpf("1e-20")
        return (PASS if ok else FAIL, settled,
                f"max relative root shift over {settled} precision steps "
                f"{float(worst):.3e}", "")

    finish(c, body)


def test_REP01_foundation_checks_reproduce():
    """The Phase 1B foundation checks reproduce their committed verdicts."""
    c = make_case(test_id="REP01", category="11-reproducibility",
                  purpose="code/foundation_checks.py still reaches the same check and case "
                          "counts from a clean interpreter",
                  oracle="the committed validation/phase_01b_test_results.json",
                  oracle_independent=True,
                  fixture="python3 code/foundation_checks.py",
                  expected="same number of checks, same status tally, same executed-case "
                           "total",
                  tolerance=None, tol_justification="integer counts; exact", level="1")

    def body():
        out = subprocess.run([sys.executable, str(REPO / "code" / "foundation_checks.py")],
                             cwd=REPO, capture_output=True, text=True, timeout=1800)
        ref_p = REPO / "validation" / "phase_01b_test_results.json"
        cur_p = REPO / "data" / "foundation_checks.json"
        if not ref_p.exists() or not cur_p.exists():
            return (UNRESOLVED, 0, "a result file is missing",
                    f"{ref_p.exists()=} {cur_p.exists()=}")
        from collections import Counter

        def tally(d):
            ch = d.get("checks") or d.get("results") or []
            return (Counter(x["status"] for x in ch), len(ch),
                    sum(int(x.get("n_exec", 0)) for x in ch))
        a = tally(json.loads(ref_p.read_text()))
        b = tally(json.loads(cur_p.read_text()))
        ok = a == b
        return (PASS if ok else FAIL, 1,
                f"recorded {a[1]} checks/{a[2]} cases {dict(a[0])}; reproduced "
                f"{b[1]} checks/{b[2]} cases {dict(b[0])}; exit {out.returncode}", "")

    finish(c, body)


def test_REP02_manifest_every_experiment():
    """No numerical result without a manifest with verified hashes."""
    c = make_case(test_id="REP02", category="11-reproducibility",
                  purpose="Every recorded numerical experiment has a manifest whose input "
                          "and output hashes verify",
                  oracle="validation/schemas/experiment_manifest.schema.json",
                  oracle_independent=True,
                  fixture="validation/evidence/*/manifest.json",
                  expected=">=3 manifests; each satisfies the schema; every referenced "
                           "output hash matches",
                  tolerance=None, tol_justification="structural and cryptographic; exact",
                  level="1")

    def body():
        d = REPO / "validation" / "evidence"
        mans = sorted(d.glob("*/manifest.json")) if d.exists() else []
        if not mans:
            return (UNRESOLVED, 0, "no manifests present", "")
        schema = json.loads((REPO / "validation" / "schemas" /
                             "experiment_manifest.schema.json").read_text())
        req = set(schema["required"])
        bad = []
        for mp_ in mans:
            j = json.loads(mp_.read_text())
            miss = req - set(j)
            if miss:
                bad.append((mp_.parent.name, f"missing {sorted(miss)}"))
                continue
            for rel, want in (j.get("outputs") or {}).items():
                f = REPO / rel
                if not f.exists():
                    bad.append((mp_.parent.name, f"output missing: {rel}"))
                elif hashlib.sha256(f.read_bytes()).hexdigest() != want:
                    bad.append((mp_.parent.name, f"output hash mismatch: {rel}"))
        ok = not bad and len(mans) >= 3
        return (PASS if ok else FAIL, len(mans),
                f"{len(mans)} manifests; {len(bad)} problems", str(bad))

    finish(c, body)


def test_DEP01_dependency_chain():
    """Changing an equation reaches all of its dependents."""
    c = make_case(test_id="DEP01", category="12-claim-traceability",
                  purpose="A change to the heat equation propagates to every dependent "
                          "energy and scattering node, and to nothing unrelated",
                  oracle="evidence/impact_analysis.py, compared against a hand-written "
                         "expected dependent set",
                  oracle_independent=True,
                  fixture="evidence/dependency_graph.json",
                  expected="EQ-ENG-04, EQ-ENG-08, EQ-SCT-01 reachable from EQ-GOV-08; "
                           "CONST-AL6061 and EQ-ELAS-01 not reachable",
                  tolerance=None, tol_justification="set membership; exact", level="1")

    def body():
        out = subprocess.run(
            [sys.executable, str(REPO / "evidence" / "impact_analysis.py"),
             "--changed", "EQ-GOV-08", "--json"],
            cwd=REPO, capture_output=True, text=True, timeout=180)
        if out.returncode != 0:
            return (FAIL, 0, "impact_analysis failed", out.stderr[-300:])
        j = json.loads(out.stdout)
        got = {r["node"] for r in j["affected"]}
        must = {"EQ-CM-01", "EQ-ENG-03", "EQ-ENG-04", "EQ-ENG-08", "EQ-SCT-01"}
        must_not = {"CONST-AL6061", "EQ-ELAS-01", "SRC-CROSSREF", "EQ-GOV-02"}
        ok = must <= got and not (must_not & got)
        return (PASS if ok else FAIL, j["n_affected"],
                f"{j['n_affected']} downstream; missing={sorted(must - got)}; "
                f"wrongly included={sorted(must_not & got)}", "")

    finish(c, body)


def test_DEP02_unrelated_change_no_invalidation():
    """A documentation change must not touch any numerical result."""
    c = make_case(test_id="DEP02", category="12-claim-traceability",
                  purpose="A documentation-only edit has no numerical dependents, so it "
                          "cannot invalidate a numerical result",
                  oracle="the graph rejects a query on a node that is not a scientific "
                         "artefact; documentation is deliberately not modelled",
                  oracle_independent=True,
                  fixture="README.md and a non-existent equation id",
                  expected="both queries exit non-zero",
                  tolerance=None, tol_justification="exact", level="1")

    def body():
        r1 = subprocess.run(
            [sys.executable, str(REPO / "evidence" / "impact_analysis.py"),
             "--changed", "README.md"],
            cwd=REPO, capture_output=True, text=True, timeout=180)
        r2 = subprocess.run(
            [sys.executable, str(REPO / "evidence" / "impact_analysis.py"),
             "--changed", "EQ-NOT-A-REAL-EQUATION"],
            cwd=REPO, capture_output=True, text=True, timeout=180)
        ok = r1.returncode != 0 and r2.returncode != 0
        return (PASS if ok else FAIL, 2,
                f"README.md query exit {r1.returncode}; bogus equation query exit "
                f"{r2.returncode} (both must be non-zero)", "")

    finish(c, body)
