"""
Material-data provenance -- task 5 of Phase 1D.

The defect this file closes: `build_materials()` used to contain hard-coded
literals while `literature/material_constants.csv` was documentation that
nothing read. Two sources for the same numbers, free to diverge, with no
guarantee that a reported result used the versioned one. A validator that only
*parses* the CSV would not fix it either, because parsing does not prevent the
literals from drifting away.

What is pinned here:
  DAT01  the table is complete and correctly labelled for every material;
  DAT02  the frozen literal fallback equals the table-loaded values EXACTLY,
         so the two sources cannot drift silently;
  DAT03  a missing critical constant raises rather than defaulting to zero --
         the check that stops an invented pyroelectric/thermoelastic/relaxation
         constant reaching a result;
  DAT04  a missing unit or an unlabelled status raises;
  DAT05  ASSUMED constants are flagged, so a quantitative claim can refuse to
         rest on them, and no ASSUMED constant is silently promoted.

Reproduce:
    cd validation && python3 -m pytest tests/test_08_provenance.py -v
"""
from __future__ import annotations

import csv
import pathlib
import shutil
import tempfile

from harness import REPO, PASS, FAIL, make_case, finish

import foundation_checks as fc

CSV = REPO / "literature" / "material_constants.csv"


# ===========================================================================
def test_PRV01_table_complete_and_labelled():
    c = make_case(
        test_id="PRV01", category="1-data-integrity",
        purpose="Every material in the production set comes from the versioned table "
                "with a unit, a source and a controlled status for every constant",
        oracle="literature/material_constants.csv read directly with the csv module, "
               "independent of foundation_checks.load_material_table",
        oracle_independent=True,
        fixture="literature/material_constants.csv",
        expected=">= 4 materials; every required symbol present for each; every row "
                 "carries a non-empty unit, source and a status in the controlled "
                 "vocabulary; 0 violations",
        tolerance=None, tol_justification="exact set and vocabulary comparison",
        level="1")

    def body():
        rows = list(csv.DictReader(CSV.open(newline="")))
        bad = []
        for r in rows:
            for field in ("unit", "source", "status"):
                if not (r.get(field) or "").strip():
                    bad.append(f"{r['material']}/{r['symbol']}: empty {field}")
            if fc._normalise_status(r["status"]) is None:
                bad.append(f"{r['material']}/{r['symbol']}: "
                           f"status {r['status']!r} not in vocabulary")
        mats = {r["material"] for r in rows}
        by_mat = {m: {r["symbol"] for r in rows if r["material"] == m}
                  for m in mats}
        missing = {m: [s for s in fc.REQUIRED_SYMBOLS if s not in by_mat[m]]
                   for m in mats}
        missing = {m: v for m, v in missing.items() if v}
        return rows, mats, bad, missing

    def judge(res):
        rows, mats, bad, missing = res
        detail = (f"{len(rows)} rows, materials: {sorted(mats)}\n"
                  f"violations: {len(bad)}\nmissing required: {missing}")
        if bad:
            return FAIL, f"{len(bad)} labelling violation(s)", detail + "\n" + "\n".join(bad[:10])
        if missing:
            return FAIL, f"incomplete materials: {sorted(missing)}", detail
        if len(mats) < 4:
            return FAIL, f"only {len(mats)} materials", detail
        return PASS, f"{len(rows)} rows, {len(mats)} materials, complete and labelled", detail

    def run():
        status, actual, detail = judge(body())
        return status, len(fc.REQUIRED_SYMBOLS), actual, detail

    finish(c, run)


# ===========================================================================
def test_PRV02_fallback_equals_table():
    """The frozen literals must equal the table-loaded values exactly.

    This is the anti-drift check. `_build_materials_frozen()` exists only so a
    checkout without ../literature can still run; the moment the two disagree,
    the mismatch means one of them was edited without the other, and a result
    could depend on which path was taken.
    """
    c = make_case(
        test_id="PRV02", category="1-data-integrity",
        purpose="The frozen literal fallback and the versioned table agree exactly, "
                "so the table is genuinely the single source of truth",
        oracle="bit-for-bit comparison of the two construction paths; a third "
               "independent statement (the CSV text itself) is parsed in DAT01",
        oracle_independent=True,
        fixture="literature/material_constants.csv vs _build_materials_frozen()",
        expected="identical float value for all 19 coupled constants of all 4 "
                 "materials; 0 mismatches",
        tolerance=0.0,
        tol_justification="both paths express the same decimal literals, so the "
                          "difference must be exactly zero; any non-zero value "
                          "means the sources have drifted",
        level="1")

    def body():
        table = fc.load_material_table()
        frozen = fc._build_materials_frozen()
        fields = ("C11", "C12", "C13", "C33", "C44", "e31", "e33", "e15",
                  "eps11", "eps33", "rho", "cE", "kappa11", "kappa33",
                  "tauq", "T0", "beta1", "beta3", "p")
        diffs, n = [], 0
        for nm in table:
            a = fc.material_from_table(nm, table)
            b = frozen[nm]
            for f in fields:
                n += 1
                va, vb = float(getattr(a, f)), float(getattr(b, f))
                if va != vb:
                    diffs.append(f"{nm}.{f}: table={va!r} literal={vb!r}")
        return diffs, n

    def judge(res):
        diffs, n = res
        detail = (f"compared {n} constant values across the two paths\n"
                  f"mismatches: {len(diffs)}")
        if diffs:
            return FAIL, f"{len(diffs)} of {n} values differ", detail + "\n" + "\n".join(diffs)
        return PASS, f"all {n} values identical", detail

    def run():
        status, actual, detail = judge(body())
        return status, 19, actual, detail

    finish(c, run)


# ===========================================================================
def _write_mutated_csv(mutate):
    """Copy the real CSV to a temp file, apply mutate(rows), return the path."""
    tmpdir = pathlib.Path(tempfile.mkdtemp(prefix="dat_"))
    dst = tmpdir / "material_constants.csv"
    with CSV.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    rows = mutate(rows)
    with dst.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return dst, tmpdir


def test_PRV03_missing_constant_blocks():

    c = make_case(
        test_id="PRV03", category="1-data-integrity",
        purpose="A missing critical constant raises instead of being filled in, so an "
                "invented pyroelectric / thermoelastic / relaxation constant cannot "
                "reach a result",
        oracle="direct construction of a corrupted table (one row deleted) and an "
               "assertion that load_material_table raises ProvenanceError",
        oracle_independent=True,
        fixture="a temporary copy of material_constants.csv with PZT-5H/p3 removed",
        expected="ProvenanceError raised, naming the missing symbol; the pristine "
                 "table still loads afterwards (negative control)",
        tolerance=None, tol_justification="exception presence, not a numeric tolerance",
        level="1")

    def body():
        # remove the ASSUMED pyroelectric row: exactly the constant the review
        # forbids inventing
        dst, tmpdir = _write_mutated_csv(
            lambda rows: [r for r in rows
                          if not (r["material"] == "PZT-5H" and r["symbol"] == "p3")])
        raised = None
        try:
            fc.load_material_table(str(dst))
        except fc.ProvenanceError as exc:
            raised = str(exc)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
        # negative control: the pristine table must still load
        pristine_ok = True
        try:
            fc.load_material_table()
        except Exception:                                       # noqa: BLE001
            pristine_ok = False
        return raised, pristine_ok

    def judge(res):
        raised, pristine_ok = res
        detail = (f"corrupted table raised: {raised!r}\n"
                  f"pristine table still loads: {pristine_ok}")
        if not pristine_ok:
            return FAIL, "the negative control failed: the real table no longer loads", detail
        if raised is None:
            return FAIL, ("a missing critical constant was accepted -- an invented "
                          "value could reach a result"), detail
        if "p" not in raised:
            return FAIL, f"raised, but the message does not name the symbol: {raised}", detail
        return PASS, "missing critical constant blocked", detail

    def run():
        status, actual, detail = judge(body())
        return status, 2, actual, detail

    finish(c, run)


# ===========================================================================
def test_PRV04_unlabelled_rejected():

    c = make_case(
        test_id="PRV04", category="1-data-integrity",
        purpose="A constant with no unit, no source, or an unrecognised status is "
                "rejected rather than silently accepted",
        oracle="three independently corrupted copies of the table, each broken in "
               "one way; the validator must reject all three",
        oracle_independent=True,
        fixture="temporary copies with (a) empty unit, (b) empty source, "
                "(c) status 'PROBABLY_FINE'",
        expected="ProvenanceError for all three corruptions; pristine table loads",
        tolerance=None, tol_justification="exception presence",
        level="1")

    def body():
        cases = {}

        def setfield(mat, sym, field, value):
            def f(rows):
                out = []
                for r in rows:
                    r = dict(r)
                    if r["material"] == mat and r["symbol"] == sym:
                        r[field] = value
                    out.append(r)
                return out
            return f

        for label, mut in (
                ("empty unit", setfield("PZT-5H", "C11", "unit", "")),
                ("empty source", setfield("PZT-5H", "C11", "source", "")),
                ("bad status", setfield("PZT-5H", "C11", "status", "PROBABLY_FINE")),
        ):
            dst, tmpdir = _write_mutated_csv(mut)
            try:
                fc.load_material_table(str(dst))
                cases[label] = "ACCEPTED (bad)"
            except fc.ProvenanceError as exc:
                cases[label] = f"rejected: {str(exc)[:70]}"
            except Exception as exc:                            # noqa: BLE001
                cases[label] = f"wrong exception {exc!r}"
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)
        return cases

    def judge(cases):
        detail = "\n".join(f"  {k}: {v}" for k, v in cases.items())
        bad = [k for k, v in cases.items() if not v.startswith("rejected")]
        if bad:
            return FAIL, f"{len(bad)} corruption(s) accepted: {bad}", detail
        return PASS, f"all {len(cases)} corruptions rejected", detail

    def run():
        status, actual, detail = judge(body())
        return status, 3, actual, detail

    finish(c, run)


# ===========================================================================
def test_PRV05_assumed_flagged():

    c = make_case(
        test_id="PRV05", category="1-data-integrity",
        purpose="ASSUMED constants are flagged on the medium so a quantitative claim "
                "can refuse to rest on them, and none is silently promoted",
        oracle="the status column of the table itself, compared against the "
               "'assumed' list attached to each Medium6mm",
        oracle_independent=True,
        fixture="literature/material_constants.csv",
        expected="the assumed list on each medium equals the set of rows marked "
                 "ASSUMED for that material; at least one material has an ASSUMED "
                 "constant, so the mechanism is exercised rather than vacuous",
        tolerance=None, tol_justification="exact set comparison",
        level="1")

    def body():
        table = fc.load_material_table()
        out = {}
        for nm, d in table.items():
            medium = fc.material_from_table(nm, table)
            got = set(getattr(medium, "assumed", []) or [])
            want = {fc._SYMBOL_ALIASES.get(s, s) for s in d["assumed"]}
            out[nm] = (got, want)
        return out

    def judge(out):
        lines, bad = [], []
        for nm, (got, want) in out.items():
            lines.append(f"  {nm}: flagged={sorted(got)} expected={sorted(want)}")
            if got != want:
                bad.append(nm)
        n_assumed = sum(len(w) for _, w in out.values())
        detail = "\n".join(lines) + f"\nassumed constants in the set: {n_assumed}"
        if bad:
            return FAIL, f"mismatch for {bad}", detail
        if n_assumed == 0:
            return FAIL, ("no material carries an ASSUMED constant, so this check "
                          "cannot detect a promotion"), detail
        return PASS, f"{n_assumed} ASSUMED constants correctly flagged", detail

    def run():
        status, actual, detail = judge(body())
        return status, 4, actual, detail

    finish(c, run)


# ===========================================================================
def test_PRV06_table_values_are_physically_admissible():
    """A unit slip in the table must break a PHYSICAL bound, not just a
    cross-source comparison.

    WHY THIS EXISTS (Phase 1D finding). Test DAT01 compares the code against the
    table. That was a real check while build_materials() hard-coded literals.
    After Phase 1D production LOADS the table, so DAT01 compares an object
    against the very rows it was built from: it is now a TAUTOLOGY and can no
    longer detect a corrupted table. Mutation M06 (a 1000x unit slip in
    PZT-5H C11) is the demonstration: it survives DAT01 by construction.

    The honest fix is not to weaken M06 but to add an independent check that
    does not consult the code: dimensional plausibility of the table itself,
    against externally known physical ranges for these material classes. This
    will not catch an arbitrary small edit -- and that limitation is stated
    rather than papered over -- but it does catch the unit-slip class, which is
    the one that silently changes every dimensionless group.
    """
    c = make_case(
        test_id="PRV06", category="1-data-integrity",
        purpose="Every tabulated constant lies inside a physically admissible range for "
                "its material class, so a unit slip cannot pass unnoticed",
        oracle="externally stated order-of-magnitude bounds for metals, "
               "piezoceramics and II-VI semiconductors, written in this test and "
               "derived from handbook ranges rather than from the code",
        oracle_independent=True,
        fixture="literature/material_constants.csv",
        expected="every value inside its declared band; C11 of each solid within "
                 "1e9..1e12 Pa; longitudinal speed within 1000..12000 m/s",
        tolerance=None,
        tol_justification="band membership is a boolean; the bands are decades wide, "
                          "so they reject unit slips (1000x) while admitting genuine "
                          "material variation",
        level="2")

    def body():
        table = fc.load_material_table()
        # decades-wide bands: [min, max] per symbol, PHYSICS not code
        BANDS = {
            "C11": (1e9, 1e12), "C12": (1e8, 1e12), "C13": (1e8, 1e12),
            "C33": (1e9, 1e12), "C44": (1e8, 1e12),
            "e31": (-1e2, 1e2), "e33": (-1e2, 1e2), "e15": (-1e2, 1e2),
            "eps11": (1e-12, 1e-7), "eps33": (1e-12, 1e-7),
            "rho": (1e2, 2e4),
            "cE": (1e2, 2e3),
            "kappa11": (1e-2, 1e3), "kappa33": (1e-2, 1e3),
            "tauq": (0.0, 1e-8),
            "T0": (1.0, 2000.0),
            "beta1": (1e2, 1e8), "beta3": (1e2, 1e8),
            "p3": (-1e-2, 1e-2),
        }
        bad, n = [], 0
        for mat, d in table.items():
            for sym, val in d["values"].items():
                if sym not in BANDS:
                    continue
                lo, hi = BANDS[sym]
                n += 1
                if not (lo <= val <= hi):
                    bad.append(f"{mat}/{sym} = {val:g} outside [{lo:g}, {hi:g}]")
            # an independent derived quantity: the longitudinal speed
            C11, C33, rho = d["values"]["C11"], d["values"]["C33"], d["values"]["rho"]
            for axis, C in (("x", C11), ("z", C33)):
                v = (C / rho) ** 0.5
                n += 1
                if not (1e3 <= v <= 1.2e4):
                    bad.append(f"{mat}: c_L({axis}) = {v:.0f} m/s implausible")
        return bad, n

    def judge(res):
        bad, n = res
        detail = f"checked {n} values against decade-wide physical bands\n" + \
                 "\n".join(bad[:15])
        if bad:
            return FAIL, f"{len(bad)} implausible value(s) -- likely a unit slip", detail
        return PASS, f"all {n} values inside their physical bands", detail

    def run():
        status, actual, detail = judge(body())
        return status, 1, actual, detail

    finish(c, run)
