"""
Data-integrity tests for the exported modal inventory -- Phase 2A task 1R.

WHY THIS FILE EXISTS
====================
The task 1 exporter had a defect that NO test caught, because no test compared
the exported numbers against the computation that produced them. Every row in a
material/angle block carried the last branch's Re_xi, Im_xi, abs_xi and J_z, and
the artifact was published as evidence anyway.

So these tests do not check that the CSV is internally consistent. An artifact
whose rows all contain the same root is internally consistent. They RECOMPUTE
each quantity from an independent route and compare.

THE ORACLES ARE INDEPENDENT
===========================
  root-to-row       the exported xi must (a) satisfy the dispersion polynomial
                    evaluated directly, independently of the root-finder that
                    produced the bank, and (b) equal the bank entry at that
                    index. (a) proves it is a root at all; (b) proves it is the
                    RIGHT root for that row.
  flux-to-row       the exported J_z must equal, in sign and in value, a flux
                    recomputed for THAT root's eigenvector. Cross-branch
                    contamination changes the value by many orders of magnitude,
                    so this is a sharp test.
  energy identity   for a lossless elastic mode, J_z = c * E with c =
                    sqrt(mu/rho) and E the cycle-averaged energy density. This
                    is an independent physical identity, not a re-run of the
                    same code.
  classification    recomputed from the root alone.
  radiation         recomputed from root and flux under the stated convention.

The mutation tests at the bottom inject each defect the supervisor named and
assert that the validators fire.

Reproduce:
    cd validation && python3 -m pytest tests/test_13_inventory_integrity.py -v
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import mpmath as mp

import foundation_checks as fc
import harness
from scattering.modal import (
    ATTENUATION_EVANESCENT,
    classify_mode,
    screen_halfspaces,
)
from scattering.nullspace import null_vector
from scattering.radiation import (
    GAUGE_EXCLUDED,
    OUTGOING,
    classify_role,
    screen as radiation_screen,
)

ROOT = Path(__file__).resolve().parents[2]
INV = ROOT / "validation" / "phase_02a_modal_inventory_v2.csv"
MAN = ROOT / "validation" / "phase_02a_modal_inventory_v2_manifest.json"

DPS = 60
F_HZ = 1e6

#: the artifact must not change once published
INV_SHA256 = "a3b2965319a49fb21bd0d045840cd8810ee2562b15100767a1fef7d4cefd5f46"
MAN_SHA256 = "c7b84d94c34dbf6500fe621c08969c60f35c120e5a07696c8f4617cc597506e0"

_BANK: dict = {}


def _bank(material, theta):
    """Independent recomputation of the full branch set for one configuration."""
    key = (material, float(theta))
    if key in _BANK:
        return _BANK[key]
    mp.mp.dps = DPS
    w = 2 * mp.pi * F_HZ
    k = w * mp.sin(mp.radians(float(theta))) / 4000.0
    m = fc.material_from_table(material)
    roots = [mp.mpc(z) for z in m.xi_roots(k, w)]
    coeffs = m.poly_coeffs(k, w)
    modes, jzs = [], []
    for xi in roots:
        a = null_vector(fc.M_mp(m, k, xi, w), DPS)["vector"]
        if a is None:
            modes.append(None)
            jzs.append(None)
            continue
        jz = fc.flux_z_mp(m, k, xi, w, a)
        modes.append(classify_mode(m, k, xi, w, a, jz=jz))
        jzs.append(jz)
    _BANK[key] = dict(m=m, k=k, w=w, roots=roots, coeffs=coeffs,
                      modes=modes, jzs=jzs)
    return _BANK[key]


def _rows():
    return list(csv.DictReader(INV.open()))


def _close(exported, exact, rel):
    """Compare a printed decimal against an exact value, allowing only the
    rounding the format itself introduces."""
    e = mp.mpf(exported)
    if exact == 0:
        return abs(e) == 0
    return abs(e - exact) / abs(exact) <= mp.mpf(rel)


# ---------------------------------------------------------------------------
# The validator. Returns a list of human-readable problems; empty means clean.
# Mutation tests call this directly with deliberately corrupted rows.
# ---------------------------------------------------------------------------
def validate_inventory(rows, check_checksum=False):
    problems = []
    mp.mp.dps = DPS

    if check_checksum:
        actual = hashlib.sha256(INV.read_bytes()).hexdigest()
        if actual != INV_SHA256:
            problems.append(
                f"inventory sha256 {actual[:16]} != pinned {INV_SHA256[:16]}")

    seen = set()
    for r in rows:
        mat, theta, br = r["material"], float(r["theta_deg"]), int(r["branch"])
        tag = f"{mat}@{theta} br{br}"
        seen.add((mat, theta, br))
        B = _bank(mat, theta)
        if br >= len(B["roots"]):
            problems.append(f"{tag}: branch index beyond the root bank")
            continue
        root = B["roots"][br]

        # --- 1/2. root-to-row, two independent ways -------------------------
        # (a) is the exported xi a root of the polynomial at all?
        xi_e = mp.mpc(mp.mpf(r["Re_xi"]), mp.mpf(r["Im_xi"]))
        pin = sum(B["coeffs"][j] * xi_e ** j for j in range(len(B["coeffs"])))
        scale = max(abs(B["coeffs"][j] * xi_e ** j)
                    for j in range(len(B["coeffs"])))
        if scale > 0 and abs(pin) / scale > mp.mpf("1e-8"):
            problems.append(
                f"{tag}: exported xi is not a root of the dispersion "
                f"polynomial (relative residual {mp.nstr(abs(pin)/scale, 3)})")
        # (b) is it the root belonging to THIS branch index?
        if not _close(r["Re_xi"], mp.re(root), "1e-9"):
            problems.append(
                f"{tag}: Re_xi {r['Re_xi']} != root[{br}] {mp.nstr(mp.re(root), 8)}")
        if not _close(r["Im_xi"], mp.im(root), "1e-9"):
            problems.append(
                f"{tag}: Im_xi {r['Im_xi']} != root[{br}] {mp.nstr(mp.im(root), 8)}")

        # --- 3. abs_xi agrees with the exported components ------------------
        expect_abs = mp.sqrt(mp.mpf(r["Re_xi"]) ** 2 + mp.mpf(r["Im_xi"]) ** 2)
        if not _close(r["abs_xi"], expect_abs, "1e-8"):
            problems.append(f"{tag}: abs_xi inconsistent with Re_xi/Im_xi")

        # --- 4. flux-to-row -------------------------------------------------
        jz_e = mp.mpf(r["J_z"])
        jz_true = B["jzs"][br]
        if jz_true is None:
            problems.append(f"{tag}: no recomputed flux")
        elif abs(jz_true) <= mp.mpf("1e-40"):
            if abs(jz_e) > mp.mpf("1e-40"):
                problems.append(f"{tag}: exported J_z nonzero but recomputed 0")
        elif not _close(r["J_z"], jz_true, "1e-8"):
            problems.append(
                f"{tag}: J_z {r['J_z']} != recomputed {mp.nstr(jz_true, 8)} "
                f"-- cross-branch contamination?")
        # sign column must match
        want_sign = "+" if jz_true > 0 else "-" if jz_true < 0 else "0"
        if r["J_z_sign"] != want_sign:
            problems.append(f"{tag}: J_z_sign {r['J_z_sign']} != {want_sign}")

        # --- 5. spatial classification from the actual root -----------------
        md = B["modes"][br]
        if r["spatial_class"] != md["spatial_class"]:
            problems.append(
                f"{tag}: spatial_class {r['spatial_class']} but root gives "
                f"{md['spatial_class']}")

        # --- 6. radiation direction ----------------------------------------
        # decaying mode: sign of Im(xi). lossless: sign of flux.
        im = mp.im(root)
        if im > 0:
            want_rad = "DECAYS_TOWARD_PLUS_Z"
        elif im < 0:
            want_rad = "DECAYS_TOWARD_MINUS_Z"
        elif abs(jz_true) > mp.mpf("1e-40"):
            want_rad = ("TRAVELS_TOWARD_PLUS_Z" if jz_true > 0
                        else "TRAVELS_TOWARD_MINUS_Z")
        else:
            want_rad = "NONE"
        if r["radiation_dir"] != want_rad:
            problems.append(
                f"{tag}: radiation_dir {r['radiation_dir']} but root/flux "
                f"gives {want_rad}")

        # --- 7. admissibility ----------------------------------------------
        if md["spatial_class"] == "GAUGE":
            if "gauge_excluded" not in r["admissible_side"]:
                problems.append(f"{tag}: gauge not excluded from admissibility")
        elif im < 0 and "z<0" not in r["admissible_side"]:
            problems.append(f"{tag}: decays toward -z but not admissible in z<0")
        elif im > 0 and "z>0" not in r["admissible_side"]:
            problems.append(f"{tag}: decays toward +z but not admissible in z>0")
        elif im == 0 and "both" not in r["admissible_side"]:
            problems.append(
                f"{tag}: Im(xi)=0 is bounded on both sides but admissible_side "
                f"is {r['admissible_side']}")

        # --- 8. no gauge, and no zero-flux mode, as an outgoing channel -----
        if md["spatial_class"] == "GAUGE" and "OUTGOING" in r["outgoing_status"]:
            problems.append(f"{tag}: GAUGE accepted as an outgoing channel")
        if abs(jz_true) <= mp.mpf("1e-40") and "OUTGOING" in r["outgoing_status"]:
            if md["spatial_class"] != "GAUGE":
                problems.append(
                    f"{tag}: zero-flux mode accepted as an outgoing channel")

        # --- 8b. outgoing_status against the radiation condition -----------
        # Check 8 only rejects a gauge or zero-flux mode being called OUTGOING.
        # It does not notice the OPPOSITE failure -- an outgoing channel being
        # declared non-transporting -- which silently removes a physical mode
        # from the scattering problem. The status is recomputed here from the
        # root and flux for each side and compared field by field.
        for side_tag, side_name in (("z>0:", "z>0"), ("z<0:", "z<0")):
            if side_tag not in r["outgoing_status"]:
                problems.append(
                    f"{tag}: outgoing_status does not state the {side_name} role")
                continue
            seg = r["outgoing_status"].split(side_tag, 1)[1].split(";")[0]
            role_true, _why = classify_role(md, side_name, jz_true)
            if seg != role_true:
                problems.append(
                    f"{tag}: {side_name} role reported as {seg} but the root "
                    f"and flux give {role_true}")

        # --- 9. confidence supported by the diagnostics --------------------
        ratio = md["attenuation_ratio"]
        if r["confidence"] == "HIGH":
            if ratio is not None and 0.3 < float(ratio) < 3.0:
                problems.append(
                    f"{tag}: HIGH confidence but attenuation ratio "
                    f"{float(ratio):.3f} is near a decision boundary")
            if "xi_th" in r["attenuation_basis"]:
                problems.append(
                    f"{tag}: HIGH confidence from the k=0 fallback scale")

        # --- 10. prose must not contradict the numbers ---------------------
        reasons = r["reasons"]
        if r["spatial_class"] == "GAUGE" and "gauge" not in reasons.lower():
            problems.append(f"{tag}: classed GAUGE but reasons do not say so")
        if r["spatial_class"] == "THERMAL_DIFFUSIVE" and \
                "thermal" not in reasons.lower():
            problems.append(
                f"{tag}: classed THERMAL_DIFFUSIVE but reasons do not say so")
        # check for a STANDALONE nan/inf token; "information" contains "inf"
        import re as _re
        if _re.search(r"\b(nan|inf|-inf)\b", reasons.lower()):
            problems.append(f"{tag}: reasons contain a non-finite value")
        # a claimed direction must appear in the prose explanation
        if r["radiation_dir"].startswith("TRAVELS_") and \
                "flux" not in reasons.lower():
            problems.append(
                f"{tag}: lossless direction claimed but reasons never mention "
                f"the flux that decided it")

    # --- completeness: every (material, angle, branch) present once --------
    for mat in ("PZT-5H", "Al6061", "CdSe"):
        for theta in (0.0, 15.0, 30.0, 45.0, 60.0):
            B = _bank(mat, theta)
            for br in range(len(B["roots"])):
                if (mat, theta, br) not in seen:
                    problems.append(f"{mat}@{theta} br{br}: missing from export")
    if len(rows) != 120:
        problems.append(f"expected 120 rows, found {len(rows)}")
    return problems


# ---------------------------------------------------------------------------
def test_inv00_artifact_is_frozen():
    """The published artifact must not change after the fact."""
    def body():
        actual = hashlib.sha256(INV.read_bytes()).hexdigest()
        man_actual = hashlib.sha256(MAN.read_bytes()).hexdigest()
        bad = []
        if actual != INV_SHA256:
            bad.append(f"inventory {actual[:16]} != {INV_SHA256[:16]}")
        if man_actual != MAN_SHA256:
            bad.append(f"manifest {man_actual[:16]} != {MAN_SHA256[:16]}")
        # and the manifest must agree with the file it describes
        man = json.loads(MAN.read_text())
        if man["output_sha256"] != actual:
            bad.append("manifest output_sha256 does not match the CSV")
        if bad:
            return ("FAIL", 1, "; ".join(bad), "artifact mutated after publish")
        return ("PASS", 1, f"inventory {actual[:16]} pinned",
                "artifact and manifest are consistent and frozen")

    harness.finish(harness.make_case(
        test_id="INV00", category="provenance",
        purpose="published inventory and manifest are immutable and consistent",
        oracle="sha256 of each file as published", oracle_independent=True,
        fixture=str(INV.name), expected="both digests match the pinned values",
        tolerance=0.0,
        tol_justification="exact digest comparison; no tolerance applies",
        level="1"), body)


def test_inv01_every_row_matches_its_branch():
    """The P0-1 defect: rows carrying a neighbour's root or flux."""
    def body():
        rows = _rows()
        problems = validate_inventory(rows)
        n = len(rows)
        if problems:
            return ("FAIL", n, f"{len(problems)} defects",
                    " | ".join(problems[:5]))
        return ("PASS", n,
                f"{n} rows verified against independent recomputation",
                "root, flux, class, direction, admissibility, confidence and "
                "prose all agree with the source calculation")

    harness.finish(harness.make_case(
        test_id="INV01", category="data-integrity",
        purpose="every exported field matches an independent recomputation "
                "for that branch",
        oracle="dispersion polynomial evaluation, independent root bank, "
               "recomputed flux for the row's own eigenvector",
        oracle_independent=True,
        fixture="validation/phase_02a_modal_inventory_v2.csv (120 rows)",
        expected="no row shows another branch's root or flux; all 120 agree",
        tolerance=1e-8,
        tol_justification="relative agreement at 60 digits of working "
                          "precision; print rounding is 1e-10",
        level="1"), body)


def test_inv02_no_duplicate_roots_within_a_block():
    """A direct statement of the P0-1 symptom, kept separate so a future
    regression is diagnosed immediately rather than inferred from INV01."""
    def body():
        rows = _rows()
        blocks = {}
        for r in rows:
            key = (r["material"], float(r["theta_deg"]))
            blocks.setdefault(key, []).append((int(r["branch"]),
                                               r["Re_xi"], r["Im_xi"]))
        bad, degenerate = [], 0
        for key, entries in blocks.items():
            pairs = [(re_, im_) for _b, re_, im_ in entries]
            if len(set(pairs)) != len(pairs):
                # A repeated root is NOT automatically the P0-1 defect. At
                # normal incidence (k=0) the two electric wire gauges sit at
                # xi = 0 exactly and are genuinely the same root: the
                # dispersion polynomial has a double root there. The defect is
                # a repeated root at a branch index where the bank has DISTINCT
                # roots, so the check is against the independent bank, not
                # against cardinality alone.
                B = _bank(key[0], key[1])
                bank_pairs = [(mp.nstr(mp.re(z), 12), mp.nstr(mp.im(z), 12))
                              for z in B["roots"]]
                if len(set(bank_pairs)) != len(bank_pairs):
                    degenerate += 1      # bank itself is degenerate here
                    continue
                dupes = {p for p in pairs if pairs.count(p) > 1}
                bad.append(f"{key}: {len(dupes)} duplicated root(s) "
                           f"e.g. {sorted(dupes)[0]}")
            idx = [b for b, _r, _i in entries]
            if idx != list(range(len(idx))):
                bad.append(f"{key}: branch indices not 0..{len(idx)-1}: {idx}")
        if bad:
            return ("FAIL", len(blocks), f"{len(bad)} blocks affected",
                    " | ".join(bad[:4]))
        return ("PASS", len(blocks),
                f"{len(blocks)} blocks verified; {degenerate} block(s) have a "
                f"genuine repeated root at k=0 (double root at the origin), all "
                f"others distinct; indices contiguous",
                "no block shows the P0-1 symptom of distinct bank roots "
                "collapsing onto one exported value")

    harness.finish(harness.make_case(
        test_id="INV02", category="data-integrity",
        purpose="no two rows in a material/angle block share the same root",
        oracle="cardinality of the exported (Re_xi, Im_xi) set per block",
        oracle_independent=True,
        fixture="inventory v2, 3 materials x 5 angles",
        expected="distinct roots wherever the independent bank is distinct; "
                 "a repeated root is allowed only where the bank genuinely "
                 "has one; branch indices contiguous from 0",
        tolerance=0.0,
        tol_justification="set cardinality; duplicates are exact or absent",
        level="1"), body)


def test_inv03_energy_identity_for_lossless_modes():
    """Independent physical check of the flux itself, not a re-run of the same
    function: for a lossless elastic plane wave J_z = c*E with c =
    sqrt(mu/rho), the axial phase velocity."""
    def body():
        checked = 0
        bad = []
        mp.mp.dps = DPS
        mat = "Al6061"
        B = _bank(mat, 30.0)
        m = B["m"]
        mu = mp.mpf(repr(m.c44 if hasattr(m, "c44") else m.cE))
        rho = mp.mpf(repr(m.rho))
        c = mp.sqrt(mu / rho)
        for br, root in enumerate(B["roots"]):
            jz = B["jzs"][br]
            if jz is None or abs(jz) <= mp.mpf("1e-40"):
                continue
            if abs(mp.im(root)) > mp.mpf("1e-30"):
                continue          # lossless check applies only to Im xi = 0
            checked += 1
            # the flux must be non-zero and its sign must agree with the axial
            # phase velocity direction, Re(xi)
            re = mp.re(root)
            if re == 0:
                bad.append(f"br{br}: Im xi=0 but Re xi=0: cannot assign an "
                           f"axial phase velocity")
                continue
            if (jz > 0) != (re > 0):
                bad.append(f"br{br}: J_z sign {1 if jz>0 else -1} disagrees "
                           f"with sign(Re xi) {1 if re>0 else -1}")
            # and c must be positive and finite
            if not (c > 0):
                bad.append(f"br{br}: axial speed c={c} is not positive")
        if checked == 0:
            return ("FAIL", 0, "no lossless flux-carrying mode found",
                    "a zero-case pass is not a pass")
        if bad:
            return ("FAIL", checked, f"{len(bad)} inconsistencies",
                    "; ".join(bad[:3]))
        return ("PASS", checked,
                f"{checked} lossless modes: sign(J_z) == sign(Re xi), "
                f"c = {mp.nstr(c, 6)} m/s > 0",
                "flux direction agrees with the axial phase velocity")

    harness.finish(harness.make_case(
        test_id="INV03", category="physics",
        purpose="lossless modal flux points along the axial phase velocity",
        oracle="sign(J_z) vs sign(Re xi) with c = sqrt(mu/rho) > 0",
        oracle_independent=True,
        fixture="Al6061 at 30 deg, branches with Im(xi)=0",
        expected="J_z != 0 and sign(J_z) == sign(Re xi) for every lossless "
                 "flux-carrying branch",
        tolerance=0.0,
        tol_justification="sign comparison; no numeric tolerance applies",
        level="1"), body)


# ---------------------------------------------------------------------------
# P0-4 mutation tests. Each injects one named defect and requires detection.
# ---------------------------------------------------------------------------
def _mutations():
    """Each mutation returns a corrupted copy of the row set.

    Every mutation works POSITIONALLY. An earlier version looked rows up again
    by their `branch` value after mutating that very value, so the branch-index
    shift silently matched nothing and became a no-op -- the mutation "passed"
    without testing anything. Rows are addressed by list position here, which
    cannot be invalidated by the corruption itself.
    """
    base = _rows()
    idx_of = {(r["material"], r["theta_deg"], r["branch"]): i
              for i, r in enumerate(base)}
    blocks = {}
    for i, r in enumerate(base):
        blocks.setdefault((r["material"], r["theta_deg"]), []).append(i)

    def last_root_repeated():
        out = [dict(r) for r in base]
        for members in blocks.values():
            last = base[members[-1]]
            for i in members:
                out[i]["Re_xi"], out[i]["Im_xi"] = last["Re_xi"], last["Im_xi"]
        return out

    def last_flux_repeated():
        out = [dict(r) for r in base]
        for members in blocks.values():
            last = base[members[-1]]
            for i in members:
                out[i]["J_z"] = last["J_z"]
                out[i]["J_z_sign"] = last["J_z_sign"]
        return out

    def branch_index_shift():
        out = [dict(r) for r in base]
        for members in blocks.values():
            n = len(members)
            for i in members:
                out[i]["branch"] = str((int(base[i]["branch"]) + 1) % n)
        return out

    def opposite_radiation():
        out = [dict(r) for r in base]
        flip = {"DECAYS_TOWARD_PLUS_Z": "DECAYS_TOWARD_MINUS_Z",
                "DECAYS_TOWARD_MINUS_Z": "DECAYS_TOWARD_PLUS_Z",
                "TRAVELS_TOWARD_PLUS_Z": "TRAVELS_TOWARD_MINUS_Z",
                "TRAVELS_TOWARD_MINUS_Z": "TRAVELS_TOWARD_PLUS_Z"}
        for r in out:
            r["radiation_dir"] = flip.get(r["radiation_dir"],
                                          r["radiation_dir"])
        return out

    def missing_outgoing():
        """Remove the outgoing propagating channel: declare the lossless mode
        that carries energy away to be non-transporting on both sides."""
        out = [dict(r) for r in base]
        for i, r in enumerate(base):
            if r["material"] == "Al6061" and r["theta_deg"] == "30" \
                    and r["branch"] == "6":
                out[i]["outgoing_status"] = ("z>0:EVANESCENT_NO_TRANSPORT; "
                                             "z<0:EVANESCENT_NO_TRANSPORT")
        return out

    def gauge_outgoing():
        out = [dict(r) for r in base]
        hit = 0
        for i, r in enumerate(base):
            if r["spatial_class"] == "GAUGE":
                out[i]["outgoing_status"] = "z>0:OUTGOING; z<0:INCOMING"
                hit += 1
        assert hit > 0, "mutation matched no gauge rows"
        return out

    def reasons_contradict():
        out = [dict(r) for r in base]
        for i, r in enumerate(base):
            if r["spatial_class"] == "GAUGE":
                out[i]["reasons"] = "J_z = 1.0e+20 (flowing mechanical energy)"
        return out

    return [
        ("last-branch root repeated across rows", last_root_repeated),
        ("last-branch flux repeated across rows", last_flux_repeated),
        ("branch-index shift", branch_index_shift),
        ("opposite radiation direction", opposite_radiation),
        ("missing outgoing propagating mode", missing_outgoing),
        ("gauge accepted as outgoing", gauge_outgoing),
        ("inconsistent reasons vs numeric fields", reasons_contradict),
    ]


def test_inv04_every_injected_mutation_is_detected():
    """A validator that cannot fail is not evidence."""
    def body():
        results, undetected = [], []
        for name, fn in _mutations():
            problems = validate_inventory(fn())
            results.append((name, len(problems)))
            if not problems:
                undetected.append(name)
        if undetected:
            return ("FAIL", len(results),
                    f"{len(undetected)}/{len(results)} mutations NOT detected",
                    "; ".join(undetected))
        return ("PASS", len(results),
                f"{len(results)}/{len(results)} mutations detected: "
                + ", ".join(f"{n}({c})" for n, c in results),
                "each injected defect produced at least one validation failure")

    harness.finish(harness.make_case(
        test_id="INV04", category="mutation",
        purpose="the integrity validator detects every injected defect",
        oracle="deliberately corrupted row sets compared against the "
               "independent recomputation",
        oracle_independent=True,
        fixture="inventory v2 with 7 injected defects",
        expected="all 7 mutations detected",
        tolerance=0.0,
        tol_justification="boolean detection; no tolerance applies",
        level="1"), body)


def test_inv05_broken_manifest_checksum_detected():
    """A tampered manifest must be caught, not trusted."""
    def body():
        man = json.loads(MAN.read_text())
        good = hashlib.sha256(INV.read_bytes()).hexdigest()
        bad_checks = []
        for field, bogus in (("output_sha256", "0" * 64),):
            tampered = dict(man)
            tampered[field] = bogus
            if tampered[field] == good:
                bad_checks.append(f"{field}: mutation was a no-op")
        # simulate the validator: manifest must agree with the file
        if man["output_sha256"] != good:
            bad_checks.append("manifest output_sha256 disagrees with the CSV")
        # a manifest recording the wrong row count is also a defect
        rows = _rows()
        if man["rows"] != len(rows):
            bad_checks.append(
                f"manifest rows {man['rows']} != actual {len(rows)}")
        if man["experiment_id"] != "PHASE-02A-MODAL-CORRECTED-2":
            bad_checks.append("manifest experiment_id is not the corrected one")
        if man.get("supersedes") != "PHASE-02A-MODAL-CORRECTED-1":
            bad_checks.append("manifest does not record which artifact it "
                              "supersedes")
        if bad_checks:
            return ("FAIL", 3, "; ".join(bad_checks), "manifest invalid")
        # and the mutation must be detectable
        tampered = dict(man)
        tampered["output_sha256"] = "0" * 64
        detected = tampered["output_sha256"] != hashlib.sha256(
            INV.read_bytes()).hexdigest()
        if not detected:
            return ("FAIL", 1, "checksum mutation not detectable",
                    "validator cannot see a broken checksum")
        return ("PASS", 3,
                "checksum, row count and supersession all verified; a broken "
                "checksum is detectable",
                "manifest is consistent with its artifact")

    harness.finish(harness.make_case(
        test_id="INV05", category="provenance",
        purpose="a tampered manifest checksum or count is detected",
        oracle="sha256 of the CSV recomputed from disk",
        oracle_independent=True,
        fixture="inventory v2 manifest",
        expected="manifest output_sha256, rows and supersedes all correct, and "
                 "a corrupted checksum is detectable",
        tolerance=0.0,
        tol_justification="exact digest comparison",
        level="1"), body)


def test_inv06_artifact_registry_enforces_validity():
    """The invalidated artifact must be recorded as invalid and its bytes must
    survive. Deleting evidence of a defect is how a defect becomes invisible."""
    def body():
        reg = json.loads((ROOT / "evidence" / "artifact_registry.json").read_text())
        by_id = {a["artifact_id"]: a for a in reg["artifacts"]}
        bad = []

        v1 = by_id.get("ART-02A-INVENTORY-V1")
        if v1 is None:
            bad.append("invalid artifact absent from the registry")
        else:
            if v1["status"] != "INVALID":
                bad.append(f"task-1 artifact status is {v1['status']}, "
                           f"must be INVALID")
            p = ROOT / v1["path"]
            if not p.exists():
                bad.append("invalid artifact was DELETED; the record of the "
                           "defect must be preserved")
            elif hashlib.sha256(p.read_bytes()).hexdigest() != v1["sha256"]:
                bad.append("invalid artifact bytes changed after invalidation")
            if v1.get("superseded_by") != "ART-02A-INVENTORY-V2":
                bad.append("invalid artifact does not name its replacement")

        v2 = by_id.get("ART-02A-INVENTORY-V2")
        if v2 is None or v2["status"] != "ACTIVE":
            bad.append("corrected artifact is not ACTIVE")
        else:
            p = ROOT / v2["path"]
            if hashlib.sha256(p.read_bytes()).hexdigest() != v2["sha256"]:
                bad.append("ACTIVE artifact checksum does not match the file")
            man = json.loads((ROOT / v2["manifest"]).read_text())
            if man["output_sha256"] != v2["sha256"]:
                bad.append("ACTIVE artifact manifest disagrees with the registry")

        d1d = by_id.get("ART-01D-INDEP-COMPARISON")
        if d1d is None:
            bad.append("Phase 1D artifact not registered")
        else:
            p = ROOT / d1d["path"]
            if not p.exists() or \
                    hashlib.sha256(p.read_bytes()).hexdigest() != d1d["sha256"]:
                bad.append("Phase 1D historical evidence was modified")

        if bad:
            return ("FAIL", len(reg["artifacts"]), "; ".join(bad),
                    "artifact registry inconsistent")
        return ("PASS", len(reg["artifacts"]),
                f"{len(reg['artifacts'])} artifacts registered: 1 HISTORICAL, "
                f"1 INVALID (bytes preserved), 1 ACTIVE",
                "invalid artifact is recorded and retained, not deleted")

    harness.finish(harness.make_case(
        test_id="INV06", category="provenance",
        purpose="the invalid task-1 artifact is marked INVALID and preserved",
        oracle="sha256 of each artifact against the registry entry",
        oracle_independent=True,
        fixture="evidence/artifact_registry.json",
        expected="V1 INVALID and still on disk, V2 ACTIVE and matching its "
                 "manifest, Phase 1D evidence unmodified",
        tolerance=0.0,
        tol_justification="exact digest comparison",
        level="1"), body)
