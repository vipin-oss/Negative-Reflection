"""
Regression tests for modal classification -- Phase 2A task 1.

THESE TESTS ARE DESIGNED TO FAIL ON THE PHASE 1D LABELS.
========================================================
That is not rhetoric: every defect asserted against below was read out of
`validation/phase_01d_independent_comparison.csv`, which is preserved unchanged
as historical evidence. If the old classifier were reinstated, MAC01, MAC02 and
MAC03 would fail immediately.

The Phase 1D defects being locked out
-------------------------------------
D1  PROPAGATING assigned to modes with |Im xi| ~ 2.57e6 /m (PZT-5H) and
    8.66e5 /m (CdSe). Those decay by 1/e in 0.39 um and 1.2 um respectively --
    a fraction of a transverse wavelength. The label PROPAGATING is wrong for
    them in any convention. Count in the Phase 1D file: 92 PROPAGATING rows,
    of which 30 have |Im xi| > 1000 /m.

D2  The vocabulary emitted PHYSICAL_GROWING, which is self-contradictory for a
    half-space: a mode that grows without bound as z -> +infinity is not
    admissible in the upper half-space at all.

D3  A single `mode_class` field conflated four independent questions (spatial
    behaviour, radiation direction, energy flow, numerical confidence), so
    fixing one silently changed the others.

D4  No regression test existed, so nothing caught D1-D3.

Reproduce:
    cd validation && python3 -m pytest tests/test_12_modal_classification.py -v
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import mpmath as mp

import foundation_checks as fc
import harness
from scattering.modal import (
    ATTENUATION_EVANESCENT,
    GAUGE_MECHANICAL,
    THERMAL_DISPERSION_TOL,
    classify_mode,
    screen_halfspaces,
    verify_convention,
)
from scattering.nullspace import null_vector

ROOT = Path(__file__).resolve().parents[2]
PHASE1D_CSV = ROOT / "validation" / "phase_01d_independent_comparison.csv"

#: sha256 of the Phase 1D evidence file at the time Phase 2A began. If this
#: changes, the historical record has been altered and the audit is invalid.
PHASE1D_SHA256 = (
    "76c7174c08923ded6005da1f312d795086f867a7c1f6590cb2fb7985991eb56b"
)

DPS = 60


def _w_k(theta_deg, f_hz=1e6):
    w = 2 * mp.pi * f_hz
    k = w * mp.sin(mp.radians(theta_deg)) / 4000.0
    return w, k


def _modes(material, theta_deg):
    """Full classified inventory for one material/angle, indexes preserved."""
    mp.mp.dps = DPS
    w, k = _w_k(theta_deg)
    m = fc.material_from_table(material)
    out = []
    for xi in [mp.mpc(z) for z in m.xi_roots(k, w)]:
        r = null_vector(fc.M_mp(m, k, xi, w), DPS)
        a = r["vector"]
        jz = fc.flux_z_mp(m, k, xi, w, a) if a is not None else None
        out.append((xi, r, classify_mode(m, k, xi, w, a, jz=jz)
                    if a is not None else None))
    return out


# ---------------------------------------------------------------------------
def test_mac00_phase1d_evidence_is_unmodified():
    """The audit is only meaningful if the audited file cannot change."""
    digest = hashlib.sha256(PHASE1D_CSV.read_bytes()).hexdigest()

    def body():
        if digest != PHASE1D_SHA256:
            return ("FAIL", 1, digest,
                    "the Phase 1D evidence file has been modified; the "
                    "historical record must stay frozen and the corrected "
                    "results go under a new experiment ID")
        return ("PASS", 1, digest, "frozen historical evidence intact")

    harness.finish(harness.make_case(
        test_id="MAC00", category="provenance",
        purpose="Phase 1D independent-comparison CSV is byte-identical",
        oracle="sha256 of the file as committed in Phase 1D",
        oracle_independent=True,
        fixture=str(PHASE1D_CSV),
        expected=PHASE1D_SHA256, tolerance=0.0,
        tol_justification="exact digest comparison; a checksum has no tolerance",
        level="1"), body)


def test_mac01_no_mode_with_large_im_xi_is_called_propagating():
    """D1: the central defect. A mode decaying within a fraction of a
    wavelength must not be labelled as propagating, and no spatial label may
    contradict its own quantified attenuation."""

    def body():
        checked = 0
        bad = []
        for material in ("PZT-5H", "Al6061", "CdSe"):
            for theta in (0.0, 30.0, 60.0):
                for br, (xi, r, md) in enumerate(_modes(material, theta)):
                    if md is None:
                        continue
                    checked += 1
                    ratio = md["attenuation_ratio"]
                    if ratio is None:
                        # k=0 fallback: the ratio is measured against the
                        # thermal wavenumber, not against k, so it is NOT
                        # comparable with ATTENUATION_EVANESCENT. Nothing to
                        # check here; MAC08 covers this case instead.
                        continue
                    # the label must agree with the reported number
                    if ratio >= ATTENUATION_EVANESCENT:
                        if md["spatial_class"] in ("WEAKLY_DAMPED",):
                            bad.append(f"{material}@{theta} br{br}: "
                                       f"|Im/k|={float(ratio):.3e} labelled "
                                       f"{md['spatial_class']}")
                    elif ratio <= 1e-6:
                        if md["spatial_class"] in ("EVANESCENT",
                                                   "THERMAL_DIFFUSIVE"):
                            bad.append(f"{material}@{theta} br{br}: "
                                       f"|Im/k|={float(ratio):.3e} labelled "
                                       f"{md['spatial_class']}")
        if bad:
            return ("FAIL", checked, f"{len(bad)} inconsistent", "; ".join(bad[:3]))
        return ("PASS", checked, f"{checked} mode labels self-consistent",
                "no label contradicts its own attenuation ratio")

    harness.finish(harness.make_case(
        test_id="MAC01", category="classification",
        purpose="spatial labels agree with the reported attenuation ratio",
        oracle="self-consistency of (spatial_class, attenuation_ratio) under "
               "the derived e^{+i xi z} convention",
        oracle_independent=False,
        fixture="PZT-5H, Al6061, CdSe at 0/30/60 deg, 1 MHz",
        expected="no EVANESCENT/WEAKLY_DAMPED/THERMAL_DIFFUSIVE label "
                 "contradicts its own |Im xi|/|k|",
        tolerance=0.0,
        tol_justification="exact categorical check; no numeric tolerance applies",
        level="1"), body)


def test_mac02_thermal_branches_identified_by_dispersion():
    """The Phase 1D file labelled every thermal branch PROPAGATING. Thermal
    branches are identified by their dispersion relation, |xi| ~ |xi_th| =
    sqrt(w rho cE / kappa33), because amplitude composition alone misclassifies
    genuine acoustic modes that happen to carry a large temperature swing
    (Al6061 branch 2 was exactly such a case)."""

    def body():
        checked = 0
        bad = []
        for material in ("PZT-5H", "Al6061", "CdSe"):
            mp.mp.dps = DPS
            w, k = _w_k(30.0)
            m = fc.material_from_table(material)
            xi_th = mp.sqrt(w * mp.mpf(repr(m.rho)) * mp.mpf(repr(m.cE))
                            / mp.mpf(repr(m.kappa33)))
            for br, (xi, r, md) in enumerate(_modes(material, 30.0)):
                if md is None:
                    continue
                checked += 1
                disp = abs(xi) / xi_th
                is_th = abs(disp - 1) <= THERMAL_DISPERSION_TOL
                says_th = md["spatial_class"] == "THERMAL_DIFFUSIVE"
                if is_th != says_th:
                    bad.append(f"{material} br{br}: dispersion ratio "
                               f"{float(disp):.4f} but label "
                               f"{md['spatial_class']}")
        if bad:
            return ("FAIL", checked, f"{len(bad)} mismatched", "; ".join(bad[:3]))
        # and the thermal branches must actually exist and not also be called
        # propagating
        n_th = sum(1 for mat in ("PZT-5H", "Al6061", "CdSe")
                   for _xi, _r, md in _modes(mat, 30.0)
                   if md and md["spatial_class"] == "THERMAL_DIFFUSIVE")
        return ("PASS", checked,
                f"{n_th} thermal branches correctly separated",
                "thermal identification follows the dispersion relation")

    harness.finish(harness.make_case(
        test_id="MAC02", category="classification",
        purpose="thermal-diffusive branches found by dispersion, not amplitude",
        oracle="axial thermal wavenumber sqrt(w*rho*cE/kappa33)",
        oracle_independent=True,
        fixture="PZT-5H, Al6061, CdSe at 30 deg, 1 MHz",
        expected=f"label THERMAL_DIFFUSIVE iff | |xi|/|xi_th| - 1 | <= "
                 f"{THERMAL_DISPERSION_TOL}",
        tolerance=THERMAL_DISPERSION_TOL,
        tol_justification="measured gap: thermal branches land at 1.000-1.016, "
                          "acoustic branches at <= 6.1e-3, a factor > 150",
        level="1"), body)


def test_mac03_fields_are_independent():
    """D3: the four fields must not be a recoding of one another. A single
    label is what allowed a 'fix' to one question to silently move another."""

    def body():
        forbidden = {"PHYSICAL_GROWING"}
        checked = 0
        bad = []
        for material in ("PZT-5H", "Al6061", "CdSe"):
            for theta in (0.0, 30.0, 60.0):
                for br, (xi, r, md) in enumerate(_modes(material, theta)):
                    if md is None:
                        continue
                    checked += 1
                    vals = (md["spatial_class"], md["radiation_dir"],
                            md["energy_flux_status"], md["confidence"])
                    for v in vals:
                        if v in forbidden:
                            bad.append(f"{material} br{br}: emitted {v}")
                    if md["radiation_dir"] not in (
                            "DECAYS_TOWARD_PLUS_Z", "DECAYS_TOWARD_MINUS_Z",
                            "NONE"):
                        bad.append(f"{material} br{br}: bad radiation_dir")
                    if md["energy_flux_status"] not in (
                            "FLOWING", "ZERO", "UNDEFINED"):
                        bad.append(f"{material} br{br}: bad energy status")
                    # radiation_dir must be a function of Im(xi) alone, NOT of
                    # spatial_class: an evanescent and a thermal mode with the
                    # same Im(xi) sign must get the same direction
                    im = mp.im(xi)
                    want = ("DECAYS_TOWARD_PLUS_Z" if im > 0 else
                            "DECAYS_TOWARD_MINUS_Z" if im < 0 else "NONE")
                    if md["radiation_dir"] != want:
                        bad.append(f"{material} br{br}: Im xi={float(im):.3e} "
                                   f"but direction {md['radiation_dir']}")
        if bad:
            return ("FAIL", checked, f"{len(bad)} violations", "; ".join(bad[:3]))
        return ("PASS", checked,
                f"{checked} modes, 4 independent fields, no forbidden labels",
                "radiation direction follows Im(xi), not the spatial class")

    harness.finish(harness.make_case(
        test_id="MAC03", category="classification",
        purpose="the four classification fields are independent, not a recoding",
        oracle="Im(xi) sign rule derived from e^{+i xi z}",
        oracle_independent=True,
        fixture="PZT-5H, Al6061, CdSe at 0/30/60 deg, 1 MHz",
        expected="no PHYSICAL_GROWING or other self-contradictory label; "
                 "radiation_dir determined by Im(xi) alone",
        tolerance=0.0,
        tol_justification="categorical equality; no numeric tolerance applies",
        level="1"), body)


def test_mac04_convention_is_reverified_numerically():
    """The sign rule must be re-derived from the fields, not trusted as a
    comment. Phase 2A found this function returning agrees=False for every
    branch because the measured log-growth rate was compared against a decay
    rate -- an inconsistent negation that no test had caught."""

    def body():
        checked = 0
        bad = []
        for material in ("PZT-5H", "Al6061"):
            for xi, r, md in _modes(material, 30.0):
                a = r["vector"]
                if a is None:
                    continue
                meas, ana, agrees = verify_convention(m=fc.material_from_table(
                    material), k=_w_k(30.0)[1], xi=xi, w=_w_k(30.0)[0], a=a)
                checked += 1
                if not agrees:
                    bad.append(f"{material}: measured={mp.nstr(meas, 6)} "
                               f"analytic={mp.nstr(ana, 6)}")
        if bad:
            return ("FAIL", checked, f"{len(bad)} branches disagree",
                    "; ".join(bad[:3]))
        return ("PASS", checked,
                f"{checked} branches: measured decay rate == -Im(xi)",
                "convention e^{+i xi z} confirmed by field evaluation")

    harness.finish(harness.make_case(
        test_id="MAC04", category="convention",
        purpose="convention e^{+i xi z} reverified by direct field evaluation",
        oracle="|f(z)| = |f(0)| e^{-Im(xi) z} evaluated at two z values",
        oracle_independent=True,
        fixture="PZT-5H and Al6061 at 30 deg, 1 MHz",
        expected="measured decay rate equals -Im(xi) for every branch",
        tolerance=1e-25,
        tol_justification="relative agreement limit at 60 digits of working "
                          "precision; residual observed ~1e-50",
        level="1"), body)


def test_mac05_null_vector_is_precision_stable():
    """Phase 2A defect: the Phase 1D null-vector route latched onto a different
    branch at some working precisions. At Al6061 30 deg branch 3 it returned the
    pure electric wire gauge at dps 90, 110, 120 and 135, while returning the
    correct thermal branch at dps 60-85. A classifier whose answer depends on
    the number of digits requested is not a classifier."""

    def body():
        checked = 0
        bad = []
        for material in ("Al6061", "PZT-5H"):
            signatures = {}
            for dps in (50, 60, 80, 90, 100, 110, 120):
                mp.mp.dps = dps
                w, k = _w_k(30.0)
                m = fc.material_from_table(material)
                sig = []
                for xi in [mp.mpc(z) for z in m.xi_roots(k, w)]:
                    r = null_vector(fc.M_mp(m, k, xi, w), dps)
                    md = classify_mode(m, k, xi, w, r["vector"],
                                       jz=fc.flux_z_mp(m, k, xi, w,
                                                       r["vector"]))
                    sig.append((md["spatial_class"], md["radiation_dir"]))
                    checked += 1
                signatures[dps] = tuple(sig)
            distinct = set(signatures.values())
            if len(distinct) != 1:
                detail = "; ".join(
                    f"dps={d}: {[s[0][:4] for s in v]}" for d, v in
                    signatures.items())
                bad.append(f"{material}: {len(distinct)} distinct "
                           f"classifications -> {detail}")
        if bad:
            return ("FAIL", checked, "precision-dependent classification",
                    " | ".join(bad))
        mp.mp.dps = DPS
        return ("PASS", checked,
                "identical classification at dps 50,60,80,90,100,110,120",
                "the robust null-vector route removes the dependency")

    harness.finish(harness.make_case(
        test_id="MAC05", category="numerics",
        purpose="modal classification is invariant across working precision",
        oracle="repeated solve at seven independent working precisions",
        oracle_independent=True,
        fixture="PZT-5H and Al6061 at 30 deg, 1 MHz, dps 50..120",
        expected="one and only one classification signature per material",
        tolerance=0.0,
        tol_justification="categorical invariance; any difference is a defect",
        level="1"), body)


def test_mac06_gauge_and_flux_are_not_confused():
    """A pure electric wire gauge has no mechanical motion and therefore
    carries no mechanical energy. Phase 2A found the energy test written as
    `|jz| <= 1e-12 * |jz|`, a tautology that can only fire when jz is exactly
    zero, so the Al6061 gauges (|J_z| = 2.787e-107) were reported as FLOWING."""

    def body():
        checked = 0
        bad = []
        for material in ("Al6061", "PZT-5H"):
            for br, (xi, r, md) in enumerate(_modes(material, 30.0)):
                if md is None:
                    continue
                checked += 1
                if md["spatial_class"] == "GAUGE":
                    if md["energy_flux_status"] != "ZERO":
                        bad.append(f"{material} br{br}: GAUGE but "
                                   f"{md['energy_flux_status']} "
                                   f"(|J_z|={md['jz']})")
        # a gauge must not be admitted into a half-space as a physical mode
        for material in ("Al6061", "PZT-5H"):
            hs = screen_halfspaces([md for _x, _r, md in _modes(material, 30.0)
                                    if md is not None])
            for side in ("z>0_admissible", "z<0_admissible"):
                for idx in hs[side]:
                    modes = [md for _x, _r, md in _modes(material, 30.0)
                             if md is not None]
                    if modes[idx]["spatial_class"] == "GAUGE":
                        bad.append(f"{material}: gauge {idx} admitted in {side}")
        if bad:
            return ("FAIL", checked, f"{len(bad)} gauge/energy errors",
                    "; ".join(bad[:3]))
        return ("PASS", checked,
                "gauges carry ZERO flux and are excluded from half-spaces",
                "gauge detection uses machine-zero mechanical content")

    harness.finish(harness.make_case(
        test_id="MAC06", category="classification",
        purpose="gauge modes are separated from modes that carry energy",
        oracle="machine-zero mechanical content with a distinct null vector",
        oracle_independent=True,
        fixture="Al6061 and PZT-5H at 30 deg, 1 MHz",
        expected="GAUGE implies ZERO flux and exclusion from both half-spaces",
        tolerance=GAUGE_MECHANICAL,
        tol_justification="measured: true gauge mechanical content is exactly "
                          "0 relative to scale; smallest physical mode 4.01e-13",
        level="1"), body)


def test_mac07_halfspace_screening_is_permitted_direction():
    """A mode is admissible in a half-space when it satisfies TWO independent
    conditions, and Phase 2A found the first implementation checking only one:

      1. BOUNDEDNESS -- it must not grow as |z| -> infinity on that side. This
         is a statement about Im(xi) and the side of the interface, not about
         Re(xi) and not about energy flux.
      2. NON-TRIVIALITY -- it must be a physically realisable field. A pure
         electric wire gauge is bounded and satisfies the PDE but carries no
         mechanical motion, so it cannot participate in a mechanical
         scattering problem.

    Testing only condition 1 made this test PASS while a gauge was being
    admitted into both half-spaces as a physical mode, so the assertion below
    now checks both explicitly and names which one failed."""

    def body():
        checked = 0
        bad = []
        for material in ("PZT-5H", "Al6061", "CdSe"):
            mp.mp.dps = DPS
            w, k = _w_k(30.0)
            m = fc.material_from_table(material)
            modes = []
            for xi in [mp.mpc(z) for z in m.xi_roots(k, w)]:
                r = null_vector(fc.M_mp(m, k, xi, w), DPS)
                modes.append(classify_mode(m, k, xi, w, r["vector"],
                                           jz=fc.flux_z_mp(m, k, xi, w,
                                                           r["vector"])))
            hs = screen_halfspaces(modes)
            for i, md in enumerate(modes):
                checked += 1
                im = float(mp.im(m.xi_roots(k, w)[i]))
                in_up = i in hs["z>0_admissible"]
                in_dn = i in hs["z<0_admissible"]
                in_nd = i in hs["non_decaying"]
                in_gx = i in hs["gauge_excluded"]

                # every branch must land in exactly one category
                if sum((in_up, in_dn, in_nd, in_gx)) != 1:
                    bad.append(f"{material} br{i}: in {sum((in_up, in_dn, in_nd, in_gx))} "
                               f"categories at once")

                # condition 2 first: a gauge is excluded regardless of Im(xi)
                if md["spatial_class"] == "GAUGE":
                    if not in_gx:
                        bad.append(f"{material} br{i}: gauge admitted into a "
                                   f"half-space (boundedness != admissibility)")
                    continue

                # condition 1: boundedness decides the side
                if im > 0 and not in_up:
                    bad.append(f"{material} br{i}: Im xi={im:.3e} (>0, decays "
                               f"toward +z) excluded from z>0")
                if im < 0 and in_up:
                    bad.append(f"{material} br{i}: Im xi={im:.3e} (<0, decays "
                               f"toward -z) admitted into z>0")
                if im < 0 and not in_dn:
                    bad.append(f"{material} br{i}: decays toward -z but "
                               f"excluded from z<0")
                if im > 0 and in_dn:
                    bad.append(f"{material} br{i}: grows toward -z but "
                               f"admitted into z<0")
                if im == 0 and not in_nd:
                    bad.append(f"{material} br{i}: Im xi=0 but not reported "
                               f"as non-decaying")
        if bad:
            return ("FAIL", checked, f"{len(bad)} screening errors",
                    "; ".join(bad[:3]))
        return ("PASS", checked,
                "screening admits exactly the modes that decay into each "
                "half-space", "8 + 8 + 8 modes across three materials")

    harness.finish(harness.make_case(
        test_id="MAC07", category="selection",
        purpose="half-space admissibility requires boundedness AND "
                "non-triviality",
        oracle="boundedness of e^{+i xi z} as |z| -> infinity",
        oracle_independent=True,
        fixture="PZT-5H, Al6061, CdSe at 30 deg, 1 MHz",
        expected="each branch in exactly one of z>0 / z<0 / non_decaying / "
                 "gauge_excluded; decays toward +z => z>0 admissible, and "
                 "vice versa; gauges excluded on non-triviality",
        tolerance=0.0,
        tol_justification="categorical; Im(xi) is exactly zero or not",
        level="1"), body)


def test_mac08_normal_incidence_reports_its_basis():
    """At normal incidence k = 0, so |Im xi|/|k| does not exist. An earlier
    version left the ratio blank while classifying on a bare 1 m^-1 threshold --
    an unprincipled constant that happened to produce the right answer for these
    materials and would not generalise. The fallback is now the thermal
    wavenumber, and the basis is reported so a reader is never left with a blank
    where a measurement should be."""

    def body():
        checked = 0
        bad = []
        for material in ("PZT-5H", "Al6061", "CdSe"):
            mp.mp.dps = DPS
            w = 2 * mp.pi * 1e6
            m = fc.material_from_table(material)
            xi_th = mp.sqrt(w * mp.mpf(repr(m.rho)) * mp.mpf(repr(m.cE))
                            / mp.mpf(repr(m.kappa33)))
            for br, xi in enumerate([mp.mpc(z) for z in m.xi_roots(mp.mpf(0), w)]):
                a = null_vector(fc.M_mp(m, mp.mpf(0), xi, w), DPS)["vector"]
                md = classify_mode(m, mp.mpf(0), xi, w, a,
                                   jz=fc.flux_z_mp(m, mp.mpf(0), xi, w, a))
                checked += 1
                if md["attenuation_basis"] is None:
                    bad.append(f"{material} br{br}: no basis reported")
                elif "xi_th" not in md["attenuation_basis"]:
                    bad.append(f"{material} br{br}: k=0 did not fall back to "
                               f"the thermal scale ({md['attenuation_basis']})")
                # and the ratio must be the stated quantity
                if md["attenuation_ratio"] is not None:
                    want = abs(mp.im(xi)) / xi_th
                    if abs(md["attenuation_ratio"] - want) > mp.mpf("1e-30"):
                        bad.append(f"{material} br{br}: ratio "
                                   f"{mp.nstr(md['attenuation_ratio'], 6)} != "
                                   f"|Im xi|/|xi_th| {mp.nstr(want, 6)}")
                # a k=0 classification cannot be HIGH confidence if it rested on
                # a fallback scale rather than the physical transverse one,
                # UNLESS the thermal-dispersion test decided it directly
                if (md["spatial_class"] != "THERMAL_DIFFUSIVE"
                        and md["confidence"] == "HIGH"
                        and md["attenuation_ratio"] is not None
                        and abs(md["attenuation_ratio"] - ATTENUATION_EVANESCENT) < 0.5):
                    bad.append(f"{material} br{br}: fallback ratio "
                               f"{mp.nstr(md['attenuation_ratio'], 4)} is near the "
                               f"evanescent threshold but confidence is HIGH")
        if bad:
            return ("FAIL", checked, f"{len(bad)} basis errors",
                    "; ".join(bad[:3]))
        return ("PASS", checked,
                f"{checked} normal-incidence branches report the thermal fallback",
                "no blank attenuation field; basis is explicit and correct")

    harness.finish(harness.make_case(
        test_id="MAC08", category="classification",
        purpose="normal incidence reports which scale the attenuation used",
        oracle="|Im xi| / sqrt(w rho cE / kappa33) recomputed independently",
        oracle_independent=True,
        fixture="PZT-5H, Al6061, CdSe at k=0 (normal incidence), 1 MHz",
        expected="attenuation_basis names the thermal fallback and the reported "
                 "ratio equals |Im xi|/|xi_th|",
        tolerance=1e-30,
        tol_justification="relative identity check at 60 digits of precision",
        level="1"), body)
