"""
Phase 2A task 2 -- scattering assembly, classical validation and mode selection.

These checks cover the equations introduced by task 2:

  EQ-02A-STATE-01    physical interface state and its dimensional scales
  EQ-02A-ELASREF-01  independent isotropic elastic P-SV reference
  EQ-02A-ASM-01      two-half-space assembly with opposite-sign continuity
  EQ-02A-SOL-01      coupled solve, interface residuals, rescaling invariance
  EQ-02A-CLASS-01    gauge versus decoupled electrostatics
  EQ-02A-RAD-02      scale-aware propagation criterion

Every check states its assumptions, its oracle and why the oracle is
independent. The coupled amplitudes are labelled DIAGNOSTIC throughout: G4 and
G5 remain OPEN and nothing here closes them.
"""
from __future__ import annotations

import sys

import mpmath as mp

from harness import PASS, FAIL, make_case, finish
import foundation_checks as fc

from scattering.assembly import (                       # noqa: E402
    COND_TRUST, assemble, assemble_reduced_elastic, matrix_rank,
)
from scattering.elastic_reference import ISO, solve_psv  # noqa: E402
from scattering.modal import classify_mode               # noqa: E402
from scattering.nullspace import null_vector             # noqa: E402
from scattering.radiation import (                       # noqa: E402
    ATTENUATION_PROPAGATING, classify_role,
)
from scattering.solver import (                          # noqa: E402
    EVIDENCE_LABEL, analyse, boundary_residuals, solve,
)

W = 2 * mp.pi * 1e6
DPS = 40


# ===========================================================================
def test_SCTA01_reduced_elastic_matches_independent_oracle():
    """The 4x4 welded assembly must reproduce the independent reference."""
    c = make_case(
        test_id="SCTA01", category="10-classical-elastic",
        purpose="the assembly's modal machinery reproduces the independently "
                "derived isotropic P-SV reference",
        oracle="elastic_reference.solve_psv, itself validated against "
               "closed-form Fresnel, Snell, acoustic and grazing limits",
        oracle_independent=True,
        fixture="Al6061 -> steel-like, 0/10/20/30/45/60/75 deg, "
                "omega = 2*pi*1e6 rad/s",
        expected="all four complex amplitudes agree to < 1e-40",
        tolerance=1e-40,
        tol_justification="both sides use 50-digit arithmetic, so agreement is "
                          "at rounding; the measured worst case is 5e-42",
        level="1")

    def body():
        mp.mp.dps = 50
        A = ISO("Al6061", 2700.0, 107.6e9, 26.5e9)
        B = ISO("steel-like", 7800.0, 281.0e9, 82.0e9)
        worst, n, bad = mp.mpf(0), 0, []
        for th in (0, 10, 20, 30, 45, 60, 75):
            asm = assemble_reduced_elastic(A, B, W, mp.mpf(th))
            ref = solve_psv(A, B, W, mp.mpf(th))
            for ch in ("rP", "rS", "tP", "tS"):
                n += 1
                d = abs(asm["amps"][ch] - ref[ch])
                worst = max(worst, d)
                if d > mp.mpf("1e-40"):
                    bad.append(f"th={th} {ch} diff {mp.nstr(d, 4)}")
            if asm["rank"] != 4:
                bad.append(f"th={th} rank {asm['rank']} != 4")
        if bad:
            return (FAIL, n, f"worst {mp.nstr(worst, 4)}", " | ".join(bad[:4]))
        return (PASS, n, f"worst difference {mp.nstr(worst, 4)}",
                "the assembly builds its own columns and selection rule; "
                "nothing is imported from solve_psv")

    finish(c, body)


# ===========================================================================
def test_SCTA02_dof_counted_not_assumed():
    """The system size must follow the modal inventory, never a hard-coded 8."""
    c = make_case(
        test_id="SCTA02", category="6-interface-conditions",
        purpose="unknowns equal the admitted-mode count and gauges are excluded",
        oracle="the assembly's own modal inventory, recomputed from the roots",
        oracle_independent=True,
        fixture="Al6061, PZT-5H and Al6061|PZT-5H at 0/30/60 deg",
        expected="n_unknowns == admitted side1 + admitted side2; no gauge admitted",
        tolerance=0.0,
        tol_justification="exact integer comparison",
        level="1")

    def body():
        mp.mp.dps = DPS
        n, bad = 0, []
        A_al = fc.material_from_table("Al6061")
        A_pz = fc.material_from_table("PZT-5H")
        for nm1, m1, nm2, m2 in (("Al6061", A_al, "Al6061", A_al),
                                 ("PZT-5H", A_pz, "PZT-5H", A_pz),
                                 ("Al6061", A_al, "PZT-5H", A_pz)):
            for th in (0, 30, 60):
                k = W * mp.sin(mp.radians(mp.mpf(th))) / mp.mpf(4000)
                S = assemble(m1, m2, k, W)
                n += 1
                want = len(S["admitted1"]) + len(S["admitted2"])
                if S["n_unknowns"] != want:
                    bad.append(f"{nm1}|{nm2}@{th}: unk {S['n_unknowns']} "
                               f"!= {want}")
                for sd, adm in ((S["side1"], S["admitted1"]),
                                (S["side2"], S["admitted2"])):
                    gauges = [b for b, md in enumerate(sd.modes)
                              if md is not None
                              and md["spatial_class"] == "GAUGE"]
                    leaked = [b for b in gauges if b in adm]
                    if leaked:
                        bad.append(f"{nm1}|{nm2}@{th}: gauge {leaked} admitted")
        if bad:
            return (FAIL, n, f"{len(bad)} defects", " | ".join(bad[:4]))
        return (PASS, n, "counts follow the inventory",
                "a square 8x8 appears only where the physics supplies it")

    finish(c, body)


# ===========================================================================
def test_SCTA03_gauge_vs_electrostatic():
    """The defect this task found: a decoupled electrostatic mode is physical."""
    c = make_case(
        test_id="SCTA03", category="6-interface-conditions",
        purpose="a mechanically decoupled mode with D_z != 0 is ELECTROSTATIC "
                "and admissible, not a GAUGE",
        oracle="D_z = -eps33*xi*phi evaluated from the medium constants and "
               "the root, independent of any component-norm comparison",
        oracle_independent=True,
        fixture="Al6061 at 30 deg (eps33 = eps0, no piezoelectric coupling)",
        expected="D_z = -eps33*k = -6.954e-9, mode admitted",
        tolerance=1e-20,
        tol_justification="D_z is computed two ways from the same exact inputs; "
                          "agreement is at rounding",
        level="1")

    def body():
        mp.mp.dps = DPS
        A = fc.material_from_table("Al6061")
        k = W * mp.sin(mp.radians(mp.mpf(30))) / mp.mpf(4000)
        S = assemble(A, A, k, W)
        sd = S["side1"]
        electro = [b for b, md in enumerate(sd.modes)
                   if md is not None and md["spatial_class"] == "ELECTROSTATIC"]
        outgoing = [b for b in electro if sd.roles[b][0] == "OUTGOING"]
        n = len(electro)
        bad = []
        if not electro:
            return (FAIL, 0, "no ELECTROSTATIC mode found",
                    "the xi = -/+ i k pair should be classified ELECTROSTATIC")
        # D_z must equal -eps33 * k in magnitude
        b = electro[0]
        dz = sd.raw_states[b][6]
        want = abs(mp.mpf(repr(8.8541878128e-12))) * k
        err = abs(abs(dz) - want)
        if err > mp.mpf("1e-20") * max(want, mp.mpf(1)):
            bad.append(f"D_z {mp.nstr(abs(dz), 6)} != eps33*k "
                       f"{mp.nstr(want, 6)}")
        # exactly one member of the pair is admissible per half-space
        admitted = [b for b in outgoing if b in S["admitted1"]]
        if len(admitted) != len(outgoing):
            bad.append(f"{len(admitted)}/{len(outgoing)} outgoing "
                       f"electrostatic modes admitted")
        if len(outgoing) != 1:
            bad.append(f"{len(outgoing)} outgoing electrostatic modes, "
                       f"expected exactly 1 (the other decays away)")
        if bad:
            return (FAIL, n, " | ".join(bad), "")
        return (PASS, n, f"D_z = {mp.nstr(abs(dz), 6)} = eps33*k",
                "classifying this as a gauge deleted the electrostatic channel "
                "of an elastic half-space while keeping the piezoelectric one")

    finish(c, body)


# ===========================================================================
def test_SCTA04_weak_attenuation_propagates():
    """Im(xi) != 0 must not be treated as a cusp for near-lossless waves."""
    c = make_case(
        test_id="SCTA04", category="7-energy",
        purpose="a mode whose attenuation is 3e-13 per wavelength is "
                "propagating, so PZT-5H has an incident channel at theta > 0",
        oracle="attenuation |Im xi|/|Re xi| computed from the root alone",
        oracle_independent=True,
        fixture="PZT-5H, 30 deg, omega = 2*pi*1e6 rad/s",
        expected="at least one INCOMING mode exists; attenuation of the "
                 "incident mode is below 1e-4 per wavelength",
        tolerance=1e-4,
        tol_justification="the affected case is 3e-13 per wavelength (decay "
                          "length 9.6e8 m for a 1.8 mm wave) and the weakest "
                          "genuinely damped mode is 6.7e-11, so 1e-4 sits ten "
                          "orders above one and four below the other",
        level="1")

    def body():
        mp.mp.dps = DPS
        m = fc.material_from_table("PZT-5H")
        k = W * mp.sin(mp.radians(mp.mpf(30))) / mp.mpf(4000)
        S = assemble(m, m, k, W)
        sd = S["side1"]
        incoming = [b for b, (r, _w) in enumerate(sd.roles) if r == "INCOMING"]
        n = len(sd.roots)
        if not incoming:
            return (FAIL, n, "no INCOMING mode",
                    "a strictly-positive Im(xi) test removes every incident "
                    "channel, so the medium cannot be excited at all")
        # attenuation of the weakest incoming mode must be in the band
        atts = []
        for b in incoming:
            xi = sd.roots[b]
            re_ = abs(mp.re(xi))
            if re_:
                atts.append(abs(mp.im(xi)) / re_)
        att = min(atts) if atts else mp.inf
        if att > ATTENUATION_PROPAGATING:
            return (FAIL, n, f"attenuation {mp.nstr(att, 4)} above the "
                             f"threshold", "")
        return (PASS, n, f"{len(incoming)} incoming, weakest attenuation "
                         f"{mp.nstr(att, 4)} per wavelength",
                "the strict rule declared this localised on the far side and "
                "removed the incident channel entirely")

    finish(c, body)


# ===========================================================================
def test_SCTA05_interface_continuity_and_rescaling():
    """Interface residuals recomputed independently, and rescaling invariance."""
    c = make_case(
        test_id="SCTA05", category="6-interface-conditions",
        purpose="the reconstructed interface field is continuous and does not "
                "move under a consistent mode rescaling",
        oracle="the interface state rebuilt from the modal expansion, not the "
               "linear-algebra residual",
        oracle_independent=True,
        fixture="Al6061|PZT-5H at 30 and 60 deg; PZT-5H|PZT-5H at 30 deg",
        expected="interface residual < 1e-40 and INVARIANT rescaling",
        tolerance=1e-40,
        tol_justification="square, well-conditioned (cond <= 3.2e4) systems "
                          "solve at rounding in 50-digit arithmetic; the "
                          "measured residual is 2.7e-51",
        level="1")

    def body():
        mp.mp.dps = 50
        A_al = fc.material_from_table("Al6061")
        A_pz = fc.material_from_table("PZT-5H")
        n, bad = 0, []
        for nm1, m1, nm2, m2, th in (("Al6061", A_al, "PZT-5H", A_pz, 30),
                                     ("Al6061", A_al, "PZT-5H", A_pz, 60),
                                     ("PZT-5H", A_pz, "PZT-5H", A_pz, 30)):
            k = W * mp.sin(mp.radians(mp.mpf(th))) / mp.mpf(4000)
            S = assemble(m1, m2, k, W)
            r = analyse(S)
            n += 1
            if not r["ok"]:
                bad.append(f"{nm1}|{nm2}@{th}: {r.get('reason')}")
                continue
            if not r["square"]:
                bad.append(f"{nm1}|{nm2}@{th}: not square "
                           f"({r['n_equations']}x{r['n_unknowns']})")
            if r["interface_residual_max"] > mp.mpf("1e-40"):
                bad.append(f"{nm1}|{nm2}@{th}: interface residual "
                           f"{mp.nstr(r['interface_residual_max'], 4)}")
            if r["rescale_verdict"] != "INVARIANT":
                bad.append(f"{nm1}|{nm2}@{th}: rescaling "
                           f"{r['rescale_verdict'][:30]}")
            if not r["cond_trust"]:
                bad.append(f"{nm1}|{nm2}@{th}: cond {mp.nstr(r['cond_equilibrated'], 4)} "
                           f"above COND_TRUST")
        if bad:
            return (FAIL, n, f"{len(bad)} defects", " | ".join(bad[:4]))
        return (PASS, n, "continuous and rescaling-invariant",
                "a plus sign instead of minus leaves the solver residual at "
                "1e-37 while violating the interface by 2x the incident state")

    finish(c, body)


# ===========================================================================
def test_SCTA06_coupled_is_labelled_diagnostic():
    """G4 and G5 are open; nothing coupled may be reported as validated."""
    c = make_case(
        test_id="SCTA06", category="7-energy-passivity",
        purpose="coupled amplitudes carry the DIAGNOSTIC label and no coupled "
                "energy or thermodynamic gate is claimed",
        oracle="the record's own label field, checked on every configuration",
        oracle_independent=True,
        fixture="three material pairs at 0/30/60 deg",
        expected="label == DIAGNOSTIC for every analysed configuration",
        tolerance=0.0,
        tol_justification="exact string comparison",
        level="1")

    def body():
        mp.mp.dps = DPS
        A_al = fc.material_from_table("Al6061")
        A_pz = fc.material_from_table("PZT-5H")
        n, bad = 0, []
        for nm1, m1, nm2, m2 in (("Al6061", A_al, "PZT-5H", A_pz),
                                 ("PZT-5H", A_pz, "PZT-5H", A_pz),
                                 ("Al6061", A_al, "Al6061", A_al)):
            for th in (0, 30, 60):
                k = W * mp.sin(mp.radians(mp.mpf(th))) / mp.mpf(4000)
                r = analyse(assemble(m1, m2, k, W))
                n += 1
                if r.get("label") != EVIDENCE_LABEL:
                    bad.append(f"{nm1}|{nm2}@{th}: label {r.get('label')}")
                if "R" in r or "T" in r:
                    bad.append(f"{nm1}|{nm2}@{th}: reports a validated R or T")
        if bad:
            return (FAIL, n, " | ".join(bad[:4]), "")
        return (PASS, n, f"{n} configurations labelled {EVIDENCE_LABEL}",
                "G4 and G5 remain OPEN; a converged boundary solve is not "
                "evidence that R, T and A are validated")

    finish(c, body)
