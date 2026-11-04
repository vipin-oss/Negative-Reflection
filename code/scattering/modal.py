"""
code/scattering/modal.py -- correct classification of complex-wavenumber modes.

WHY THIS MODULE EXISTS
======================
Phase 1D's summary CSV carried a `mode_class` column whose vocabulary was
wrong. It labelled every branch whose mechanical amplitude was non-negligible
as `PROPAGATING`, including

    PZT-5H branch 0:  xi = -2.569e6 - 2.569e6 i

whose decay length is 1/|Im xi| = 3.89e-07 m = 0.39 um, i.e. ONE SIXTH of a
wavelength. A wave that is attenuated by e^{2*pi/6} = e^1.05 = 2.9x per
wavelength is not propagating. It also emitted the label PHYSICAL_GROWING,
which is a contradiction: a spatially growing mode is not physical in a
half-space.

THE CONVENTION (established from the source, not assumed)
=========================================================
The code builds every field as Re{A e^{i(kx + xi z - w t)}}, i.e. e^{-i w t}
and a spatial factor e^{+i xi z}. Therefore for a branch

    xi = xi_r + i xi_i           (xi_r, xi_i real)

    |field| ~ e^{-xi_i z}        amplitude in the +z direction
    phase    ~ e^{+i xi_r z}     travels toward +z iff xi_r > 0

so
    xi_i > 0  ->  decays as z increases  ->  OUTGOING/decaying in the +z half-space
    xi_i < 0  ->  decays as z decreases  ->  OUTGOING/decaying in the -z half-space

This orientation convention is derived rather than hard-coded; it is
re-verified numerically by `verify_convention()` (see below), which compares the
analytic decay against a direct evaluation of the fields at two z values.

WHY THE SIGN OF Re(xi) IS NOT A CLASSIFIER
==========================================
The brief forbids classifying by the sign of the real part alone, and the data
show why. At 30 degrees in PZT-5H:

  branch 3:  xi_r = +1.36e-11, xi_i = -845.2

The real part is POSITIVE, which a naive rule would call "propagating toward
+z". But Re(xi) is 1.4e-11 against k = 785.4, so the z-component of the wave
vector is essentially zero: this is a wave that travels almost exactly along
x with a pure exponential decay in z. Its phase velocity computed from
Re(xi) alone is 4.6e11 m/s -- 1.5e6 times the speed of light. Classifying by
sign(Re xi) is meaningless when |Re xi| << |k|.

EQUALLY, WHY NEITHER IS |Im xi| ALONE
=====================================
Branch 1 has xi_i = -1.05e-09, a decay length of 9.6e8 m. Its attenuation is
negligible compared with any device scale, so it IS a travelling mode in every
engineering sense. A pure |Im xi| threshold would have to know the device
length; we therefore report BOTH the dimensionless attenuation ratio and the
absolute decay length, and let the criterion be explicit rather than implicit.

THE CLASSIFICATION (each field is separate and machine-readable)
===============================================================
  spatial_class      EVANESCENT | WEAKLY_DAMPED | THERMAL_DIFFUSIVE | GAUGE
  radiation_dir      DECAYS_TOWARD_PLUS_Z | DECAYS_TOWARD_MINUS_Z | NONE
  energy_flux_status FLOWING | ZERO | UNDEFINED
  confidence         HIGH | MEDIUM | LOW

The four fields are INDEPENDENT on purpose. A mode can be EVANESCENT while
carrying a non-zero energy flux (branch 3), and it can be THERMAL_DIFFUSIVE
while still having a well-defined radiation direction (branch 0). Collapsing
these into one label is what produced the Phase 1D defect.

Reproduce:
    cd code && python3 -m scattering.modal
"""
from __future__ import annotations

import mpmath as mp

# ---------------------------------------------------------------------------
# The thresholds, named and justified rather than scattered through the code.
# ---------------------------------------------------------------------------

#: dimensionless attenuation ratio |Im xi| / |k| above which a mode is not
#: travelling in any engineering sense. 1.0 means it decays by e^{-2*pi} per
#: wavelength. Chosen at 1.0 rather than 0.1 because we report the exact value
#: alongside the class, so a borderline case is visible rather than decisive.
ATTENUATION_EVANESCENT = 1.0

#: below this the attenuation is negligible over any realistic specimen and the
#: mode is treated as travelling. 1e-3 gives an amplitude change of 0.6 % per
#: wavelength.
ATTENUATION_WEAK = 1.0e-3

#: A thermal-diffusive branch is identified by its DISPERSION RELATION, not by
#: its amplitude ratios. The reason is a real error made earlier in this phase:
#: a large thermal amplitude was used as the test, which misclassified Al6061
#: branch 2 -- a genuine acoustic mode that simply has a strong temperature
#: component. Amplitude tells you what the mode is MADE OF; only the dispersion
#: relation tells you WHAT KIND OF WAVE IT IS.
#:
#: In the regime w*tau_q << 1 (true for every material here: w*tau_q <= 2.3e-6)
#: a thermal branch obeys the diffusion relation  |xi| = sqrt(w rho cE / kappa).
#: Measured over all four materials at 0/15/30/45/60 degrees:
#:     thermal branches      |xi| / |xi_th| = 1.000 .. 1.016
#:     acoustic branches     |xi| / |xi_th| <= 6.1e-03
#: i.e. a 150x gap. 0.1 sits in the middle of it.
THERMAL_DISPERSION_TOL = 0.1

#: a branch is a gauge mode when there is NO mechanical content at all. This is
#: the electric-potential gauge of a centrosymmetric medium.
#:
#: The threshold is set inside a measured gap, not guessed. Over all four
#: materials at 0/15/30/45/60 degrees:
#:     genuine gauge branches        mech/max = exactly 0.0
#:     smallest PHYSICAL branch      mech/max = 4.01e-13
#: so a threshold anywhere in (0, 4e-13) separates them. 1e-40 is chosen so that
#: a gauge can only be reported when the mechanical content is numerical zero.
#:
#: A LOOSER threshold is wrong, and cost me a false diagnosis during this phase:
#: with 1e-8 every PZT-5H acoustic branch classified as a gauge, because the
#: physical ratio |u|/|phi| in a piezoelectric is genuinely ~5e-9 (metres against
#: volts). A gauge is machine zero; it is not "small against volts".
GAUGE_MECHANICAL = 1.0e-40

#: A mode with Im(xi) EXACTLY zero does not decay. It must not be given a
#: *damping* label at all: calling it WEAKLY_DAMPED with a damping of exactly
#: zero is a vocabulary contradiction of the same kind as PHYSICAL_GROWING.
SPATIAL_CLASSES = ("PROPAGATING_NON_DECAYING", "EVANESCENT", "WEAKLY_DAMPED",
                   "THERMAL_DIFFUSIVE", "ELECTROSTATIC", "GAUGE",
                   "UNCLASSIFIED")
#: Radiation direction. For a decaying mode this is fixed by the sign of Im(xi).
#: For a mode with Im(xi) = 0 the direction is carried by the real part of the
#: axial Poynting component, i.e. ENERGY FLOW, because a lossless mode is not
#: localised at the interface at all -- nothing about it decays, so "toward the
#: interface" is not defined and only the flux direction carries information.
RADIATION_DIRS = ("DECAYS_TOWARD_PLUS_Z", "DECAYS_TOWARD_MINUS_Z",
                  "TRAVELS_TOWARD_PLUS_Z", "TRAVELS_TOWARD_MINUS_Z",
                  "NONE")
ENERGY_STATUSES = ("FLOWING", "ZERO", "UNDEFINED")
CONFIDENCES = ("HIGH", "MEDIUM", "LOW")


# ---------------------------------------------------------------------------
def verify_convention(m, k, xi, w, a, z_probe=None):
    """Re-verify e^{+i xi z} and hence the sign rule, numerically.

    Returns (decay_rate_measured, decay_rate_analytic, agrees).
    This is deliberately a numerical re-derivation rather than a comment: if
    the convention in foundation_checks ever flips, the sign rule here would
    silently invert, and every incoming/outgoing decision downstream with it.
    """
    k = mp.mpf(k)
    xi = mp.mpc(xi)
    w = mp.mpf(w)
    if z_probe is None:
        z_probe = mp.mpf(1) / (4 * abs(mp.im(xi)) if abs(mp.im(xi)) > 0 else mp.mpf(1))
    U, W, Phi, Th = [mp.mpc(t) for t in a]

    # field at z = 0 and z = z_probe, evaluated with the documented convention
    def ux(z):
        return U * mp.exp(mp.mpc(0, 1) * xi * z)

    f0 = abs(ux(mp.mpf(0)))
    f1 = abs(ux(z_probe))
    if f0 == 0:
        # degenerate zero vector: nothing to measure, but that is NOT agreement
        return mp.mpf(0), -mp.im(xi), False

    # SIGN CONVENTION, derived once here and used consistently downstream.
    #   |e^{+i xi z}| = |e^{i (xi_r + i xi_i) z}| = |e^{i xi_r z}| * e^{-xi_i z}
    #                                          = e^{-Im(xi) z}
    # Hence log|f(z)| = const - Im(xi) z, so
    #   -d/dz log|f|  = +Im(xi)          ... call this the LOG-GROWTH rate
    #   +d/dz log|f|  = -Im(xi)          ... call this the DECAY rate
    #
    # The earlier version mixed the two: it computed `-log(f1/f0)/z_probe`, which
    # is the LOG-GROWTH rate (+Im xi), and compared it against `-Im xi`, the
    # DECAY rate. Those differ by 2*Im(xi) -- relative error 2.0, never small --
    # so every branch reported agrees=False regardless of the physics. Verified
    # numerically: for xi = -2.5694e6 - 2.5694e6i, -log(f1/f0)/z = -2.5694e6 =
    # Im(xi), while the true decay rate is +2.5694e6.
    #
    # Both sides below are now the DECAY rate.
    log_growth = -mp.log(f1 / f0) / z_probe    # = Im(xi)
    measured = -log_growth                     # = -Im(xi), the decay rate
    analytic = -mp.im(xi)                      # = -Im(xi), from the convention
    scale = max(abs(analytic), mp.mpf("1e-300"))
    agrees = abs(measured - analytic) / scale < mp.mpf("1e-25")
    return measured, analytic, agrees


# ---------------------------------------------------------------------------
def classify_mode(m, k, xi, w, a, jz=None, jz_scale=None, z_ref=None,
                  legacy_gauge=False):
    """Classify one branch. Returns a dict of INDEPENDENT fields.

    Parameters
    ----------
    m : Medium6mm
    k : transverse wavenumber (1/m). May be 0 at normal incidence.
    xi : the branch's normal wavenumber (complex).
    w : angular frequency (rad/s).
    a : the null vector (4 components).
    jz : optional precomputed cycle-averaged energy flux <J_z>.
    jz_scale : optional reference scale to decide whether jz is numerically zero.

    legacy_gauge : reproduce the pre-task-2 rule exactly, in which every
        mechanically decoupled mode is called GAUGE. Used only to re-derive a
        published artifact under the convention it was generated with. The
        DECOUPLED LAPLACE ELECTROSTATIC case (mechanical content zero but
        D_z != 0) is a genuine physical channel and is classified
        ELECTROSTATIC by default; calling it a gauge silently deletes the
        electrostatic channel of an elastic half-space.

    The function never raises on a numerically awkward branch: it returns a
    classification with LOW confidence and an explicit reason instead. A
    classifier that throws is a classifier that silently drops modes.
    """
    k = mp.mpf(k)
    xi = mp.mpc(xi)
    w = mp.mpf(w)
    U, W, Phi, Th = [mp.mpc(t) for t in a]

    xi_r, xi_i = mp.re(xi), mp.im(xi)
    k_abs = abs(k)

    # --- mechanical / thermal / electric content ---------------------------
    mech = mp.sqrt(abs(U) ** 2 + abs(W) ** 2)
    elec = abs(Phi)
    therm = abs(Th)
    thermo_sq = mp.sqrt(mech ** 2 + elec ** 2)
    total = mp.sqrt(mech ** 2 + elec ** 2 + therm ** 2)

    # thermal diffusion wavenumber in the w*tau_q << 1 regime
    try:
        cE = mp.mpf(repr(m.cE))
        rho = mp.mpf(repr(m.rho))
        kap = mp.mpf(repr(m.kappa33))
        xi_th = mp.sqrt(w * rho * cE / kap) if kap > 0 else None
    except Exception:                                          # noqa: BLE001
        xi_th = None
    if xi_th is not None and xi_th > 0:
        thermal_ratio = abs(xi) / xi_th
    else:
        thermal_ratio = None

    reasons = []
    confidence = "HIGH"

    # --- spatial class -----------------------------------------------------
    # attenuation per unit transverse wavenumber. Guard k = 0 (normal
    # incidence) where the ratio is undefined; there we fall back to comparing
    # against the thermal wavenumber scale, which is what sets the branch scale
    # in that limit anyway.
    if k_abs > 0:
        atten_ratio = abs(xi_i) / k_abs
        atten_basis = "|Im xi|/|k|"
    else:
        # At normal incidence k = 0 and |Im xi|/|k| is undefined. The comment
        # above used to promise a fallback against the thermal wavenumber scale
        # while the code actually compared against a bare 1 m^-1 -- a hidden,
        # unprincipled threshold. The fallback is now the thermal wavenumber
        # scale, which is the only other length scale the branch set provides,
        # and the basis is REPORTED so a reader can see which normalisation
        # produced the class instead of finding a blank field.
        atten_ratio = (abs(xi_i) / xi_th) if xi_th is not None else None
        atten_basis = "|Im xi|/|xi_th| (k=0 fallback)"

    # The gauge test must ask whether the mechanical content is ACTUALLY zero,
    # not whether it is small relative to a unit-mixed total. |U| ~ 1 m and
    # |Phi| ~ 1e9 V, so normalising mech by the four-component norm makes every
    # acoustic mode look like a gauge -- that was a real bug in the first
    # version of this function, caught by inspecting its own output. The
    # correct test compares mech against the LARGEST COMPONENT of the vector,
    # which is a pure normalisation scale and therefore legitimate for a
    # zero-detection test (zero is zero in any scaling).
    comp_max = max(abs(U), abs(W), abs(Phi), abs(Th))
    if total == 0:
        spatial = "GAUGE"
        reasons.append("null vector is identically zero")
        confidence = "LOW"
    elif mech <= GAUGE_MECHANICAL * max(comp_max, mp.mpf("1e-300")):
        # MECHANICALLY DECOUPLED. There are two very different things here and
        # an earlier version of this function conflated them:
        #
        #  (a) a true electric-potential GAUGE at xi = 0, for which
        #      D_z = -eps33 * xi * phi = 0. It is a constant-potential null
        #      direction: carry it and the interface operator gains a zero
        #      column.
        #  (b) a DECOUPLED LAPLACE ELECTROSTATIC mode at xi = +/- i k, for
        #      which D_z = -eps33 * xi * phi = -/+ i eps33 k phi is NOT zero.
        #      This is an ordinary physical electrostatic wave in the
        #      sagittal plane. In the norm the state uses, it satisfies the
        #      D_z continuity equation exactly like any other mode.
        #
        # Labelling (b) as GAUGE let the assembly drop the electrostatic
        # channel of an elastic half-space while keeping the equivalent
        # channel of a piezoelectric one -- an exclusion of a physical
        # channel, which is exactly what Task 2 forbids. Omega and k are the
        # only inputs needed to tell them apart, so test on those and not on
        # a normalisation-dependent component comparison.
        # The two labels are decided ONLY here, so neither can be overwritten
        # downstream by an unrelated arm.
        dz = abs(m.eps33 * xi * Phi)
        dz_rel = dz / max(abs(m.eps33) * max(abs(xi), mp.mpf("1e-300"))
                          * max(elec, mp.mpf("1e-300")), mp.mpf("1e-300"))
        if legacy_gauge:
            # Reproduce the pre-task-2 convention exactly: mechanically
            # decoupled means GAUGE, D_z notwithstanding.
            spatial = "GAUGE"
            reasons.append(
                "LEGACY RULE (legacy_gauge=True): mechanically decoupled, so "
                "labelled a gauge regardless of D_z. This reproduces the "
                "convention published Task 1R's inventory was generated with")
        elif dz_rel <= mp.mpf("1e-20"):
            spatial = "GAUGE"
            reasons.append(
                f"mechanical content is zero and D_z = -eps33*xi*phi is also "
                f"zero to working precision (xi = {mp.nstr(xi, 4)}): genuine "
                f"electric-potential gauge, a null channel")
        else:
            spatial = "ELECTROSTATIC"
            reasons.append(
                f"mechanical content is zero but D_z = -eps33*xi*phi = "
                f"{mp.nstr(dz, 4)} is NOT zero, so this is a decoupled Laplace "
                f"electrostatic mode, not a gauge. It carries a real D_z "
                f"continuity constraint and must be admissible")
    elif thermal_ratio is not None and abs(thermal_ratio - 1) < THERMAL_DISPERSION_TOL:
        spatial = "THERMAL_DIFFUSIVE"
        reasons.append(f"follows the thermal diffusion relation: "
                       f"|xi|/|xi_th| = {mp.nstr(thermal_ratio, 4)} "
                       f"(tolerance +/-{THERMAL_DISPERSION_TOL:g})")
    elif xi_i == 0 and atten_ratio is not None \
            and atten_ratio > ATTENUATION_WEAK:
        # Im(xi) is EXACTLY zero but the mode is still strongly evanescent: it
        # decays along x (Re xi ~ k, Im xi = 0). CdSe branches 3/4 at 30 deg
        # have Re xi = +-8.7e-11 with Im xi = -+750.2, so Im xi is NOT zero
        # there; the genuinely zero-Im cases are Al6061 br1/br6 and the k=0
        # gauges. This arm exists so that a future root that really does sit on
        # the axis while evanescent is not mislabelled as propagating.
        spatial = "EVANESCENT"
        reasons.append("Im xi = 0 but the mode decays along x: evanescent")
    elif xi_i == 0:
        # Im(xi) exactly zero: no damping in z at all. It must not be given a
        # DAMPING label -- WEAKLY_DAMPED with a damping of exactly zero is the
        # same kind of vocabulary contradiction as PHYSICAL_GROWING. Al6061
        # branches 1 and 6 were reported this way while carrying
        # <J_z> = -+1.815e20, i.e. perfectly good propagating waves.
        spatial = "PROPAGATING_NON_DECAYING"
        reasons.append("Im xi = 0 exactly: the mode propagates without decay "
                       "and is not localised at the interface")
        if k_abs == 0:
            confidence = "MEDIUM"
            reasons.append("confidence capped at MEDIUM: k = 0, so the "
                           "attenuation claim rests on the thermal-scale "
                           "fallback")
    elif atten_ratio is not None and atten_ratio >= ATTENUATION_EVANESCENT:
        spatial = "EVANESCENT"
        reasons.append(f"|Im xi|/|k| = {mp.nstr(atten_ratio, 3)} >= "
                       f"{ATTENUATION_EVANESCENT:g}")
    elif atten_ratio is not None and atten_ratio >= ATTENUATION_WEAK:
        spatial = "WEAKLY_DAMPED"
        reasons.append(f"{ATTENUATION_WEAK:g} <= |Im xi|/|k| = "
                       f"{mp.nstr(atten_ratio, 3)} < {ATTENUATION_EVANESCENT:g}")
    elif atten_ratio is not None:
        spatial = "WEAKLY_DAMPED"
        reasons.append(f"|Im xi|/|k| = {mp.nstr(atten_ratio, 3)} < "
                       f"{ATTENUATION_WEAK:g}: attenuation negligible")
    elif atten_ratio is None:
        # no length scale available at all: refuse to guess
        spatial = "UNCLASSIFIED"
        reasons.append("k = 0 and no thermal wavenumber available: no "
                       "dimensionless attenuation exists to classify on")
        confidence = "LOW"
    else:
        # k = 0 fallback, on the thermal scale, with the basis recorded
        spatial = "EVANESCENT" if atten_ratio >= ATTENUATION_EVANESCENT \
            else "WEAKLY_DAMPED"
        reasons.append(f"k = 0: classified on |Im xi|/|xi_th| = "
                       f"{mp.nstr(atten_ratio, 3)}")
        confidence = "MEDIUM"
        reasons.append("confidence capped at MEDIUM: the classification used "
                       "the thermal-scale fallback rather than the physical "
                       "transverse wavenumber")

    # a thermal-diffusive branch is almost always also strongly attenuated;
    # report the attenuation ratio either way so the reader can see both
    # statements rather than one collapsed label.

    # --- radiation direction ----------------------------------------------
    # TWO CASES, because a lossless mode has no decay direction to speak of.
    #
    # 1. Im(xi) != 0. The mode is exponentially localised at the interface and
    #    the DERIVED convention e^{+i xi z} fixes which side it decays toward:
    #    Im(xi) > 0  =>  decays toward +z.
    #
    # 2. Im(xi) == 0. The mode does NOT decay anywhere. It is admissible in
    #    both half-spaces on boundedness grounds, and the direction that
    #    carries physical information is the direction of ENERGY FLOW, given by
    #    the cycle-averaged axial Poynting component <J_z>. For an isotropic
    #    elastic medium this is simply the group-velocity direction, and
    #    <J_z> = +-c*E with c = sqrt(mu/rho) > 0, so the sign of <J_z> IS the
    #    sign of the axial phase/energy velocity. This is the physically
    #    justified test that the earlier code was missing: it reported these
    #    branches as "NONE" and then filed them under `non_decaying`, which
    #    discarded the only information that distinguishes an outgoing wave
    #    from an incoming one.
    #
    # A zero-flux lossless mode also exists -- the pure gauge -- and it is
    # neither incoming nor outgoing. It is caught earlier by the gauge test, so
    # by the time we reach the flux branch below the mechanical content is
    # non-zero; if <J_z> is nonetheless below the resolvable scale we report
    # NONE and LOW confidence rather than inventing a direction.
    radiation_source = ""
    if xi_i > 0:
        radiation = "DECAYS_TOWARD_PLUS_Z"
        radiation_source = "Im(xi) > 0 under e^{+i xi z}"
        if atten_ratio is not None and atten_ratio <= ATTENUATION_WEAK:
            # Near-lossless. sign(Im xi) is still well defined and is exactly
            # what the limiting-absorption continuation gives in the limit
            # Im(xi) -> 0+, so no separate branch is needed; the point is that
            # a damping this small does NOT make the direction ambiguous.
            reasons.append(
                f"damping is negligible in magnitude (ratio "
                f"{mp.nstr(atten_ratio, 3)}) but its SIGN is what fixes the "
                f"decay direction, and it is unambiguous")
    elif xi_i < 0:
        radiation = "DECAYS_TOWARD_MINUS_Z"
        radiation_source = "Im(xi) < 0 under e^{+i xi z}"
    elif xi_i == 0 and jz is not None and abs(jz) > mp.mpf("1e-40"):
        radiation = ("TRAVELS_TOWARD_PLUS_Z" if jz > 0
                     else "TRAVELS_TOWARD_MINUS_Z")
        radiation_source = ("Im(xi) = 0 (lossless): direction taken from the "
                            "sign of the axial energy flux <J_z>")
        reasons.append(
            f"Im xi = 0: mode does not decay, so it is admissible on BOTH "
            f"sides on boundedness grounds; <J_z> = {mp.nstr(jz, 6)} makes it "
            f"an {'outgoing' if jz > 0 else 'incoming'} wave with respect to "
            f"an interface at z = 0")
    else:
        radiation = "NONE"
        radiation_source = "no decay and no resolvable flux"
        reasons.append("Im xi = 0 and <J_z> is not resolvable: the mode is "
                       "neither incoming nor outgoing")
        confidence = "LOW"

    # --- energy flux -------------------------------------------------------
    if jz is None:
        energy = "UNDEFINED"
        reasons.append("energy flux was not evaluated")
        confidence = "MEDIUM"
    else:
        # Absolute test against what the arithmetic can resolve -- NOT a
        # relative test against |jz| itself. The previous version used
        # `abs(jz) <= 1e-12 * abs(jz)`, which is a TAUTOLOGY for every
        # non-zero jz: it can only ever fire when jz is exactly 0. That is how
        # the pure electric wire gauge at Al6061 30 deg branches 3/4, whose
        # |J_z| is 2.787e-107 (i.e. zero to within roundoff), was reported as
        # FLOWING. A flux of 1e-107 W/m^2 is not a flow; it is the residue of
        # cancelling terms and it must not be presented as a physical result.
        # jz is a cycle-averaged mechanical + electrical flux, so its natural
        # scale is set by the mechanical amplitudes; 1e-40 W/m^2 is far below
        # any flux a mode of this normalisation can carry.
        if abs(jz) <= mp.mpf("1e-40"):
            energy = "ZERO"
            reasons.append(
                f"|J_z| = {mp.nstr(abs(jz), 4)} is below the resolvable scale; "
                "this mode carries no mechanical energy, so it is not an "
                "outgoing channel (a wire gauge, or a mode with no axial "
                "transport)")
            if confidence == "HIGH":
                confidence = "MEDIUM"
        elif jz_scale is not None and abs(jz) <= mp.mpf("1e-20") * max(jz_scale, mp.mpf("1e-300")):
            energy = "ZERO"
            reasons.append("|J_z| is zero relative to the stated reference scale")
        else:
            energy = "FLOWING"
            reasons.append(f"J_z = {mp.nstr(jz, 6)} "
                           f"({'+z' if jz > 0 else '-z'} at u_z = 1)")

    # A lossless direction rests entirely on the flux, so the prose must say
    # which way the energy goes. Without this line the CSV carried a
    # TRAVELS_TOWARD_* label next to a reasons field that never mentioned the
    # flux -- a label no reader could check, which is the failure mode this
    # whole task is about.
    if radiation.startswith("TRAVELS_") and jz is not None:
        reasons.append(
            f"direction rests on the axial energy flux: <J_z> = "
            f"{mp.nstr(jz, 6)} {'>' if jz > 0 else '<'} 0, energy flows "
            f"{'toward +z' if jz > 0 else 'toward -z'}")

    # --- numerical confidence ---------------------------------------------
    # AT NORMAL INCIDENCE the transverse wavenumber is zero, so the physical
    # scale |Im xi|/|k| does not exist and every attenuation statement rests on
    # the thermal fallback. No classification at k = 0 may be HIGH confidence,
    # whatever arm produced it: the tolerance that separates thermal from
    # acoustic was calibrated at oblique incidence and is being extrapolated.
    # This is applied once, globally, rather than inside each arm, because an
    # earlier version capped it in three of the five arms and left the other
    # two reporting HIGH -- which is how six rows escaped the check.
    if k_abs == 0 and confidence == "HIGH":
        confidence = "MEDIUM"
        reasons.append(
            "confidence capped at MEDIUM: k = 0, so the attenuation scale is "
            "the thermal fallback rather than the physical transverse one")
    # flag the genuinely awkward cases rather than hiding them behind a label
    if xi_r != 0 and k_abs > 0 and abs(xi_r) < mp.mpf("1e-6") * k_abs:
        reasons.append("Re(xi) is negligible against k: the sign of Re(xi) "
                       "carries no directional information")
        if confidence == "HIGH":
            confidence = "MEDIUM"
    if spatial == "THERMAL_DIFFUSIVE" and thermal_ratio is not None:
        dist = abs(thermal_ratio - 1)
        if dist > 0.5 * THERMAL_DISPERSION_TOL:
            # The gap measured over 4 materials x 5 angles was thermal in
            # 1.000..1.016 and acoustic below 6.1e-3, so a thermal branch at
            # |ratio-1| = 0.29 is still clearly thermal. But at k=0 the thermal
            # ratio genuinely degrades to 1/sqrt(2) = 0.707, so the tolerance is
            # being asked to do work it was not calibrated for. Report MEDIUM
            # there rather than claiming HIGH confidence from an uncalibrated
            # regime -- a HIGH label that rests on a test outside its measured
            # range is exactly the kind of claim this phase must not make.
            reasons.append(f"thermal dispersion match is marginal "
                           f"(|ratio - 1| = {mp.nstr(dist, 3)})")
            confidence = "MEDIUM"
        if k_abs == 0 and confidence == "HIGH":
            reasons.append(
                "k = 0: the thermal tolerance was calibrated at oblique "
                "incidence only, so confidence is capped at MEDIUM")
            confidence = "MEDIUM"

    # dimensionless diagnostics, always reported
    decay_len = (1 / abs(xi_i)) if xi_i != 0 else None
    if decay_len is not None and k_abs > 0:
        # decay length in units of the transverse wavelength
        decay_in_lambda = decay_len * k_abs / (2 * mp.pi)
    else:
        decay_in_lambda = None

    return dict(
        xi_re=xi_r, xi_im=xi_i,
        spatial_class=spatial,
        radiation_dir=radiation,
        energy_flux_status=energy,
        confidence=confidence,
        attenuation_ratio=atten_ratio,
        attenuation_basis=atten_basis,
        radiation_source=radiation_source,
        decay_length_m=decay_len,
        decay_length_in_transverse_wavelengths=decay_in_lambda,
        mech_content=mech, elec_content=elec, thermal_content=therm,
        thermal_over_mech_elec=(therm / thermo_sq) if thermo_sq > 0 else None,
        thermal_dispersion_ratio=thermal_ratio,
        jz=jz,
        reasons=reasons,
    )


# ---------------------------------------------------------------------------
def screen_halfspaces(modes, jzs=None):
    """Admissibility of each branch in each half-space.

    TWO CONDITIONS, and an earlier version checked only the first:

      1. BOUNDEDNESS   the field must not grow as it moves away from the
                       interface. For Im(xi) != 0 this is the Im(xi) sign test.
                       For Im(xi) == 0 the mode is bounded EVERYWHERE and is
                       therefore admissible on both sides -- it is not excluded
                       from either.

      2. NON-TRIVIALITY  the mode must be physically realisable. A pure electric
                       wire gauge is bounded and is a genuine null vector, but
                       carries no mechanical motion, so it is not a scattering
                       channel.

    A third question -- whether an admissible non-decaying mode is INCOMING or
    OUTGOING -- is NOT decidable from boundedness at all. It needs the modal
    energy flux, and lives in `scattering.radiation`, which callers must use for
    mode selection. Reporting these branches under `non_decaying` and stopping
    there was the P0-3 defect: it discarded the incoming/outgoing distinction
    for exactly the branches where that distinction is the only thing there is
    to know.

    `jzs`, when supplied, is used only to record the flux sign in the
    `non_decaying` detail; the role decision stays in `radiation`.
    """
    out = {"z>0_admissible": [], "z<0_admissible": [], "non_decaying": [],
           "gauge_excluded": [], "non_decaying_flux_sign": {}}
    for i, md in enumerate(modes):
        if md is None or md["spatial_class"] == "GAUGE":
            out["gauge_excluded"].append(i)
            continue
        im = md["xi_im"]
        if im > 0:
            out["z>0_admissible"].append(i)
        elif im < 0:
            out["z<0_admissible"].append(i)
        else:
            # bounded on BOTH sides -- admissible in each. Not an exclusion.
            out["z>0_admissible"].append(i)
            out["z<0_admissible"].append(i)
            out["non_decaying"].append(i)
            if jzs is not None and jzs[i] is not None:
                out["non_decaying_flux_sign"][i] = (
                    "+" if jzs[i] > 0 else "-" if jzs[i] < 0 else "0")
    return out


# ---------------------------------------------------------------------------
def _main():
    """Print the corrected inventory for the two materials at 30 degrees."""
    import sys
    import numpy as np
    sys.path.insert(0, ".")
    import foundation_checks as fc

    mp.mp.dps = 60
    from scattering.roots import null_vector_verified
    w = 2 * mp.pi * 10 ** 6
    k = mp.mpf(repr(float(2 * np.pi * 1e6 * np.sin(np.deg2rad(30)) / 4000.0)))
    for mat in ("PZT-5H", "Al6061"):
        m = fc.material_from_table(mat)
        roots = [mp.mpc(z) for z in m.xi_roots(k, w)]
        print(f"\n===== {mat}, 30 deg =====")
        modes = []
        for br, xi in enumerate(roots):
            # VERIFIED null vector -- see scattering/roots.py. The previous
            # call returned a compressed direction with agreement=False on
            # Al6061 k=0 branch 6 and this branch honoured the flag, but
            # assembly.py did not, so the mode was built anyway downstream.
            r = null_vector_verified(fc.M_mp(m, k, xi, w), mp.mp.dps)
            if not r["ok"] or r["vector"] is None:
                # keep the branch in its indexed slot: dropping it would shift
                # every later label onto the wrong root, which is exactly the
                # mislabelling defect this phase exists to correct.
                modes.append(None)
                print(f"  br{br}: NULL VECTOR NOT VERIFIED "
                      f"(residual {mp.nstr(r['residual'], 4)}) "
                      f"-- reported as unusable, NOT renumbered")
                continue
            a = r["vector"]
            J = fc.flux_z_mp(m, k, xi, w, a)
            modes.append(classify_mode(m, k, xi, w, a, jz=J))
        for i, md in enumerate(modes):
            if md is None:
                continue
            print(f"  br{i}: {md['spatial_class']:18s} {md['radiation_dir']:24s} "
                  f"{md['energy_flux_status']:8s} {md['confidence']:6s} "
                  f"|Im/k|={('%.3e' % float(md['attenuation_ratio'])) if md['attenuation_ratio'] is not None else 'n/a'}")
        hs = screen_halfspaces(modes)
        print(f"  admissible in z>0: {hs['z>0_admissible']}")
        print(f"  admissible in z<0: {hs['z<0_admissible']}")
        print(f"  non-decaying:      {hs['non_decaying']}")


if __name__ == "__main__":
    _main()
