"""
Radiation conditions for half-space mode selection -- Phase 2A task 1R (P0-3).

THE FAILURE THIS REPLACES
=========================
`screen_halfspaces` assigned a mode by the sign of Im(xi) and filed every branch
with Im(xi) = 0 under `non_decaying`. That is wrong twice over:

  * Im(xi) = 0 does NOT mean inadmissible. For Al6061 at 30 deg branches 1 and 6
    have Im(xi) exactly zero and carry <J_z> = -+1.815e20 W/m^2. They are
    perfectly good propagating modes, one incoming and one outgoing. Filing them
    under "non-decaying" discarded the only information that distinguishes them.

  * Boundedness is NOT admissibility. A pure electric wire gauge has Im(xi) = 0,
    is bounded everywhere, and is a genuine null vector of the operator, but it
    carries no mechanical energy and cannot be excited by a mechanical wave.

THE PHYSICAL TEST
=================
For an isotropic elastic medium the axial energy flux of a plane wave is

    <J_z> = +-c E          with   c = sqrt(mu/rho) > 0

where E is the energy density and +-c is the axial phase velocity. So for a
lossless propagating mode,

    sign(<J_z>) = sign(axial phase velocity) = direction of energy transport

and the outgoing condition at an interface z = 0 is simply <J_z> > 0 for the
region z > 0 / <J_z> < 0 for the region z < 0, i.e. energy flows AWAY from the
interface into the region of interest.

LIMITING ABSORPTION / CONTINUATION
==================================
Checking `sign(<J_z>)` directly is only valid while the mode is genuinely
lossless. Two things can break it:

  (a) a branch sitting ON the real axis but infinitesimally below or above it,
      where sign(Im xi) and sign(<J_z>) can momentarily disagree near a cutoff
      (PZT-5H branch 3 has Re xi = +1.4e-11 and <J_z> = -1.5e10 -- opposite
      signs), and

  (b) a branch whose decay is below the arithmetic's resolving power, so that
      Im(xi) is reported as exactly 0 when it is actually +-eps.

The principled resolution is the limiting-absorption argument: analytically
continue Im(omega) -> 0+ (equivalently, add an infinitesimal loss) and take the
limit. A mode is OUTGOING if its flux is directed away from the interface in
that limit. For a lossless mode this is unambiguous and agrees with
sign(<J_z>); for a nearly-lossless mode the continuation removes the tie.

We implement the check in the order that keeps each step defensible:

  1. GAUGE          -> never a physical channel, on non-triviality grounds.
  2. Im(xi) != 0    -> OUTGOING iff it decays away from the interface.
  3. Im(xi) == 0    -> OUTGOING iff <J_z> points away from the interface,
                       subject to <J_z> being resolvable; otherwise UNRESOLVED
                       and reported as such rather than guessed.

Step 3 is where a naive implementation goes wrong in the other direction: it is
tempting to add "and if the flux is zero, treat it as outgoing anyway because it
is bounded". That is exactly the gauge-as-outgoing error.

Reproduce:
    cd code && python3 -m scattering.radiation
"""
from __future__ import annotations

import mpmath as mp

#: below this |<J_z>| the flux is not resolvable at the working precision and no
#: direction may be claimed. Matches the classifier's absolute zero test.
FLUX_RESOLVABLE = mp.mpf("1e-40")

#: Outgoing/incoming role of one mode with respect to an interface at z = 0.
OUTGOING = "OUTGOING"
INCOMING = "INCOMING"
EVANESCENT_EITHER = "EVANESCENT_NO_TRANSPORT"
GAUGE_EXCLUDED = "GAUGE_EXCLUDED"
UNRESOLVED = "UNRESOLVED_NO_RADIATION_TEST"

ROLES = (OUTGOING, INCOMING, EVANESCENT_EITHER, GAUGE_EXCLUDED, UNRESOLVED)


def classify_role(classification, side, jz):
    """Decide incoming/outgoing for one mode, in the region `side`.

    Parameters
    ----------
    classification : dict from `classify_mode`
    side : "z>0" or "z<0" -- the region the field occupies
    jz : the cycle-averaged axial flux for THIS mode (not a neighbour's)

    Returns (role, reason). Every branch returns one of ROLES plus a sentence
    stating which test decided it -- never a bare label.
    """
    if classification is None:
        return UNRESOLVED, "no classification available"
    if classification["spatial_class"] == "GAUGE":
        return (GAUGE_EXCLUDED,
                "pure electric gauge: bounded, but carries no mechanical "
                "motion, so it is not a physical scattering channel")

    im = classification["xi_im"]
    # the direction pointing AWAY from an interface at z = 0, into `side`
    away_is_plus = (side == "z>0")

    if im != 0:
        decays_toward_plus = im > 0
        decays_away = (decays_toward_plus == away_is_plus)
        if not decays_away:
            return (EVANESCENT_EITHER,
                    f"Im(xi) = {mp.nstr(im, 6)}: the mode decays AWAY from "
                    f"{side}, so it grows there and is not admissible -- it "
                    f"belongs to the other half-space")
        # decays into `side`; outgoing iff it also transports energy into it
        if jz is not None and abs(jz) > FLUX_RESOLVABLE:
            flux_away = (jz > 0) == away_is_plus
            return (OUTGOING,
                    f"decays into {side} (Im xi = {mp.nstr(im, 6)}) and "
                    f"<J_z> = {mp.nstr(jz, 6)} carries energy "
                    f"{'into' if flux_away else 'out of'} {side}")
        return (OUTGOING,
                f"decays into {side}; flux below the resolvable scale, so the "
                f"role rests on the decay direction alone")

    # --- Im(xi) == 0: lossless, so only the flux carries direction ----------
    if jz is None or abs(jz) <= FLUX_RESOLVABLE:
        # a bounded, non-decaying, zero-flux, non-gauge mode still cannot be
        # declared outgoing. Continuation does not help when the flux is zero
        # at every order: there is nothing to continue.
        return (UNRESOLVED,
                "Im(xi) = 0 and <J_z> is not resolvable: boundedness alone "
                "cannot decide the radiation condition, and a bounded mode is "
                "NOT thereby an outgoing mode")

    flux_away = (jz > 0) == away_is_plus
    return (OUTGOING if flux_away else INCOMING,
            f"lossless (Im xi = 0) in {side}: role decided by the limiting-"
            f"absorption continuation, which for a lossless mode reduces to "
            f"sign(<J_z>) = {mp.nstr(jz, 6)}; energy flows "
            f"{'away from' if flux_away else 'toward'} the interface")


def screen(modes, jzs, sides=("z>0", "z<0")):
    """Screen a full branch set for one material.

    Returns, per side, the outgoing set, the incoming set, and every exclusion
    with its stated reason. `modes` and `jzs` must be the SAME LENGTH and
    index-aligned with the roots -- the whole point of task 1R is that a value
    must never be read from a different branch, so this function refuses to
    proceed on a length mismatch rather than zipping to the shorter one.
    """
    if len(modes) != len(jzs):
        raise ValueError(
            f"branch/mode length mismatch: {len(modes)} modes vs {len(jzs)} "
            f"fluxes -- index alignment cannot be assumed")
    out = {}
    for side in sides:
        outgoing, incoming, excluded = [], [], {}
        for br, md in enumerate(modes):
            role, why = classify_role(md, side, jzs[br])
            if role == OUTGOING:
                outgoing.append(br)
            elif role == INCOMING:
                incoming.append(br)
            else:
                excluded[br] = (role, why)
        out[side] = {"outgoing": outgoing, "incoming": incoming,
                     "excluded": excluded}
    return out


def _main():
    import sys
    sys.path.insert(0, ".")
    import foundation_checks as fc
    from scattering.nullspace import null_vector
    from scattering.modal import classify_mode

    mp.mp.dps = 60
    w = 2 * mp.pi * 1e6
    for mat in ("Al6061", "PZT-5H", "CdSe"):
        m = fc.material_from_table(mat)
        for theta in (0.0, 30.0):
            k = w * mp.sin(mp.radians(theta)) / 4000.0
            modes, jzs = [], []
            for xi in [mp.mpc(z) for z in m.xi_roots(k, w)]:
                a = null_vector(fc.M_mp(m, k, xi, w), 60)["vector"]
                jz = fc.flux_z_mp(m, k, xi, w, a)
                modes.append(classify_mode(m, k, xi, w, a, jz=jz))
                jzs.append(jz)
            res = screen(modes, jzs)
            print(f"\n=== {mat} @ {theta} deg ===")
            for side in ("z>0", "z<0"):
                r = res[side]
                print(f"  {side}: outgoing={r['outgoing']} "
                      f"incoming={r['incoming']}")
                for br, (role, why) in sorted(r["excluded"].items()):
                    print(f"     br{br} excluded: {role}")


if __name__ == "__main__":
    _main()
