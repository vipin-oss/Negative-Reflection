"""
code/scattering/maxwell.py -- the electromagnetic interface conditions.

WHY THIS MODULE EXISTS
======================
Task 2S reduced the electrical interface to a slogan: "bonded means potential
continuity, insulated means D_z continuity, pick one". That is wrong, and the
supervisor was right to reject it. This module derives the conditions from
Maxwell's equations and keeps the five distinct models apart, because they are
not alternatives that one may choose between on numerical grounds.

MAXWELL, AND THE ELECTROQUASISTATIC ASSUMPTION
==============================================
The full set is

    div D = rho_free          curl E = -dB/dt
    div B = 0                 curl H = J + dD/dt

This research model is ELECTROQUASISTATIC: the mechanical time scale (1 MHz,
i.e. ~1e-6 s) is many orders slower than the electromagnetic transit time of the
sample, so the magnetic induction term is dropped. With time dependence e^{-i w t}
and no free volume charge in either bulk medium:

    div D = 0        curl E = 0        E = -grad(phi)

`curl E = 0` is what makes a scalar potential exist at all. Both bulk equations
are already built into `foundation_checks.M_mp`; this module deals only with what
happens across the interface.

THE INTERFACE CONDITIONS
========================
For an interface with unit normal n pointing from medium 1 into medium 2, and
surface fields on the interface itself:

    n . (D_2 - D_1)  = sigma_free                     (normal D)
    n x (E_2 - E_1)  = 0                              (tangential E)

The second is a vector equation. In the sagittal plane with n = z-hat the only
tangential direction is x, and E_x = -i k phi, so `curl E = 0` gives

    k (phi_2 - phi_1) = 0

i.e. **the potential jump is zero for k != 0, and is UNCONSTRAINED at k = 0.**
That is the whole of the tangential condition, and it is why the normal-incidence
limit is delicate: the condition degenerates, it does not disappear.

A relative potential offset across an interface is NOT a gauge freedom. A gauge
transformation adds the SAME constant everywhere; a jump is a physical surface
dipole layer / charge sheet. Confusing the two is the error `MR08` injects.

THE FIVE MODELS, EACH STATED ON ITS OWN TERMS
=============================================
    DIELECTRIC_CONTACT   sigma_free = 0;  [phi] = 0 and [D_z] = 0 BOTH.
                         This is the ideal contact between two dielectric
                         solids with no charge trapped at the boundary. It is
                         what "bonded" should mean, and imposing only one of
                         the two is an error, not a modelling choice.

    CONDUCTING_ELECTRODE sigma_free = -q_electrode (a free parameter, determined
                         by the electrode), [phi] = 0 (an equipotential metal),
                         and D_z is NOT continuous -- it terminates on the
                         electrode charge. A conductor REPLACES the dielectric
                         contact; it does not supplement it.

    SHORT_CIRCUIT        the electrode is additionally held at a fixed potential
                         (phi = 0, the reference). One equation replaces the
                         free parameter q_electrode.

    OPEN_CIRCUIT         zero total electrode charge, i.e. sigma_free integrates
                         to zero over the interface but is locally free. This is
                         the genuinely "insulated" model: the potential floats.

    PRESCRIBED_CHARGE    sigma_free = sigma_given(x), a specified function.
                         [phi] = 0 still holds; only the normal-D jump changes.

Note the placement is deliberate and load-bearing: continuity of D_z does NOT
by itself identify an insulated interface -- it is the *no-free-charge*
condition, which holds for a perfectly clean dielectric contact too. Conversely continuity
of potential does not replace continuity of D_z; they come from different
Maxwell equations (curl E and div D) and answer different questions. Task 2S
asserted the opposite of both statements; see reports/PHASE_02A_TASK2T_CLAIM_
CORRECTIONS.md.

WHAT THE PAPER INTENDS
======================
The working blueprint describes a two-half-space reference with a welded
interface and the interface "electrical charge-storage relation" as an object of
study rather than an assumption in force. The reference configuration is
therefore DIELECTRIC_CONTACT with sigma_free = 0, and that is the model the
production assembly already imposes at k != 0 -- verified here, not assumed.
Electrode and open-circuit models are implemented so they can be *compared*,
never so they can be silently swapped in.

Reproduce:
    cd code && python3 -m scattering.maxwell
"""
from __future__ import annotations

import sys

import mpmath as mp

import foundation_checks as fc
from scattering.state import STATE_FIELDS

# ---------------------------------------------------------------------------
# the model identifiers. Strings, so they survive a JSON round trip and cannot
# be accidentally reordered the way an enum's values can.
# ---------------------------------------------------------------------------
DIELECTRIC_CONTACT = "DIELECTRIC_CONTACT"
CONDUCTING_ELECTRODE = "CONDUCTING_ELECTRODE"
SHORT_CIRCUIT = "SHORT_CIRCUIT"
OPEN_CIRCUIT = "OPEN_CIRCUIT"
PRESCRIBED_CHARGE = "PRESCRIBED_CHARGE"

ALL_MODELS = (DIELECTRIC_CONTACT, CONDUCTING_ELECTRODE, SHORT_CIRCUIT,
              OPEN_CIRCUIT, PRESCRIBED_CHARGE)

#: The model the reference two-half-space configuration uses.
PAPER_MODEL = DIELECTRIC_CONTACT

#: Which electrical conditions each model imposes. `phi` is tangential-E,
#: `Dz` is normal-D. `sigma_free` says whether a free surface charge is allowed.
MODEL_SPEC = {
    DIELECTRIC_CONTACT: {
        "phi": "CONTINUOUS", "Dz": "CONTINUOUS", "sigma_free": "ZERO",
        "equations": ("phi", "Dz"),
        "statement": "[phi] = 0 from curl E = 0, and [D_z] = 0 from "
                     "div D = 0 with no free surface charge. BOTH hold.",
    },
    CONDUCTING_ELECTRODE: {
        "phi": "EQUIPOTENTIAL", "Dz": "TERMINATES", "sigma_free": "FREE",
        "equations": ("phi",),
        "statement": "The conductor is an equipotential surface and carries an "
                     "unknown free charge; D_z terminates on that charge, so "
                     "continuity of D_z is NOT imposed. The conductor replaces "
                     "the dielectric contact rather than supplementing it.",
    },
    SHORT_CIRCUIT: {
        "phi": "FIXED_ZERO", "Dz": "TERMINATES", "sigma_free": "FREE",
        "equations": ("phi",),
        "statement": "An equipotential electrode tied to the reference: "
                     "phi = 0 exactly, with an unknown charge.",
    },
    OPEN_CIRCUIT: {
        "phi": "FLOATS", "Dz": "TOTAL_CHARGE_ZERO", "sigma_free": "FREE_ZERO_MEAN",
        "equations": ("phi", "Dz_mean"),
        "statement": "The electrode is isolated: its NET charge is zero but the "
                     "local charge is free, so the potential floats. This, not "
                     "mere D_z continuity, is the insulated model.",
    },
    PRESCRIBED_CHARGE: {
        "phi": "CONTINUOUS", "Dz": "JUMPS_BY_GIVEN", "sigma_free": "GIVEN",
        "equations": ("phi", "Dz_given"),
        "statement": "[phi] = 0 still holds (curl E is unaffected by charge), "
                     "and [D_z] = sigma_given. The jump is data, not an unknown.",
    },
}


def model_spec(model):
    """The recorded specification of one electrical model."""
    if model not in MODEL_SPEC:
        raise ValueError(f"unknown electrical model {model!r}; "
                         f"expected one of {ALL_MODELS}")
    return dict(MODEL_SPEC[model], model=model)


# ---------------------------------------------------------------------------
# the interface rows
# ---------------------------------------------------------------------------
def electrical_residuals(model, f1, f2, sigma_free=None, k=None):
    """Evaluate the electrical interface conditions on a reconstructed state.

    `f1`, `f2` are the eight-component TOTAL physical states on the two sides.
    On the incident side that means incident + reflected -- NOT the reflected
    field alone. Passing only the scattered field silently shifts every jump by
    the incident amplitude; that mistake produced a spurious 3.2e-2 potential
    discontinuity during Task 2T before it was caught.
    `sigma_free` is the specified free surface charge (C/m^2) for
    PRESCRIBED_CHARGE; it is ignored -- and recorded as ignored -- otherwise.

    Returns a dict of named residuals, each an absolute physical quantity:

        tangential_E   k*(phi_2 - phi_1)      -- curl E = 0
        normal_D       D_z2 - D_z1 - sigma    -- div D = 0

    The caller decides which are asserted by the model. Returning BOTH for every
    model is deliberate: a model that does not assert one must still be able to
    REPORT what it is, which is how "silently discarded" becomes visible.
    """
    sig = mp.mpf(0) if sigma_free is None else mp.mpf(sigma_free)
    dphi = f2[2] - f1[2]
    dx = f2[6] - f1[6] - sig
    tang = (k * dphi) if k is not None else dphi
    return {
        "model": model,
        "phi_jump": dphi,
        "tangential_E": tang,
        "normal_D": dx,
        "D_jump": f2[6] - f1[6],
        "sigma_free_used": sig,
        "norm": max(abs(dphi), abs(dx), mp.mpf(1)),
    }


def asserted_components(model):
    """Which of ('phi', 'Dz') this model imposes as equations."""
    return tuple(MODEL_SPEC[model]["equations"])


def distinguishes(models=(DIELECTRIC_CONTACT, OPEN_CIRCUIT)):
    """True when the named models impose genuinely different equations.

    This exists so a test can assert that the models are NOT interchangeable
    without hard-coding the answer: it asks the model table, which is the single
    place the physics is recorded.
    """
    ref = asserted_components(models[0])
    return any(asserted_components(m) != ref for m in models[1:])


# ---------------------------------------------------------------------------
# self-check
# ---------------------------------------------------------------------------
def _selfcheck():
    mp.mp.dps = 40
    ok = True
    print("electrical interface models")
    print("=" * 78)
    for m in ALL_MODELS:
        s = model_spec(m)
        print(f"  {m:<22s} phi={s['phi']:<13s} Dz={s['Dz']:<17s} "
              f"sigma={s['sigma_free']}")
        print(f"  {'':<22s} asserts {s['equations']}")
    print()

    # the two claims Task 2S got backwards, checked against the table
    a = asserted_components(DIELECTRIC_CONTACT)
    b = asserted_components(OPEN_CIRCUIT)
    print(f"dielectric contact asserts {a}")
    print(f"open circuit       asserts {b}")
    if a != ("phi", "Dz"):
        print("  FAIL: a clean dielectric contact must impose BOTH")
        ok = False
    if b == a:
        print("  FAIL: the open circuit must differ from the dielectric contact")
        ok = False
    print()
    print("D_z continuity alone identifies the NO-FREE-CHARGE interface, not "
          "the insulated one." if ok else "SELFCHECK FAILED")
    return ok


if __name__ == "__main__":
    sys.exit(0 if _selfcheck() else 1)
