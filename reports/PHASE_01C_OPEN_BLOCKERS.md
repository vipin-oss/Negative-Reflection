# Phase 1C — Open Blockers

**As of commit `127eba77f397be6baaf7b24fc5b8607a903d402b` (Phase 1B state)**

This file is the readable companion to `PHASE_01C_OPEN_BLOCKERS.json`, which is
the machine-readable register that the test harness reads when deciding whether
promotion is allowed.

The distinction this file exists to protect:

> **A passing test suite is not evidence that a blocker is closed.** Several of
> the blockers below are open *precisely because* no test can close them: they
> need a derivation, a source, or a full-text review. Infrastructure that runs
> cleanly says nothing about any of them.

---

## Status summary

| id | severity | blocker | promotion-blocking |
|---|---|---|---|
| BLK-U3 | P0 | coupled amplitude normalisation / null-vector direction | yes |
| BLK-G4 | P0 | Gate G4 — LS energy identity not passed | yes |
| BLK-G5 | P0 | Gate G5 — interface dissipation and numerical balance | yes |
| BLK-G0 | P0 | Gate G0 — full-text prior art not reviewed | yes |
| BLK-G1 | P0 | Gate G1 — incompletely sourced constants | yes |
| BLK-U5 | P1 | half-factor bookkeeping between mechanical and thermal contact terms | yes |
| BLK-LS-STORAGE | P1 | **open scientific question**: vanishing cycle-averaged time derivative vs vanishing stored energy in Lord–Shulman | yes |
| BLK-DOC-01 | P1 | withdrawn τ_q causality floor still printed in `governing_equations.tex` | no |

**Promotion is refused.** `run_*.py --promotion` exits 1 while any of the seven
promotion-blocking entries is open, regardless of how many tests pass.

---

## BLK-U3 — amplitude normalisation (P0)

The null-vector *direction* of the 6 mm characteristic matrix is accurate only to
~1e-5, because the matrix has condition number ~1e22 once the displacement,
potential and temperature columns are mixed. The local energy-balance identity is
invariant under rescaling the null vector, so the 1.86e-5 residual **cannot** be a
magnitude error: it is direction error.

Confirmed independently: the same identity closes to 7.1e-17 when the null
direction is obtained by a dps=60 Cramer solve instead of the production routine.
The residual gap is explained quantitatively by the missing transverse
dissipation term scaling as (κ/ξ)² = 1.2e-7.

**What would close it:** a well-conditioned formulation — most likely
equilibrating the unknown columns before the null-space extraction — combined with
a boundary-data normalisation derived from the flux rather than from the state
vector. Loosening the tolerance would *not* close it.

**Do not** close this by declaring the residual acceptable.

---

## BLK-G4 — Gate G4 (P0)

The corrected flux functional is derived and verified against the Second Law in an
independent evaluation, but the identity does not close *through the production
code path* (see BLK-U3). Gate G4 asks for it to hold in the solver that will
actually be used.

---

## BLK-G5 — Gate G5 (P0)

No coupled scattering solver exists. Consequently:

- R, T and A have never been computed for any piezothermoelastic configuration;
- the composite power balance (EQ-ENG-08) has never been evaluated;
- interface dissipation has never been verified numerically.

Blueprint §7.4 forbids reporting R, T or A until this closes. Every reflection-
suppression figure in any future manuscript depends on it.

---

## BLK-G0 — Gate G0 (P0)

All ten literature records are metadata-only (Crossref). No full text was read.
No compatibility check with any competitor is possible, so no novelty claim can be
made and the six candidate research gaps remain *questions*, not findings.

---

## BLK-G1 — Gate G1 (P0)

`beta1`, `beta3` and the pyroelectric coefficient `p` are `ASSUMED` for PZT-5H and
BaTiO₃. The CdSe relaxation time is an order-of-magnitude value
(`VERIFIED_TO_ORDER_OF_MAGNITUDE`). Any result whose magnitude depends on the
thermal or pyroelectric channel is quantitatively unsupported.

Note the Phase 1C finding that makes this worse than it looked: until now the code
did not read `literature/material_constants.csv` at all, so an `ASSUMED` label in
the table could not influence or be checked against any computation. That is now
fixed by `test_DAT01`, which compares all 72 comparable constants across the two
sources.

---

## BLK-U5 — half-factor bookkeeping (P1)

The mechanical contact term carries a ½ from the cycle average; the thermal term
is written for a real temperature jump. This is a **convention**, not a physical
uncertainty, but applying it inconsistently produces a 50 % error in the thermal
channel alone. It must be fixed before any numerical energy balance is assembled.

---

## BLK-LS-STORAGE — the open scientific question (P1)

**This is the most important entry in this file and it is deliberately unresolved.**

Phase 1 described `τ_q q·∂_t q / κ` as a "relaxation storage" term, then as a
dissipation. Both readings are unsupported, and for a simple reason:

> For a monochromatic field, the cycle average of `q·∂_t q` is **identically
> zero**, because it is the time average of a quantity multiplied by its own time
> derivative over a full period.

Measured: `3.27e-6` against a dissipation scale of `1e34` — zero to machine
precision. A term that vanishes identically can be neither storage nor
dissipation *in the harmonic regime*.

What that does **not** establish, and what must not be inferred from it:

- it does not show that Lord–Shulman theory has no stored-energy contribution;
- it does not show that the relaxation time is thermodynamically inert;
- it says nothing about the transient regime, where the term is not zero.

The distinction that matters scientifically is between **a vanishing
cycle-averaged time derivative** and **a vanishing stored-energy term**. Those are
different statements and the second does not follow from the first. Settling it
requires a complete thermodynamic derivation of the Lord–Shulman internal energy
and of what the correct harmonic expression for its stored part is — not a
numerical test, and certainly not a tolerance adjustment.

**Preserved as a scientific question.** No test in this phase closes it, and none
should be written that appears to.

---

## BLK-DOC-01 — the document still prints a withdrawn rule (P1, not blocking)

`theory/governing_equations.tex` still contains the section "Causality floor on
τ_q", which asserts that requiring the thermal wave not to outrun `c_L` gives a
lower bound on `τ_q` and that parameters below it are rejected. That rule was
withdrawn in Phase 1B as non-physical: it conflates a constitutive relaxation time
with a signal speed.

The code was fixed. The report was fixed. **The `.tex` was not**, so the document
and the code currently disagree. It is recorded rather than quietly patched, and
`CLM-H1-02` is filed `INVALID` so the withdrawn claim remains visible.

---

## Blockers that this phase legitimately *did* reduce

For completeness, and to avoid overstating the position in either direction:

| item | before | now | what actually changed |
|---|---|---|---|
| U-11 claim audit | open | **closed** | every claim now has a status, a location and, where applicable, an executed basis |
| false PASS detection | none | **closed** | the passivity test can now fail; a negative control proves it |
| U-11-1 (C1 null space) | claimed | **claim downgraded, calculation done** | `test_EIG03` computes the null space; the claim stays `UNVERIFIED` because the *original* function still does not |
| U-11-3 (A2 scope) | overstated | **restated** | claim narrowed to its true scope; `test_DIM01` covers all 33 terms |
| code–table divergence | invisible | **detected** | `test_DAT01`, 72 constants compared cross-source |

That is progress in *knowing what is true*, which is the only kind of progress this
phase claims. It is not progress in *the physics being right*, and it does not
close a single scientific gate.
