# 04 — Research Gap Evidence

**Document status:** Phase 2 deliverable
**Rule applied throughout:** an absence of results in a search index is recorded as an *absence of
evidence*, not as *evidence of absence*. Every gap below is graded accordingly, and the grading is
carried forward into `06_Candidate_Gap_Matrix.md`.

---

## 1. What "gap" is allowed to mean here

Negative reflection of elastic guided waves (mechanism M1, see `03_Negative_Reflection_Taxonomy.md`)
is a solved problem in the isotropic case. The question is therefore not "does negative reflection
exist" but:

> **Which specific, checkable claim about M1 is currently unsupported by any published result, and
> would change how the phenomenon is understood or used if it turned out either way?**

Four candidate claim-types survived screening. Each is stated below in falsifiable form, with
(a) what is established, (b) what is not, (c) how I searched, (d) the strength of the gap.

---

## GAP-1 — No modal scattering calculation at a free edge of a *coupled* plate

**Claim tested.** *"The efficiency of forward↔backward conversion at a free edge in a
piezoelectric (or piezomagnetic, or magneto-electro-elastic) plate is governed by the same
near-degeneracy condition as in the isotropic case, and the electrical/magnetic surface condition
does not change it."*

**(a) Established.** M1 at a free edge of an isotropic plate, including the reflection-coefficient
matrix as a function of frequency and angle: Germano 2012 (`10.1103/PhysRevB.85.012102`),
Gérardin 2016 (`10.1121/1.4959024`), Veres 2016 (`10.1103/PhysRevB.93.174304`),
Gérardin 2019 (`10.1038/s41598-019-38480-3`). Backward branches with ZGV points in piezoelectric
plates: Smirnov et al. 2023 (`10.3390/s23020648`) — experimentally established for YX-LiNbO₃.

**(b) Not established.** No verified source computes edge reflection coefficients for a coupled
plate. The piezoelectric literature establishes that the *dispersion relation* has backward branches
and ZGV points; it never computes what happens at an edge. The negative-reflection literature
establishes the edge calculation but only for a passive isotropic solid.

**(c) How searched.** Query family 1 crossed with {piezoelectric, piezomagnetic,
magneto-electro-elastic, flexoelectric}; query family 5 in full; forward citations of the four M1
papers screened for any coupled-medium follow-up. Result: zero matches combining an edge-scattering
calculation with a coupled constitutive law.

**(d) Gap strength: PLAUSIBLE, bordering on VERIFIED.** The absence is well supported — the two
literatures are large enough that a bridging paper would very likely have surfaced. But I cannot
exclude a paper in a venue not indexed by OpenAlex or a non-English-language journal.

**Why it might be scientifically interesting rather than incremental.** Two non-obvious things can
happen and neither is obvious in advance:

- The pair (S₁, S₂ᵦ) is defined by *elastic* branch repulsion. Piezoelectric coupling shifts and
  reorders branches, so in a coupled plate the ZGV pair may not even be the same pair of modes.
  Whether a *usable* forward/backward pair still exists near the same k is an open question.
- The short-circuit condition φ = 0 adds a constraint that a traction-free edge does not have. The
  edge problem in a piezoelectric plate is therefore a *different* boundary-value problem, not the
  same one with modified coefficients. Whether the near-degeneracy condition (C2 of the taxonomy)
  survives that extra constraint is genuinely unknown.

**Evidence already generated in this project, and why it matters for this gap.** A limited
premise test was run (`code/check_02_flux_and_zgv.py`, `reports/figures/fig2`, `fig3`):

- The S₁/S₂ᵦ ZGV of a PZT-5H plate poled through the thickness lies at **f = 1.12975 MHz,
  k = 808.95 1/m** for a 1.5 mm plate, versus **1.89038 MHz, 1061.0 1/m** for the isotropic
  Duralumin reference — a **−35.4 %** frequency shift. A usable fold does exist, but is not where
  the elastic intuition puts it.
- The shift is **identical at 0.05 mm, 0.5 mm and 1.5 mm** thickness. Piezoelectric coupling alone
  produces **no intrinsic size dependence** — only the geometric f_ZGV·h = const scaling. This is a
  direct, computed refutation of the naive "piezoelectric ⇒ size-dependent" premise.
- The electrical surface condition (open vs short) shifts this particular ZGV by
  **−0.0006 %**. Per-branch sensitivity is 10³–10⁵ times larger on other branches at the same
  wavenumber (`reports/figures/fig3a`), so the insensitivity is specific to this mode, not a
  coding artefact. The recovered potential of the mode is **antisymmetric through the thickness**
  with φ = 0 at the mid-plane (`reports/figures/fig3b`), which is consistent with both faces being
  effectively charge-free.

**What this does to the gap.** It narrows it usefully. The interesting question is no longer
"does piezoelectricity matter" (it does, by 35 %) but "**which** mode pair carries the fold, and does
the extra electrical boundary condition at the edge destroy the ρ ≈ −1 conversion that the elastic
case enjoys". That is a sharper and more defensible question, and it is the one carried into P1.

---

## GAP-2 — No treatment of higher-order (strain-gradient / nonlocal) boundary conditions in edge scattering

**Claim tested.** *"A strain-gradient or nonlocal internal length changes the forward↔backward
conversion efficiency at a free edge, and does so through the boundary conditions rather than
through the dispersion relation alone."*

**(a) Established.** Size-dependent dispersion for coupled plates, extensively: e.g. Chen, Guo & Pan
2017 (`10.1016/j.jsv.2017.04.001`), the gradient electro-mechanics Stroh formalism for Lamb waves
(`10.1063/5.0196873`), and the review `10.1080/19475411.2022.2091058`. Dr. Gupta's own nonlocal
dispersion work (`10.1080/15376494.2023.2291804`, `10.1108/hff-10-2023-0615`, and others) is in this
class.

**(b) Not established.** In gradient theories the traction boundary condition is augmented by a
**higher-order traction** (a double force / moment traction, depending on the variant), and the
variational principle yields additional essential conditions on the normal derivative of the
displacement. In a nonlocal (integral) theory the boundary conditions become **nonlocal** — the
stress at an interior point depends on the strain in a neighbourhood, and the correct consistent
edge condition is itself an open modelling question in the literature. **Every verified size-dependent
plate paper I found solves a dispersion problem (infinite plate, no edge) or a bending/vibration
problem. None solves a semi-infinite plate with a free edge and a modal scattering expansion.**

**(c) How searched.** Query family 7 in full; supplemented by searches on {strain gradient,
free edge, boundary conditions, plate}, {nonlocal, edge condition, wave scattering}, and the
reference lists of the gradient-electro-mechanics and flexomagnetic papers.

**(d) Gap strength: VERIFIED for the specific absence; PLAUSIBLE for the consequence.**
The absence of edge-scattering calculations in gradient plates is well supported. That such a
calculation would *change the answer* is a physical hypothesis, not an established fact — and it
carries a real risk (documented in `08`), because the extra boundary conditions may turn out to be
satisfiable in a trivial way that leaves the reflection matrix essentially unchanged.

**Why it is not merely "add a length scale".** The point is structural. The M1 mechanism works
because the boundary condition can be satisfied by *two* modes. A gradient theory needs *more*
conditions at the edge (a higher-order traction, or an integral condition). Two modes cannot
generally satisfy three or four conditions. So either (i) additional evanescent modes must supply
the extra boundary layer, in which case the conversion coefficient is necessarily degraded by an
amount that depends on the ratio of the internal length to the wavelength, or (ii) the theory is
degenerate at the edge and the problem is ill-posed. **Deciding which of these happens is a real
result either way**, and it is a result Dr. Gupta is unusually well placed to obtain, since he has
already spent years on nonlocal dispersion and on energy-ratio consistency at interfaces.

This is the strongest and most original of the four gaps.

---

## GAP-3 — Imperfect (spring-type) *edges* have never been considered as NR control parameters

**Claim tested.** *"A finite tangential/normal stiffness at the reflecting boundary tunes the
forward↔backward conversion continuously between the free-edge and clamped limits, and negative
reflection survives over a finite stiffness range."*

**(a) Established.** Spring-type imperfect interfaces are standard in the layered-media literature:
membrane-type MEE interfaces (`10.1016/j.compstruct.2022.115661`), nonlocal MEE nanoplates with
imperfect contacts (`10.1080/15376494.2024.2328753`), and Dr. Gupta's own work on layer
imperfections in cylindrical piezoelectric composites (`10.1016/j.compstruct.2025.118997`) and on an
imperfect fluid–porous interface (`10.47974/jios-1608`). Imperfect *joints* of plates have been
studied for edge resonance, but not for NR.

**(b) Not established.** No study treats the *reflecting edge itself* as an imperfect boundary with
finite stiffness in the negative-reflection context. The M1 papers all use a perfectly free edge.

**(c) How searched.** Query family 8 in full; combined with {edge, reflection coefficient, Lamb,
backward}. No matches.

**(d) Gap strength: PLAUSIBLE.** The absence is real, but the scientific payoff is uncertain in
advance: physically, a spring edge is a continuous interpolation between free and clamped, and the
interesting question is whether NR *degrades monotonically* (boring) or *survives up to a threshold
and then switches* (interesting). The latter is plausible — the conversion relies on a delicate
degeneracy, so a small stiffness may be enough to break it — but this is a hypothesis. It is also
the candidate most exposed to the criticism "this is a parameter study of a known effect".

---

## GAP-4 — No frequency-selective NR criterion stated in terms of the flux zero

**Claim tested.** *"The bandwidth over which NR is efficient is set by the curvature of the
dispersion fold, and can be predicted from the flux identity v_e = P/(2⟨T⟩) without solving the
edge problem."*

**(a) Established.** The ZGV point and the angular/frequency range of efficient conversion are known
empirically from the four M1 papers, and Veres 2016 shows broad-angle operation.

**(b) Not established.** I found no paper that derives the NR bandwidth from the *dispersion
curvature* and states an explicit predictive criterion. The M1 papers compute the reflection matrix
numerically and read the bandwidth off the result.

**(c) How searched.** Forward and backward citation screening of the M1 corpus; searches on
{negative reflection, bandwidth, criterion}; {zero group velocity, curvature, reflection}.

**(d) Gap strength: PLAUSIBLE but LOW VALUE on its own.** It is a useful analytical contribution
and it is nearly free once the dispersion solver exists (it was in fact validated here: the flux
zero reproduces the ZGV to 10⁻⁶, `code/check_02_flux_and_zgv.py`). But as a standalone paper it is
thin. It belongs as a *section* of a larger study, not as the study.

---

## 2. Gaps considered and rejected

For completeness, and because the brief asks for critical evaluation rather than a list of
opportunities, these were examined and set aside.

**"Negative reflection in a non-Hermitian / PT-symmetric plate."** Attractive, but the reflection
mechanism there is exceptional-point coalescence, not forward↔backward M1 conversion. That is a
*different* mechanism, it requires gain (which is not passive), and PT-symmetric elasticity is a
well-developed field with its own community (`10.1016/j.jmps.2021.104325`). The novelty claim would
have to be about EPs, not about NR. Rejected as a *negative-reflection* problem; retained only as a
possible comparison case.

**"Time-modulated plate producing temporal negative reflection."** Category error risk is high, the
phenomenon is temporal, and the term "negative reflection" does not apply in its spatial sense.
Rejected.

**"Negative reflection in a phononic-crystal or metasurface-coated plate."** Mechanism M3. This is a
design/fabrication problem, and it sits far from Dr. Gupta's analytical strengths. Rejected on fit.

**"Negative reflection with viscoelastic loss."** Interesting and cheap, but the physics is a
straightforward attenuation of both branches; the main effect is that the degeneracy is lifted and
ρ never reaches −1. Rejected as insufficient novelty on its own, though it is retained as a
**necessary robustness check** inside the chosen problems (any real device is lossy).

**"Negative reflection in a rotating plate (Coriolis)."** Rotation breaks reciprocity, which is
interesting, but the effect on the elastic fold is second order and the parameter regime needed for
an observable shift is not physically accessible in a laboratory plate. Rejected on parameter
feasibility.

---

## 3. Consolidated gap table

| ID | Gap | Strength | Novelty if true | Risk | Carried into |
|----|-----|----------|-----------------|------|--------------|
| GAP-1 | Edge scattering in a *coupled* plate | Plausible→verified absence; premise quantified here | High (first modal scattering matrix for a coupled plate) | Medium (may reduce to isotropic result with shifted parameters) | P1, P2 |
| GAP-2 | Higher-order BCs in gradient edge scattering | **Verified absence** | Very high (structural result either way) | High (problem may be ill-posed or trivially unchanged) | P3 (lead), P1 |
| GAP-3 | Imperfect edge as an NR control parameter | Plausible | Medium | Medium (may be a monotone parameter study) | P4 |
| GAP-4 | Bandwidth criterion from dispersion curvature | Plausible, low value | Low–medium | Low | Section of P2 |

Plus one gap that emerged from the *numerical* work rather than the literature search and is
recorded here for honesty:

| ID | Gap | Strength | Note |
|----|-----|----------|------|
| GAP-5 | No published check that "backward ⇒ negative reflection" holds in a coupled medium | **Verified absence** | The identity v_e = dω/dk holds for the lossless elastic case (verified to 10⁻⁷ here). In a piezoelectric plate the correct flux is electromechanical and the identity must be re-derived. No paper found that does this for a plate. Small result, but it is a necessary lemma for anything else, and getting it wrong invalidates a whole paper. | 

---

## 4. What would falsify the whole programme

Stated plainly, so the committee knows the failure mode:

> If a coupled plate's forward/backward pair near its ZGV has **dissimilar stress–displacement
> fields** — i.e. if coupling destroys the near-degeneracy condition (C2) of the taxonomy — then
> the free edge cannot satisfy the boundary condition with the pair alone, ρ will not approach −1,
> and there will be **no usable negative reflection in that material**, regardless of how the
> boundary conditions are formulated.

This is a real possibility and it is tested cheaply, before any heavy derivation, by computing the
mode shapes of the pair and their overlap. That test is scheduled as the first task of the proposed
work (see `10_RESEARCH_PROPOSAL_FOR_APPROVAL.md`, §11, Step 1). If it fails, the programme stops
there — and that outcome should be reported, because it would be a genuine negative result about
coupled-media negative reflection.
