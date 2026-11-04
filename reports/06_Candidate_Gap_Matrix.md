# 06 — Candidate Gap Matrix

**Document status:** Phase 4 deliverable
**Scope:** the ten candidate families A–J in the brief, plus four additional candidates that the
literature search surfaced.
**Labelling convention (from the brief):**
`VERIFIED GAP` = the absence is well supported by a systematic search and the required physics is
demonstrably available · `PLAUSIBLE GAP` = absence supported, but either the physics or the novelty
is uncertain · `UNVERIFIED GAP` = could not be established either way from accessible sources ·
`ALREADY INVESTIGATED` = substantially equivalent work exists.

Every row was assessed against the same nine questions demanded by the brief. The answers are
tabulated rather than narrated, because the point of this document is to make the rejections
auditable.

---

## 0. How each question was answered

| Symbol in the table | Question | How it was decided |
|---|---|---|
| **Q1** | Do the required forward/backward branches actually exist? | Computed, or read from a verified experimental paper. Not assumed. |
| **Q2** | Can physically permissible boundary mode conversion occur? | Requires an interface/edge and a mode pair with coincident wavenumber; checked via the fold. |
| **Q3** | Does the material coupling change the relevant mechanism? | Distinguishes "moves the numbers" from "changes the structure of the boundary-value problem". |
| **Q4** | Can the governing equations be derived consistently? | Checks for known pathologies (nonlocal singular limits, gradient BC completeness). |
| **Q5** | Are realistic parameter values available? | Named data source required; "not sourced" is a fail. |
| **Q6** | Can the phenomenon be demonstrated numerically? | Judged against what the solver in this project can already do. |
| **Q7** | Does a useful application follow? | Stated, or marked none. |
| **Q8** | Do substantially equivalent studies already exist? | Graded as none / adjacent / equivalent. |
| **Q9** | Label | As defined above. |

---

## 1. Summary matrix

| ID | Candidate | Q1 branches | Q2 conversion | Q3 coupling changes mechanism | Q4 derivable | Q5 parameters | Q6 numerics | Q7 application | Q8 equivalent work | **Q9 Label** |
|----|-----------|-------------|---------------|-------------------------------|--------------|---------------|-------------|----------------|--------------------|--------------|
| **A** | NR of Lamb waves in size-dependent piezoelectric plates | **Yes**, computed (fold at 1.130 MHz for PZT-5H 1.5 mm) | Yes (free edge) | **No** if "size-dependent" = piezoelectric only; **Yes** if = flexoelectric/strain-gradient | Yes for PE; **contested** for gradient | Yes (PZT-5H, LiNbO₃ well characterised) | Yes | Tunable/acoustic devices, sensing | Adjacent (M1 isotropic; backward waves in piezo plates) | **PLAUSIBLE GAP** (as PE-only it is ALREADY INVESTIGATED in substance) |
| **B** | NR in FG piezomagnetic / MEE media | **Likely**, not computed here | Yes, but needs a second medium or an edge | Yes: grading changes both dispersion and the mode shapes | Yes (pseudo-Stroh / power-series in x₂) | Yes (CoFe₂O₄, BaTiO₃, Terfenol-D) | Yes | Gradient-tailored wave control | None found combining FG + MEE + NR | **VERIFIED GAP** |
| **C** | Nonlocal strain-gradient effects on forward–backward conversion | Yes (elastic fold exists); **not** established that gradient terms preserve it | **Unclear** — extra BCs may prevent it | **Yes, qualitatively**: changes number of edge conditions | **Contested**: consistent nonlocal BCs are an open modelling issue | ℓ poorly characterised for many materials | Yes | Size-tuned devices; fundamental | None found | **VERIFIED GAP** (highest-risk, highest-value) |
| **D** | NR at imperfect interfaces of coupled waveguides | Yes | Yes, and spring stiffness is a *continuous* control | Yes: interface law changes the scattering matrix directly | Yes (standard spring model) | Yes (interfacial stiffness tabulated for bonded joints) | Yes | Bond-quality monitoring; tunable mirrors | Adjacent (imperfect joints for edge resonance; spring interfaces for dispersion) | **PLAUSIBLE GAP** |
| **E** | Thermally controlled NR in generalized thermoelastic plates | Yes in the elastic limit; **thermal coupling splits and damps** the pair | Yes | Partially: thermal terms are largely perturbations on the elastic fold | Yes (MGT/MDD well posed) | Yes (thermal relaxation times tabulated for a few metals) | Yes | Temperature-tuned devices | Adjacent (Gupta's own MGT corpus, but not for NR) | **PLAUSIBLE GAP**, low novelty |
| **F** | NR and back-refocusing in anisotropic coupled media | **Yes** (anisotropic ZGV established: Prada et al. 2009) | Yes | **Yes**: group velocity not parallel to k; fold becomes a curve, not a point | Yes, but algebraically heavy | Yes (crystal cuts) | Moderate–hard | Direction-dependent focusing | Adjacent (anisotropic ZGV; isotropic NR) | **PLAUSIBLE GAP**, high algebraic cost |
| **G** | Frequency-selective NR near ZGVs | Yes | Yes | No — this is a property of any fold | Yes | N/A | Yes | Filter design | Adjacent (Veres 2016 shows broad-angle operation) | **PLAUSIBLE GAP**, thin on its own |
| **H** | NR in flexoelectric / flexomagnetic waveguides | Yes (flexo shifts dispersion) | Unclear (same BC issue as C) | Yes | **Contested** (see C) | f ~ 1–10 µC/m; d ~ 10⁻¹⁰ N/(A·m) — order-of-magnitude only, needs sourcing | Hard (fourth-order) | Nano-electromechanical sensing | Adjacent (flexo Lamb dispersion; flexomagnetic reflection) | **PLAUSIBLE GAP** |
| **I** | Effects of porosity, viscoelasticity and grading on NR | Yes for grading; **viscoelasticity damps and lifts the degeneracy** | Yes | Yes: loss makes the criterion flux-based, not dispersion-based | Yes | Porous/viscoelastic constants available for ceramics | Yes | Damping-aware design | Adjacent (Gupta's porosity corpus; no NR) | **PLAUSIBLE GAP**, best as a robustness section |
| **J** | Passive NR and wave localization in layered composites | Yes | Yes, but "layered" shifts it from reflection to sequential refraction | Yes | Yes | Yes | Yes | Localization, energy trapping | **ALREADY INVESTIGATED in substance** for periodic/layered plates | **ALREADY INVESTIGATED** |
| **K** (new) | Correct **electromechanical** flux and energy budget for NR | Yes | N/A | Yes — necessary lemma | Yes | N/A | Yes (done here for the reduced system) | Underpins every other claim | **None found** | **VERIFIED GAP**, but a lemma, not a paper |
| **L** (new) | NR in a **polarized** (biased) piezoelectric plate: DC bias as a control | Yes | Yes | Yes — bias changes stiffening directly | Yes | Yes (bias-dependent constants measured for PZT, LiNbO₃) | Yes | Electrically tunable mirror | None found | **VERIFIED GAP** (replaces the closed avenue in 05 §3.1) |
| **M** (new) | Thickness-graded (FG) *elastic* plate: does grading preserve the ρ ≈ −1 condition? | Yes (graded ZGV exist, `10.1121/1.4983296`) | Yes | Yes — grading changes mode shapes smoothly | Yes | Yes (two-material FG systems) | Yes | **Cheapest entry point** to the whole programme | Adjacent (graded ZGV; no NR) | **VERIFIED GAP**, low glamour |
| **N** (new) | Non-Hermitian / EP-based reflection in a plate | Yes | Yes | Yes, but it is a *different mechanism* | Yes | Requires gain | Yes | Sensing | Substantial equivalent literature | **ALREADY INVESTIGATED** as an NR problem (it is an EP problem) |

---

## 2. Detailed assessment of each candidate

### A — Negative reflection of Lamb waves in size-dependent piezoelectric plates

**What the brief hoped for.** That "size-dependent piezoelectric" would supply both a new material
and a new mechanism.

**What the computation shows** (`code/check_02_flux_and_zgv.py`, `reports/figures/fig2`):

- A forward/backward fold **does** exist in PZT-5H: f_ZGV = 1.130 MHz, k_ZGV = 809 1/m at 1.5 mm.
- Piezoelectric coupling shifts the ZGV by **−35.4 %** relative to the same elastic tensor with
  coupling switched off. That is a large, real effect.
- The shift is **identical at 0.05 mm, 0.5 mm and 1.5 mm**. There is **no intrinsic size
  dependence**: the only thickness variation is the geometric f_ZGV ∝ 1/h law that every plate obeys.

**Why the "size-dependent" framing fails for pure piezoelectricity.** The piezoelectric contribution
to stiffness is e²/ε — a ratio of material constants with no length in it (here 3.25 GPa, 2.6 % of
c11). A theory with no internal length cannot produce a size effect. This is a dimensional argument,
not a numerical one, and it is decisive.

**Therefore:** as stated ("size-dependent piezoelectric"), the candidate reduces to *"run the known
isotropic M1 calculation with piezoelectric constants"*. That is **in substance already investigated**
— the backward branches in piezoelectric plates are established (Smirnov et al. 2023) and the M1 edge
scattering is established (four papers). The combination is a materials substitution.

**Salvage.** Two genuine openings survive, and they are what P1 is built on:
(i) the fold in a *coupled* plate may be a **different mode pair** than in the elastic case, so the
near-degeneracy condition must be re-established rather than assumed; (ii) a **DC bias** (candidate
L) supplies a real, continuous electrical control of the stiffness and hence of the fold — unlike the
open/short switch, which the computation shows is inert for this mode.

**Label: PLAUSIBLE GAP, but only in the form "which mode pair carries the fold in a coupled plate,
and is the near-degeneracy preserved?" — not in the form "size-dependent piezoelectric".**

### B — Negative reflection in functionally graded piezomagnetic / magneto-electro-elastic media

**Strengths.** Grading is a genuine mechanism, not a parameter: it changes the mode shapes
continuously through the thickness, which is precisely what condition (C2) depends on. Both
Dr. Gupta's own work (`10.1016/j.euromechsol.2025.105786` on FG piezomagnetic media;
`10.1038/s41598-026-46812-3` on an FG-viscoelastic-piezomagnetic trilayer) and the MEE literature
provide all the constitutive and numerical apparatus. Graded plates are known to retain ZGV points
(`10.1121/1.4983296`).

**Weaknesses.** Requires a second material or a free edge to produce conversion; a purely graded
infinite plate has no reflection. Magneto-electro-elastic constants for CoFe₂O₄/BaTiO₃ composites are
well tabulated, which is a plus, but the piezomagnetic constants are small and the magnetic
contribution to the fold will be smaller than the piezoelectric one.

**Label: VERIFIED GAP.** No paper was found combining functional grading, a coupled constitutive law
and negative reflection.

### C — Nonlocal strain-gradient effects on forward–backward guided-wave conversion

**This is the scientifically strongest candidate and the highest risk.**

**Strengths.**
- The *structural* argument is compelling and was developed in `05` §3.3: a gradient theory adds
  natural boundary conditions (higher-order tractions) and possibly essential conditions on normal
  derivatives. The M1 mechanism works because **two** modes suffice to satisfy the **two** classical
  traction conditions. A gradient theory needs more. Either an edge boundary layer appears, with a
  quantifiable penalty on ρ, or the formulation degenerates.
- Either outcome is a publishable, novel statement — this is the rare case where the negative result
  is as informative as the positive one.
- It is squarely in Dr. Gupta's methodological history (nonlocal dispersion, energy-ratio
  consistency), while requiring him to do something new (edge scattering with higher-order BCs).

**Risks, stated plainly.**
- Consistent boundary conditions for Eringen-type *integral* nonlocal elasticity are genuinely
  unsettled in the literature. A referee can attack the formulation itself.
- The physically accessible regime requires ℓ/h to be non-negligible, i.e. micrometre plates and
  high frequencies (§6 of `05`). The demonstration will be numerical, not experimental.
- There is a real chance the extra conditions are satisfiable by the existing evanescent set with
  negligible effect, in which case the paper's conclusion is "no significant change".

**Label: VERIFIED GAP** (the absence of edge-scattering calculations in gradient plates is well
supported). Novelty: very high. Feasibility: medium.

### D — Negative reflection at imperfect interfaces of coupled elastic waveguides

**Strengths.** The spring-interface model is standard, well posed and directly supplies a *control
parameter*: the ratio of interfacial stiffness to the elastic impedance. Dr. Gupta has published on
imperfect interfaces repeatedly. Conversion at a *joint* between two plates (rather than a free
edge) is a natural extension of the M2 thickness-step geometry.

**Weaknesses.** "Imperfect interface" is a crowded field; the novelty claim must be specifically
about the *negative-reflection efficiency* as a function of stiffness, not about the interface model.
Risk that the result is a monotone degradation curve.

**Label: PLAUSIBLE GAP.**

### E — Thermally controlled negative reflection in generalized thermoelastic plates

**Assessment.** Thermal coupling in a generalized thermoelastic plate is, for the *elastic* fold,
essentially a perturbation plus damping. MGT and memory-dependent theories are Dr. Gupta's home
territory, but the physics here is thin: temperature will shift and broaden the fold, and the
interesting question (does the ρ ≈ −1 condition survive thermal damping?) is largely answered by
"no, but only weakly".

**Label: PLAUSIBLE GAP with low novelty.** Best used as a robustness section, not as a paper.

### F — Negative reflection and back-refocusing in anisotropic coupled media

**Strengths.** In an anisotropic plate the ZGV locus is a curve, not a point, and the group velocity
is not parallel to the wave vector. This makes "negative reflection" genuinely more subtle — the sign
test must be applied to in-plane components — and it is a real gap.

**Weaknesses.** The algebra and the numerical cost are both much higher, and the marginal physical
insight over the isotropic case is uncertain. Realistic for a strong group; heavy for a single
investigator.

**Label: PLAUSIBLE GAP.**

### G — Frequency-selective negative reflection near zero-group-velocity modes

**Assessment.** As a standalone paper this is a bandwidth study of a known effect. The dispersion
curvature does predict the bandwidth (this follows from the flux identity), but the M1 papers already
report the operating bandwidth experimentally. Weak on novelty.

**Label: PLAUSIBLE GAP, thin.** Folded into P2 as a section.

### H — Negative reflection in flexoelectric or flexomagnetic waveguides

**Assessment.** Same structural boundary-condition issue as C, plus much weaker and less well
characterised material constants. Flexoelectric coefficients for a given material span orders of
magnitude in the literature. Best regarded as a variant of C, not a separate problem.

**Label: PLAUSIBLE GAP**, subsumed by C.

### I — Effects of porosity, viscoelasticity and material grading on negative reflection

**Assessment.** Porosity and viscoelasticity both damp; viscoelasticity additionally lifts the
degeneracy so ρ never reaches −1. This is important as a *reality check* (every real device is
lossy) but it is not a headline result. Dr. Gupta's porosity corpus makes it cheap for him.

**Label: PLAUSIBLE GAP**, best as an ablation.

### J — Passive negative reflection and wave localization in layered composite structures

**Assessment.** For a *periodic* or *layered* plate, sequential refraction, band formation and
localization are heavily studied in the elastic metamaterial literature (`10.1088/1361-665X/ad0393`),
and the "cloaking/trapping" line (`10.1038/s41598-021-03146-6`) already does exactly the
"arrange regions to localize waves" idea. A layered-composite NR paper would be read as a
re-derivation.

**Label: ALREADY INVESTIGATED in substance.**

### K — Correct electromechanical flux and energy budget for NR (new)

A lemma, not a paper: but every other claim depends on it, and §2 of `05` shows that getting it
wrong produces exactly the factor-2π class of error that occurred in this very project.

**Label: VERIFIED GAP.** To be included as a mandatory preliminary section of whichever problem is
chosen, not as a standalone submission.

### L — DC-biased piezoelectric plate: bias as the control parameter (new)

**Why it matters.** The computation in `05` §3.1 closes the "open vs short" avenue for the S₁/S₂ᵦ
fold (−0.0006 %). A DC bias, by contrast, changes the *effective stiffness itself* through the
converse piezoelectric effect and is a continuous, experimentally accessible knob. Nothing in the
searched literature applies a bias to a negative-reflection problem.

**Label: VERIFIED GAP.** This is the credible replacement for the avenue that was closed, and it is
carried into P1.

### M — Thickness-graded elastic plate: does grading preserve ρ ≈ −1? (new)

**Why it matters.** It is the **cheapest possible entry** into the whole programme: no new
constitutive physics, a well-posed eigenvalue problem, and a direct test of the near-degeneracy
condition (C2). It is also a natural extension of Dr. Gupta's grading experience.

**Label: VERIFIED GAP**, low glamour, high probability of a clean result. Carried into P2.

### N — Non-Hermitian / exceptional-point reflection (new)

**Label: ALREADY INVESTIGATED** as an NR problem. It is an exceptional-point problem; the mechanism
is coalescence, not forward↔backward conversion, and it requires gain. Excluded from the five
problems, retained in the taxonomy for completeness.

---

## 3. What the matrix implies

Three things fall out of the table, and they drive Phase 5:

1. **The "size-dependent piezoelectric" framing (A) does not survive on its own.** The computation
   shows piezoelectric coupling is size-independent, and the backward-branch existence in
   piezoelectric plates is already published. A proposal must therefore be about **which pair carries
   the fold and whether the near-degeneracy survives** (A, re-scoped), or about **bias** (L), or about
   **a genuine internal length** (C, H, M).

2. **Candidate C is the only one where the coupling changes the *structure* of the problem** rather
   than its coefficients. That is the difference between a mechanism paper and a parameter paper, and
   it is why C leads.

3. **Several candidates (E, G, I) are ablation sections, not papers.** Including them as separate
   problems would inflate the list and weaken it. The five problems in Phase 5 are chosen so that
   each has a distinct *mechanism-level* claim, not a distinct material.

---

## 4. Explicit statement of what is *not* claimed

- I do **not** claim that candidates B, C, D, F, L or M are unexplored in the absolute sense. I claim
  that no paper was found in the accessible literature that does them, and I grade them accordingly.
  Scopus and Web of Science were not reachable (§5 of `01`), so a systematic search of the
  subscription indexes has **not** been performed and must be performed by the investigator before
  committing.
- I do **not** claim that negative reflection exists in any specific coupled material. It is
  established for isotropic plates (four papers). For coupled plates, **the existence of the fold is
  computed here, but its usefulness for negative reflection is not** — that requires the edge
  calculation and the near-degeneracy check, both of which are future work.
- I do **not** claim that the gradient extension will yield an observable effect. That is the
  central open risk.

---

## 5. Reproducibility of this phase

| Claim | Where it can be checked |
|---|---|
| "no paper combines coupled media with edge NR" | search log, `03` §6; query families 1×{coupled}, 5, 7, 8 |
| PZT-5H fold at 1.130 MHz / 809 1/m | `data/check_02_flux_and_zgv.json`, key `C_table` |
| coupling shifts ZGV by −35.4 %, identically at three thicknesses | same file, keys `C_coupling_pct_*` |
| open/short shifts it by −0.0006 % | same file, keys `C_open_short_pct_*` |
| mode potential is antisymmetric | `reports/figures/fig3b`, `code/plate_spectral.py::eigen_full` |
| piezoelectric stiffening is size-independent (e²/ε) | `05` §3.1, arithmetic from the PZT-5H constants |
