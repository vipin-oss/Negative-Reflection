# 10 — Research Proposal for Approval

> ## STATUS: AWAITING INDEPENDENT RESEARCH REVIEW — NOT APPROVED

**Submitted for review, not for execution.** Nothing in this document has been approved, funded,
scheduled, or started. No manuscript has been written. No experimental work is proposed as
already-performed. All numerical results cited here were produced by the scripts listed in §9 and
are reproducible from them.

**Submitted by:** research-discovery exercise for Dr. Vipin Gupta, Department of Mathematics,
Gurugram University, Gurugram, India (ORCID 0000-0002-5238-1963)
**Date:** 2026-11-04
**Basis:** `01`–`09` of this report set

---

## 1. Recommendation

**Primary: P1 — *Bias-tuned negative reflection in a polarized piezoelectric plate.***
**Companion: P4 — *Stiffness-controlled negative reflection at an imperfectly bonded joint.***

The two are recommended together because they share one solver, one validation framework and one
manuscript-grade dataset, and because they answer the same practical question from two directions:
P1 asks whether the *material* can tune the effect, P4 asks whether the *structure* can. Together they
constitute one coherent 18-month programme; separately, each is a single paper.

**P3 (strain-gradient)** is **not** recommended for immediate execution. It is named as the
**escalation path**: it should be scoped as the follow-on project, entered only after P1's Gate 1
result establishes that the near-degeneracy survives coupling in the classical case — because if it
does not, P3's question is moot.

**P2** is held as the **fallback**: if P1's bias-constant sourcing fails, P2 is the low-risk
substitute with no structural conditions.

---

## 2. Why P1 first — the case, stated honestly

Three arguments, in decreasing order of strength.

**(i) It starts from work already done and already validated.** This is unusual and it is the strongest
argument. Before this proposal was written, the following were established in this project
(`code/check_02_flux_and_zgv.py`, `data/check_02_flux_and_zgv.json`):

- The forward/backward fold in a coupled PZT-5H plate is at **f_ZGV = 1.12975 MHz, k_ZGV = 809.0 m⁻¹**
  for 1.5 mm thickness — located, not hypothesised.
- Coupling shifts it by **−35.4 %** relative to the elastically identical uncoupled plate, identically
  at 0.05, 0.5 and 1.5 mm.
- The open/short electrical switch shifts it by **−0.0006 %**, while other branches at the same
  wavenumber are 10³–10⁵ times more sensitive — so the insensitivity is a property of this mode.
- The mode's potential is **antisymmetric through the thickness**, φ = 0 at the mid-plane.
- The flux identity v_e = P/(2⟨T⟩) = dω/dk holds to **5.3 × 10⁻⁷** for the lossless plate, giving a
  tracking-free ZGV detector.

A cold start would spend 6 months getting to this point. P1 spends them on the actual research.

**(ii) It closes a direction that looked promising, and opens the right one.** The obvious idea —
"switch the electrical boundary condition and tune the negative reflection" — is **dead**, killed by
computation before a paper was written about it. The replacement — DC bias, which changes the
*stiffness* rather than the *boundary condition* — is physically different and is untouched in the
literature. Reporting the dead end is part of the contribution, and it is a genuinely useful service
to the community.

**(iii) Its claim is a mechanism claim, not a materials claim.** The question is *whether
electromechanical coupling preserves the near-degeneracy that negative reflection depends on* — and
the answer is not implied by anything published, because the negative-reflection literature has never
touched a coupled medium and the piezoelectric backward-wave literature has never computed an edge.

**What is honestly weak:** originality 7.0/10 (not 9), because the phenomenon and its interpretation
are unchanged — only its carrier moves. And the bias-dependent constitutive data are not yet verified
to exist. A reviewer who values only novelty should prefer P3, and §7 says so.

---

## 3. P1 — Full design

### 3.1 Title
**Electromechanical coupling and the near-degeneracy of forward/backward Lamb modes: a modal
scattering theory for the DC-biased piezoelectric plate edge.**

### 3.2 One-sentence statement
*Determine whether the forward/backward mode pair that produces negative reflection at a free edge
survives electromechanical coupling, and whether a DC bias can tune the negative-reflection window
continuously.*

### 3.3 Research questions

**RQ1.** Does the S₁/S₂ᵦ fold in a piezoelectric plate carry the *same* mode pair as in the elastically
identical uncoupled plate, or does coupling select a different pair?

**RQ2.** Is the near-degeneracy of the pair's stress–displacement fields — the condition that makes the
free-edge reflection coefficient approach −1 — preserved under coupling?

**RQ3.** Does a DC bias, acting through the converse piezoelectric effect, tune the fold's frequency
continuously while preserving that near-degeneracy, and over what range?

**RQ4.** What is the correct electrical boundary condition on a *cut* edge, and does the answer to RQ2
depend on it?

### 3.4 Hypotheses (each independently falsifiable)

| ID | Hypothesis | Falsified if |
|----|-----------|--------------|
| H1 | Coupling changes the *identity* of the mode pair carrying the fold | The pair and its mode shapes are unchanged at the coupled fold |
| H2 | The near-degeneracy survives coupling (ρ can still approach −1) | The traction-norm overlap of the pair is not near-degenerate |
| H3 | DC bias tunes f_ZGV monotonically and continuously | Bias moves f_ZGV non-monotonically, or destroys the degeneracy before useful tuning |
| H4 | The cut-edge electrical condition affects ρ measurably | Both extremes (φ=0, D_n=0) give the same ρ to numerical precision |

**H2 is the gate.** If H2 is falsified, the programme terminates with a negative result that is itself
publishable (see §3.12).

### 3.5 Physical mechanism claimed, and the distinction from specular reflection

The mechanism is **forward↔backward mode conversion at a traction-free edge** (M1,
`03_Negative_Reflection_Taxonomy.md`), modulated by **converse-piezoelectric stiffening**.

The distinction from ordinary specular reflection is operational, not verbal. For a computed reflection
matrix, ordinary reflection means the tangential Poynting flux of the reflected field is opposite to
that of the incident field. Negative reflection means the reflected mode is **backward** — its energy
flux is anti-parallel to its wave vector — so the tangential flux *keeps its sign* even though the beam
returns along the incident path. The test is applied to the **electromechanical flux**, not to Re(ρ):
a negative-valued reflection coefficient is explicitly *not* evidence of negative reflection, and no
claim in this proposal rests on one.

### 3.6 Governing equations

Quasi-static piezoelectricity, time-harmonic e^{−iωt}, plate |x₂| ≤ h:

$$\rho\,\omega^2 u_i + \big(c_{ijkl}u_{k,l}+e_{kij}\varphi_{,k}\big)_{,j}=0,\qquad
\big(e_{ijk}u_{j,k}-\varepsilon_{ij}\varphi_{,j}\big)_{,i}=0 .$$

**Static condensation.** Because Gauss's law has no time derivative, the potential is eliminated
exactly: **A**_{φφ}**φ** = −**A**_{φu}**u**, leaving a standard Hermitian generalised eigenproblem in
(u₁,u₂). This is what makes the flux identity valid for the reduced system (and *only* for the reduced
system — the identity is not inherited from the full multiphysics problem, and is re-derived in the
project).

**Electrical conditions on the faces** x₂ = ±h: φ = 0 (short) or D₂ = 0 (open), as specified; the bias
case additionally requires a bias circuit condition, which must be stated explicitly because a biased
plate with blocked charge is not the same as an open circuit.

**Edge** x₁ = 0: σ₁ⱼ = 0, plus the cut-face electrical condition — **both extremes compared**, which is
RQ4.

**Bias.** The linearised effective constants about the bias point: either from measured bias
dependence of c^E (preferred) or from third-order constants if available. The choice and its
justification are a month-1 deliverable, and if neither is available the problem switches material
(§3.14, risk 1).

### 3.7 Dimensionless groups

$$\Pi_1=\frac{e_{31}^2}{\varepsilon_{33}c_{11}}\approx 0.026,\qquad
\Pi_2=\frac{\omega h}{c_T},\qquad
\Pi_3=\frac{E_{\text{bias}}\,e_{33}}{c_{11}} .$$

**Π₁ contains no length.** This is the rigorous basis for the project's most important negative result:
piezoelectric coupling cannot produce a size effect, and the numerically identical −35.4 % shift at
three thicknesses confirms it. Any claim of "size-dependent piezoelectric negative reflection" is
dimensionally inadmissible on this ground alone, and saying so is part of the paper.

### 3.8 Materials and parameters

PZT-5H: c₁₁ = 126, c₁₂ = 79.5, c₁₃ = 84.1, c₃₃ = 117, c₄₄ = 23.0 GPa; e₃₁ = −6.5, e₃₃ = 23.3,
e₁₅ = 17.0 C/m²; ε₁₁ = 1700ε₀, ε₃₃ = 1470ε₀; ρ = 7500 kg/m³. Thickness 0.5 mm and 1.5 mm.
Bias 0 to ±2 kV/mm (to be confirmed against the measured breakdown and depoling limits of the chosen
material — **this must be sourced, not assumed**).
Fallback material: YX-LiNbO₃, whose backward/ZGV behaviour is experimentally established
(`10.3390/s23020648`).

### 3.9 Analytical programme

1. Re-derive modal orthogonality and the **electromechanical** flux for the biased plate.
2. Establish completeness of the mixed propagating/evanescent spectrum at the edge (needed for the
   expansion; classical result must not be transplanted — it is re-derived for the coupled operator).
3. Derive the **near-degeneracy condition as an explicit overlap integral** in the traction norm —
   this is the mathematical core of the paper and the object H2 is a statement about.
4. Perturbation theory in Π₃ for the fold shift.
5. Resolve the cut-edge electrical condition analytically where possible.

### 3.10 Numerical programme

4. Extend `code/plate_spectral.py` with bias-dependent constants and a bias circuit condition.
5. Compute dispersion across the fold and the pair's mode shapes; **the overlap gate** (§3.12).
6. Build the edge scattering matrix by modal projection with a truncated pseudo-inverse; publish the
   singular-value spectrum.
7. Sweep bias, thickness and the cut-edge electrical condition.
8. **Energy audit of every reported reflection matrix**: Σ|ρ_n|²(P_n/P_inc) = 1 to machine precision
   for the lossless case, and the residual reported alongside every number.

### 3.11 Validation

| Check | Reference | Tolerance |
|---|---|---|
| Cutoff frequencies | exact n·c/(4h) | 10⁻⁶ (already 6.3 × 10⁻⁷) |
| Dispersion | independent Rayleigh–Lamb root finder | 10⁻⁶ (already 5.3 × 10⁻⁷) |
| ZGV | independent exact-theory tracker | 10⁻⁶ (already 9.1 × 10⁻⁷) |
| Mesh convergence | N = 14 → 34 | change < 10⁻⁷ (already 8.0 × 10⁻⁹) |
| Bias → 0 | unbiased coupled result | exact recovery |
| Π₁ → 0 | decoupled elastic result | exact recovery |
| Reciprocity | ρ_ij = ρ_ji in the flux-normalised basis | machine precision |
| Energy | Σ\|ρ\|²(P/P_inc) = 1 | machine precision |

Two-sided validation is the standard aimed at. Where it is unavailable (the edge matrix itself), the
energy and reciprocity residuals substitute, and this substitution is stated as a limitation.

### 3.12 Go/no-go gate (month 2) — the decisive test

Compute the two modes of the fold and their overlap in the traction norm. **If the pair is not
near-degenerate — i.e. if the two modes cannot approximately satisfy the free-edge traction condition
between them — the project terminates and the negative result is written up.**

This gate is cheap (one eigenvalue run across the fold), it is scheduled before any scattering code is
written, and its outcome is pre-registered as publishable either way. It exists because the entire
premise of negative reflection in a coupled medium rests on an assumption nobody has checked.

### 3.13 Expected results, both outcomes pre-registered

**If H2 holds (positive):**
- The first modal scattering matrix for a coupled piezoelectric plate edge.
- A continuously bias-tunable negative-reflection window, with a quantified tuning range.
- A criterion for when electromechanical coupling destroys the effect (the overlap integral).
- The rigorous refutation of the "size-dependent piezoelectric NR" premise via Π₁.

**If H2 fails (null — and this counts as a result):**
- The statement that electromechanical coupling is a *dispersion* modifier, not a *reflection*
  modifier — which bounds the whole coupled-media negative-reflection programme, including the
  proposal's own follow-on.
- Reported prominently, not buried. A null here is more valuable than a marginal positive, because it
  tells the field not to build devices on the assumption.

### 3.14 Risk register (ranked, with mitigations)

| # | Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|---|
| 1 | Bias-dependent elastic constants for PZT-5H unavailable or unreliable | **High** | Programme-defining | Fall back to LiNbO₃; if still unavailable, present the bias treatment as a linearised model with stated assumptions and make the bias-free coupled edge result the paper |
| 2 | H2 fails — no near-degeneracy under coupling | Medium | Converts to a null result | Gate at month 2; pre-registered as publishable |
| 3 | Cut-edge electrical condition is not uniquely defined | Medium | Moderate | Compare both extremes; if the result depends on the choice, that dependence is reported as a result (RQ4) |
| 4 | Evanescent-basis ill-conditioning corrupts ρ | Medium | Moderate | Truncated pseudo-inverse, column scaling, published singular values, energy residual as the diagnostic — this exact failure mode already occurred once in this project |
| 5 | Effect too small to be interesting | Low–medium | Moderate | The −35.4 % dispersion shift already makes a null unlikely to be trivial |
| 6 | Reviewer judges the work incremental | Medium | High (publication) | Lead with the mechanism question (H2) and the negative result (Π₁), not with the material |

### 3.15 Timeline (18 months)

| Months | Work | Deliverable |
|---|---|---|
| 1–2 | Bias constitutive data; cut-edge electrical condition; **Gate: overlap test** | Gate report (go/no-go) |
| 3–5 | Modal orthogonality + electromechanical flux + completeness, re-derived | Analytical draft §2 |
| 6–9 | Edge scattering code; validation suite | Solver + validation table |
| 10–12 | Bias and thickness sweeps; energy audit; reciprocity | Full dataset |
| 13–15 | Manuscript | Submitted paper |
| 16–18 | Revision; extension to LiNbO₃; solver release | Revised paper + open code |

**Solver, data and figures released with the paper**, in the repository layout already used here
(`code/`, `data/`, `reports/figures/`).

---

## 4. P4 — Companion design (abbreviated; full specification in `07`)

**Title.** *Stiffness-controlled negative reflection at a spring-type joint between two elastic
waveguides.*

**Why paired with P1.** It reuses the same modal-projection machinery, the same energy audit, and the
same validation harness. The marginal cost of the second paper is roughly 5 months, not 18. It also
supplies the *application* that P1 lacks: a quantitative bond-quality probe.

**The claim.** The conversion efficiency is **not** a monotone function of interfacial stiffness — it
holds near unity up to a threshold κ* and then collapses, because M1 depends on a delicate
near-degeneracy. κ* defines a practical tolerance for real bonded joints, and therefore a measurable,
deployable diagnostic.

**Dimensionless group.** κ = K_∥/(ρc_T²/h), swept over 10⁻²–10².

**Validation — the strongest of any problem here.** Both limits are independently known: κ → 0 must
recover the free edge (the M1 result) and κ → ∞ the rigid bond. The result is bracketed from both
sides by known answers.

**Gate (month 4 of the combined programme).** If the efficiency varies smoothly and monotonically with
κ, H is false and the result becomes "negative reflection degrades gracefully" — weaker, still
reportable, and it converts the paper from a diagnostics paper into a robustness paper.

**Risk.** That it is a parameter study. The threshold hypothesis is the only thing standing between
this and an engineering report, and if the threshold does not exist, the paper is smaller.

---

## 5. P3 — Named escalation path (not recommended now)

Recorded here because the committee should see the strongest science alongside the safest bet.

**The question.** Negative reflection works because two modes satisfy two traction conditions. A
strain-gradient theory needs more conditions. Can they still be satisfied, and at what cost?

**Why it is the best science.** It is the only candidate where coupling changes the **structure** of
the boundary-value problem rather than its coefficients, and therefore the only one whose answer bounds
the phenomenon rather than re-parameterising it. Both outcomes — a quantified degradation law, or a
demonstration that the gradient formulation cannot support negative reflection — are new statements of
general value to the size-dependent-waves community, which currently treats only dispersion and
implicitly assumes this question away.

**Why not now.** A singular ℓ → 0 limit would invalidate the formulation; the higher-order conditions
must be derived rather than quoted; and the only fully independent validation is an expensive
finite-element comparison. The audit scores feasibility at 5.5/10 — the lowest of the five.

**Entry condition.** Enter P3 only after P1's Gate establishes that the classical near-degeneracy
survives coupling. If it does not, P3's question is moot: there is nothing for the gradient terms to
destroy.

---

## 6. What this programme is *not*

Stated explicitly, because these are the failure modes the brief warned about, and each has been
guarded against somewhere in `01`–`09`:

- **Not** a claim that negative reflection exists in any piezoelectric, piezomagnetic or
  magneto-electro-elastic material. It is established for isotropic plates. For coupled plates, the
  *fold* is computed here; the *negative reflection* is not, and the edge calculation is future work.
- **Not** a claim of a size effect from piezoelectricity. The opposite is proved by Π₁ and confirmed
  numerically.
- **Not** a claim that any simulation has been run beyond those in `code/`. The scripts and their JSON
  outputs are listed in §9.
- **Not** an experimental proposal. No experiment is claimed as performed, and P3 is explicitly
  flagged as not experimentally accessible.
- **Not** a claim about journal tier or acceptance probability. No journal has been assessed; no
  Q-rank, impact factor or acceptance likelihood is asserted anywhere in this work.
- **Not** a claim that any gap is definitively unexplored. Scopus and Web of Science were unreachable;
  all gaps are graded as degrees of evidence, and the investigator must repeat the search on the
  subscription indexes before committing.
- **Not** a conflation of backward-wave propagation with negative reflection; `03` exists specifically
  to keep them apart, and the flux test is applied throughout.
- **Not** a claim of mechanical-energy conservation in a multiphysics system; the electromechanical
  flux is stated explicitly and its residual reported.

---

## 7. Decision requested from the reviewer

1. **Approve, reject, or redirect** the P1 + P4 pairing.
2. **Confirm or overrule** the exclusion of P3 from the first 18 months. The argument for including it
   is that it is the best science; the argument against is that feasibility carries the same weight as
   originality under the agreed scoring, and P3 is the weakest on it.
3. **Confirm the material choice** (PZT-5H with LiNbO₃ fallback), or nominate a material for which
   bias-dependent elastic constants are known to exist — this is the single highest-risk input in the
   entire programme.
4. **Decide the standing of a null result.** The proposal pre-registers two null outcomes (Gate H2,
   Gate P4-κ) as publishable. If the reviewer does not accept a null result as a deliverable, this
   must be settled before month 2, because it changes the go/no-go logic.
5. **Confirm whether an experimental pathway is required.** If it is, P3 should not be scoped at all
   and the fallback P2 should replace it.

---

## 8. Summary status table

| Problem | Weighted score | Rank | Recommendation |
|---|---|---|---|
| **P1** bias-piezo | **7.80** | 1 | **Primary — execute** |
| **P4** spring joint | 7.68 | 2 | **Companion — execute** |
| P2 grading | 7.38 | 3 | Fallback if bias data fail |
| P3 strain-gradient | 7.30 | 4 | Escalation path, not now |
| P5 anisotropy | 6.88 | 5 | Not recommended at this time |

**Awards:** BEST OVERALL **P1** · BEST HIGH-RISK/HIGH-REWARD **P3** · BEST MATHEMATICALLY
TRACTABLE **P2** (P4 effectively tied).

---

## 9. Reproducibility statement

Every numerical claim in this proposal can be regenerated from this repository.

| Path | Contents |
|---|---|
| `code/plate_spectral.py` | Chebyshev spectral plate solver; isotropic and 6mm anisotropic; energy assembly; static condensation of the potential; `energy_velocity` (flux identity); `eigen_full` (mode + potential recovery) |
| `code/check_01_elastic_validation.py` | validation against exact Rayleigh–Lamb theory |
| `code/check_02_flux_and_zgv.py` | flux identity; ZGV location; piezoelectric / electrical-BC premise test; scale analysis |
| `code/gen_references.py` | builds `data/references.csv` and `data/excluded_records.csv` |
| `code/make_figures.py` | Figures 1–3 |
| `data/check_01_elastic_validation.json` | validation results |
| `data/check_02_flux_and_zgv.json` | all Phase-3 numbers, including the piezoelectric table |
| `data/references.csv` | 60 verified records with DOIs and verification status |
| `data/excluded_records.csv` | 3 rejected records, with reasons |
| `data/verification_log.md` | DOI rejection log |
| `reports/figures/fig1_dispersion_ve.*` | dispersion with the flux sign, showing the fold and the sign change |
| `reports/figures/fig2_validation_premise.*` | validation and the coupled-media premise |
| `reports/figures/fig3_premise_limits.*` | per-branch electrical sensitivity; the mode's antisymmetric potential |

Validation tolerances achieved, all reproducible by running the scripts: 6.3 × 10⁻⁷ (cutoffs),
5.3 × 10⁻⁷ (dispersion), 9.1 × 10⁻⁷ (ZGV vs independent tracker), 8.0 × 10⁻⁹ (mesh convergence),
5.3 × 10⁻⁷ median / 9.2 × 10⁻⁷ max (flux identity).

Three real bugs found and fixed during validation are recorded in `05` §7 — a piezoelectric sign
error, an incorrect Rayleigh–Lamb reference function, and branch-tracking failure at the
near-degeneracy — because their discovery is part of the evidence that the audit is genuine.

---

> ## STATUS: AWAITING INDEPENDENT RESEARCH REVIEW — NOT APPROVED
>
> This document submits a research proposal for external review. It is **not** an approved project,
> **not** a funded programme, and **not** a manuscript. No work beyond the feasibility computations
> documented in `code/` has been performed. The author of this document has no authority to approve
> the proposal and does not claim approval. Independent review is required before any execution,
> submission, or claim of results.
