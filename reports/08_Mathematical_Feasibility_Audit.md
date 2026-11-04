# 08 — Mathematical Feasibility Audit

**Document status:** Phase 6 deliverable
**Scope:** the five problems from `07`, tested against the 15 required checks.
**Rule:** a problem failing a *fundamental* requirement is marked **REJECTED** regardless of its
score in Phase 7. A problem failing a *serious but addressable* requirement is marked
**CONDITIONAL** and the condition is stated explicitly.

---

## 1. The fifteen checks

| # | Check | What "pass" means |
|---|-------|-------------------|
| F1 | Mathematical well-posedness | Existence and uniqueness of the solution of the stated BVP; no unresolved degeneracy |
| F2 | Derivable governing equations | Can be obtained from a variational principle or a stated balance law, from first principles, in this project |
| F3 | Consistent boundary conditions | All natural and essential conditions derived, not quoted; count matches the order of the system |
| F4 | No unjustified equation transplant | Every identity used in the new setting is re-derived there |
| F5 | Numerical tractability | Solvable with the solver that exists or a documented extension of it |
| F6 | Convergence and stability | A convergence path exists and is checkable |
| F7 | Validation against independent reference | An independent analytical or numerical reference exists |
| F8 | Conditioning | The linear algebra is not hopelessly ill-conditioned, or is handled by a documented method |
| F9 | Parameter availability | All material and geometric data exist and are sourced |
| F10 | Physical admissibility | The predicted regime is physically realisable |
| F11 | Energy consistency | A correct conservation statement exists and can be checked |
| F12 | Dimensional consistency | All governing groups dimensionless and complete |
| F13 | Limiting-case recovery | Known limits recover known results |
| F14 | Statistical/uncertainty treatment | Parameter uncertainty can be propagated (not applicable to all) |
| F15 | Scope and time | Achievable by one investigator in 18 months |

---

## 2. Audit table

Legend: **P** = pass · **C** = conditional (condition stated) · **F** = fail.

| Check | P1 (bias, piezoelectric) | P2 (grading) | P3 (strain gradient) | P4 (spring joint) | P5 (anisotropy) |
|---|---|---|---|---|---|
| F1 well-posedness | P | P | **C** | P | P |
| F2 derivable equations | P | P | P | P | P |
| F3 consistent BCs | **C** | P | **C** | P | P |
| F4 no transplant | P | P | **C** | P | P |
| F5 numerical tractability | P | **P (best)** | **C** | P | **C** |
| F6 convergence | P | P | **C** | P | P |
| F7 independent validation | P | **P (best)** | **C** | **P (best)** | P |
| F8 conditioning | **C** | P | **C** | **C** | **C** |
| F9 parameters | **C** | P | **C** | P | P |
| F10 physical admissibility | P | P | **C** | P | P |
| F11 energy consistency | P | P | **C** | P | P |
| F12 dimensional consistency | **P (diagnostic value)** | P | P | P | P |
| F13 limit recovery | P | P | **C** | **P (best)** | P |
| F14 uncertainty | P | P | P | P | P |
| F15 scope/time | P | **P (best)** | **C** | P | **C** |
| **Verdict** | **CONDITIONAL — PASS** | **PASS** | **CONDITIONAL — PASS, high risk** | **PASS** | **CONDITIONAL — PASS** |

---

## 3. Check-by-check findings

### F1 — Well-posedness

**P1, P2, P4, P5 pass.** Each reduces to a standard elliptic (or self-adjoint) eigenvalue problem in
x₂ with a well-studied modal expansion at the edge. Uniqueness follows from the standard argument once
completeness of the mixed spectrum is established, and completeness for the classical elastic plate is
established in the literature.

**P3 is conditional.** A fourth-order (or higher) system requires more boundary conditions, and there
are three distinct hazards:
1. The strain-gradient operator may be non-self-adjoint in some formulations, which destroys the
   orthogonality relation the modal expansion depends on.
2. The number of admissible evanescent modes may not match the number of remaining boundary
   conditions, a classic source of spurious roots.
3. The ℓ → 0 limit may be singular (a boundary layer of vanishing thickness that does not vanish in
   the boundary conditions).

**Condition for P3:** the variational formulation must be fixed in month 1, self-adjointness verified
explicitly, and the mode count matched to the condition count before any scattering calculation.
Gate 1 of P3 exists for exactly this.

### F2 — Derivability

**All five pass**, and P3 passes most strongly: strain-gradient elasticity has a genuine variational
structure, which is precisely why it was chosen over Eringen-type integral nonlocality (whose
consistent boundary conditions are an open modelling question in the literature). This is a
deliberate design choice, not an accident.

### F3 — Boundary conditions

**P2, P4, P5 pass** — classical conditions, well established.

**P1 is conditional.** The electrical condition on a *cut* edge is genuinely non-obvious. A cut face
in a poled piezoelectric ceramic has a depolarised surface layer, surface charge, and a condition
partway between open and short. No source was found stating the correct condition for the
modal-scattering problem. **Condition:** the problem must explicitly compare the two extremes
(φ = 0 and D_n = 0 on the cut face) and show whether the result depends on the choice; if it does,
that dependence is itself a result and must be reported.

**P3 is conditional** for the reason in F1: the higher-order conditions must be derived from the
variational principle in this project, not quoted. This is a hard requirement in the brief ("do not
transplant equations between theories without verifying compatibility") and it is the single biggest
analytical task in P3.

### F4 — No unjustified transplant

**P1, P2, P4, P5 pass.**

**P3 is conditional.** The classical modal-expansion-plus-pseudo-inverse method is standard for
second-order systems. Applying it to a fourth-order system requires re-establishing the completeness
relation, which is done for second-order operators by the standard argument (the operator is
self-adjoint with compact resolvent) and does **not** automatically carry over. **Condition:** the
completeness relation must be re-derived and the derivation included in the paper, not cited.

### F5 — Numerical tractability

**P2 passes best.** The existing spectral solver (`code/plate_spectral.py`) handles variable
coefficients with a trivial extension (Chebyshev differentiation applied to c(x₂)·u′), and the
dispersion and mode-shape machinery is already validated. This is why P2 is the low-risk member.

**P1 passes.** The solver already does piezoelectric plates with open/short conditions and has
already been used to locate the coupled fold and recover the potential
(`reports/figures/fig3`). The bias extension is a constitutive change, not a structural one.

**P4 passes.** The spring condition is a two-block coupling, algebraically standard.

**P5 is conditional** — cost. A full angular sweep with anisotropy at each (ω, θ) is expensive, and
the ZGV locus extraction requires a 2-D root search in (k, θ) rather than a 1-D one. Tractable, but it
is the most computationally demanding of the five.

**P3 is conditional, for a structural reason.** The state dimension doubles (displacement and
displacement gradient), the boundary rows must be built from derived higher-order conditions, and the
eigenvalue problem is no longer standard. It is tractable (Chebyshev handles it), but it is easily
the most error-prone.

### F6 — Convergence and stability

Established paths exist for P1, P2, P4, P5: spectral order N, retained-mode count M, and frequency
step are all independently checkable, and N-convergence is already demonstrated here
(mesh-to-mesh change 8.0 × 10⁻⁹ at the ZGV).

**P3 is conditional:** convergence must be checked in N, M, the internal length, *and* the ℓ → 0
limit simultaneously, and the limit is the hard part.

### F7 — Independent validation

**P2 and P4 pass best.**
- P2: the graded dispersion can be checked against a Thomson–Haskell compound-matrix solution with
  many thin sublayers — a completely independent method, and a very strong check.
- P4: the κ → 0 (free) and κ → ∞ (rigid bond) limits are both independently known, which brackets the
  entire result.

**P1 passes** with the bias → 0 limit plus the existing classical-plate validation.

**P5 passes** with the isotropy limit and the published LiNbO₃ backward-wave frequency.

**P3 is conditional:** the only fully independent check is a finite-element solution at one parameter
point, which is expensive; the rest are self-consistency checks (reduction, energy). This is a real
weakness of P3 and is reflected in its risk score.

### F8 — Conditioning

**All five are conditional**, because all five use the modal projection method at an edge, whose
matrix of evanescent-mode traction vectors is **inherently ill-conditioned** — the evanescent modes
decay exponentially and their columns become nearly linearly dependent as the retained set grows.

This is the standard, well-known difficulty of the field, and it has standard remedies:
- Moore–Penrose pseudo-inverse with singular-value truncation;
- scaling each column by the mode's norm at the edge;
- orthogonalisation of the retained basis;
- independent verification by energy conservation (an ill-conditioned inversion typically violates
  energy, so the energy residual is a *diagnostic* of conditioning failure, which is useful).

**This is the same difficulty that produced the spurious branch-tracking results in this project**
(`05` §7), so it is not theoretical: it has already bitten once.

**Condition (all five):** publish the singular-value spectrum of the projection matrix and the
energy residual alongside every reported reflection coefficient.

### F9 — Parameters

**P2, P4, P5 pass** — all materials are standard and well tabulated.

**P1 is conditional.** Bias-dependent elastic constants for PZT-5H under a DC field are needed but
were **not verified to exist in a usable form** during this audit. If they do not, the fallback is
LiNbO₃ (better characterised under bias) or a linearised analytical bias model with stated
assumptions.

**P3 is conditional.** Internal lengths ℓ for strain-gradient elasticity span an order of magnitude
across sources and are not agreed. **Condition:** the dimensionless scan Π_ℓ ∈ [0, 0.5] must be the
core result, with material-specific values presented as an illustration with uncertainty bars. This
is the honest treatment, and it also makes the result robust to better ℓ data later.

**P5 passes** — LiNbO₃ constants are standard, and the choice of material is justified by published
experiment.

### F10 — Physical admissibility

**P1, P2, P4 pass.** All regimes are ordinary: MHz frequencies, mm-to-sub-mm plates, achievable bias
fields.

**P5 passes** — μm-to-mm LiNbO₃ plates at MHz are standard SAW/BAW territory.

**P3 is conditional.** The *numerical* regime is accessible but the *experimental* one is not:
ℓ/h ~ 0.1 requires h ~ 10–100 µm at ℓ ~ 1–10 µm, i.e. hundreds of MHz, which is possible but beyond
routine laboratory guided-wave measurement. **Condition:** P3 must be presented as a theory-and-
numerics paper, explicitly not claiming experimental accessibility. If the committee requires an
experimental pathway, P3 should be dropped in favour of P2.

### F11 — Energy consistency

**P1, P2, P4, P5 pass.** The classical / piezoelectric flux is well defined and the residual is
checkable.

**P3 is conditional but in a *good* way:** the gradient energy includes the double-stress contribution,
which makes the conservation statement stricter and consequently a *stronger* validation of the
derived boundary conditions. If the derived conditions are wrong, energy will not balance. This turns
the hardest part of P3 into a self-checking one.

### F12 — Dimensional consistency, and one diagnostic result

**All pass.** But P1 carries a specific, already-obtained diagnostic value: the ratio
Π₁ = e₃₁²/(ε₃₃c₁₁) = 0.026 contains **no length**, which is a complete and rigorous argument that
piezoelectric coupling alone cannot produce a size effect. This was confirmed numerically (identical
−35.4 % shift at 0.05, 0.5 and 1.5 mm). Any proposal that claims otherwise is dimensionally wrong,
and P1 is the problem that establishes this cleanly.

### F13 — Limiting-case recovery

| Problem | Limits available | Strength |
|---|---|---|
| P1 | bias → 0; Π₁ → 0 (decoupled) | strong |
| P2 | Π_g → 0 (homogeneous) | very strong |
| P3 | ℓ → 0 (classical M1) | **the hard one — may be singular** |
| P4 | κ → 0 (free); κ → ∞ (rigid) | **strongest of all five** |
| P5 | isotropy (cubic → isotropic limit) | strong |

### F14 — Uncertainty

All five can propagate parameter uncertainty (material constants, ℓ, κ, incidence angle) by direct
resampling of deterministic solves — no stochastic method needed. **Pass for all.**

### F15 — Scope and time

**P2 passes best** — it is a single solver extension plus a parameter study with an analytical
criterion, comfortably 18 months for one investigator.

**P4 passes** — a matrix formulation with two known limits; the risk is scientific thinness, not time.

**P1 passes** — the solver exists and the fold is already located; the new work is the edge matrix.

**P5 is conditional** — the angular sweep is expensive and the analytical part is algebraically
heavy; realistic but tight.

**P3 is conditional** — the analytical derivation alone (variational conditions, self-adjointness,
boundary layer) is 6–9 months, leaving 9–12 for numerics and writing. Achievable for a strong
investigator, tight for one working largely alone.

---

## 4. Rejections

**No problem is REJECTED by this audit.** This is stated deliberately and should be read critically:

- P3 was the closest to rejection, on F1/F4/F7. It is retained as CONDITIONAL because
  (i) the variational structure is genuine, (ii) two independent gates exist and are cheap, and
  (iii) the failure mode is informative rather than fatal.
- If Gate 1 of P3 (the ℓ → 0 limit) fails, **P3 is automatically REJECTED** and must not be pursued;
  that condition is now on the record.

Two candidates from Phase 4 **were** rejected before reaching this audit, and the reasons belong here
because they are mathematical, not merely editorial:

1. **"Nonlocal (Eringen integral) edge scattering" — REJECTED.** The consistent boundary conditions
   for integral nonlocal elasticity are not settled; a modal expansion requires a well-defined
   boundary-value problem. Pursuing it would have violated the brief's rule against transplanting
   equations between theories without verification. This is why P3 uses a gradient (variational)
   theory instead.
2. **"Exceptional-point-based negative reflection" — REJECTED as an NR problem.** The reflection
   mechanism is coalescence, not forward↔backward conversion; the phenomenon requires gain; and the
   mathematical content is that of non-Hermitian eigenvalue theory, not of waveguide scattering. It
   would not be a negative-reflection paper.

---

## 5. Summary verdict

| Problem | Verdict | Single biggest mathematical risk | Mitigation |
|---|---|---|---|
| **P1** | PASS (conditional) | Edge electrical condition not unique | Compare both extremes; report dependence |
| **P2** | **PASS (cleanest)** | Effect may be marginal | Analytical bandwidth criterion is the deliverable |
| **P3** | PASS (conditional, high risk) | ℓ → 0 limit may be singular | Gate 1, month 3 |
| **P4** | **PASS (cleanest)** | Thin science, not thin maths | Threshold hypothesis |
| **P5** | PASS (conditional) | Cost + possible vacuity | Gate, month 5 |

**Mathematical feasibility ranking** (for Phase 7's 20 % weight): **P2 > P4 > P1 > P5 > P3.**

The audit does not reject anything outright, but it does establish something important for the
ranking: **P2 and P4 are the only problems with two-sided, independent validation** (homogeneous
limit + independent method for P2; free and rigid limits for P4). P3's mathematical content is the
deepest but its validation is the weakest — a genuine and uncomfortable trade-off that the ranking
must reflect rather than hide.
