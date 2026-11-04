# 09 — Comparative Research Ranking

**Document status:** Phase 7 deliverable
**Weighting (fixed by the brief, not by me):** originality 20 % · mathematical feasibility 20 % ·
physical significance 20 % · Gupta alignment 15 % · numerical reproducibility 10 % · application
10 % · scope 5 %.
**Hard rule:** any problem that fails a *fundamental* feasibility requirement is **REJECTED**
regardless of its total. Per `08`, no problem is rejected; but P3 carries a standing auto-rejection
condition.

---

## 1. Scoring anchors

Scores are only meaningful if the scale is fixed before the numbers. These are the anchors used.

| Criterion | 10 | 7 | 5 | 3 |
|---|---|---|---|---|
| **Originality** | Establishes a new mechanism | First treatment of a known mechanism in a structurally new setting | Known effect, new parameter or material | Straight material substitution |
| **Mathematical feasibility** | All 15 audit checks pass cleanly; two-sided independent validation | One or two addressable conditionals | A structural conditional (well-posedness, BC count, or singular limit) | Ill-posed or method not transferable |
| **Physical significance** | Changes how the phenomenon is understood or bounded | Confines the phenomenon usefully | Adds a design knob | Cosmetic |
| **Gupta alignment** | Squarely inside his verified strength *and* extends it | Adjacent to his strengths | Overlaps one subfield | Outside his record |
| **Numerical reproducibility** | Fully reproducible with already-validated code + an independent reference method at every parameter point | Reproducible with a documented extension | Substantial new numerical infrastructure | Not reproducible |
| **Application** | Deployable measurement or device with an identified user | Plausible device concept | Long-term / speculative | None |
| **Scope** | Comfortably 18 months for one investigator | Achievable but tight | Needs a second investigator or 3 years | Not achievable |

---

## 2. Scores

| Criterion (weight) | P1 bias-piezo | P2 grading | P3 strain-gradient | P4 spring joint | P5 anisotropy |
|---|---|---|---|---|---|
| Originality (0.20) | **7.0** | 5.0 | **9.0** | 5.0 | 7.0 |
| Mathematical feasibility (0.20) | 7.0 | **9.0** | 5.5 | **9.0** | 6.5 |
| Physical significance (0.20) | **8.0** | 6.0 | **9.0** | 6.0 | 7.0 |
| Gupta alignment (0.15) | 9.0 | 8.5 | 8.5 | **9.5** | 7.0 |
| Numerical reproducibility (0.10) | 8.5 | **9.5** | 5.5 | 9.0 | 7.0 |
| Application (0.10) | 8.0 | 7.0 | 5.0 | **9.0** | 7.0 |
| Scope (0.05) | 8.0 | **9.0** | 5.5 | **9.0** | 6.5 |
| **WEIGHTED TOTAL** | **7.80** | 7.38 | 7.30 | 7.68 | 6.88 |
| **Rank** | **1** | 3 | 4 | **2** | 5 |

### Arithmetic (for audit)

- **P1** = 0.20(7.0) + 0.20(7.0) + 0.20(8.0) + 0.15(9.0) + 0.10(8.5) + 0.10(8.0) + 0.05(8.0)
  = 1.40 + 1.40 + 1.60 + 1.35 + 0.85 + 0.80 + 0.40 = **7.80**
- **P2** = 0.20(5.0) + 0.20(9.0) + 0.20(6.0) + 0.15(8.5) + 0.10(9.5) + 0.10(7.0) + 0.05(9.0)
  = 1.00 + 1.80 + 1.20 + 1.275 + 0.95 + 0.70 + 0.45 = **7.375 → 7.38**
- **P3** = 0.20(9.0) + 0.20(5.5) + 0.20(9.0) + 0.15(8.5) + 0.10(5.5) + 0.10(5.0) + 0.05(5.5)
  = 1.80 + 1.10 + 1.80 + 1.275 + 0.55 + 0.50 + 0.275 = **7.30**
- **P4** = 0.20(5.0) + 0.20(9.0) + 0.20(6.0) + 0.15(9.5) + 0.10(9.0) + 0.10(9.0) + 0.05(9.0)
  = 1.00 + 1.80 + 1.20 + 1.425 + 0.90 + 0.90 + 0.45 = **7.675 → 7.68**
- **P5** = 0.20(7.0) + 0.20(6.5) + 0.20(7.0) + 0.15(7.0) + 0.10(7.0) + 0.10(7.0) + 0.05(6.5)
  = 1.40 + 1.30 + 1.40 + 1.05 + 0.70 + 0.70 + 0.325 = **6.875 → 6.88**

---

## 3. Justification of the contested scores

The arithmetic is trivial; the defensibility of the numbers is not. These are the scores a referee
would challenge, and why they stand.

### P1 originality = 7.0 (not 9)

It is *not* a new mechanism. M1 is established in isotropic plates. What is new is **structural**: in
a coupled plate the fold may be carried by a *different* mode pair, and the edge problem gains an
extra condition (the electrical one), so the near-degeneracy that M1 relies on is not inherited — it
must be re-established. That is a first treatment in a structurally new setting: **7**.

It is *not* 9 because the phenomenon and its interpretation are unchanged; only its carrier moves.

### P3 originality = 9.0 and mat. feasibility = 5.5

**9.0** because the question is whether the mechanism survives *at all* once the boundary-value
problem changes order. That is a mechanism-level question, and either answer is a new statement.

**5.5** because three of the fifteen checks are structural, not cosmetic: a possible singular ℓ → 0
limit (F1, F13), the need to derive rather than quote higher-order conditions (F3, F4), and the
absence of a fully independent validation method (F7). The audit's own conclusion was that P3 has
"the deepest mathematical content but the weakest validation" — and the weights reward feasibility
exactly as much as originality, so this cannot be waved away.

### P4 Gupta alignment = 9.5

This is the highest single score in the table, and it is justified by specificity: **imperfect
interfaces and spring-type contacts are literally his published specialism**, applied here to a
new observable. If alignment alone decided the ranking, P4 would win.

### P2 numerical reproducibility = 9.5

The only problem whose validation reference (Thomson–Haskell compound matrix for a graded plate) is
**independent at every parameter point**, not merely at two endpoints. Combined with an existing
validated solver, this is the most reproducible of the five by a clear margin.

### P5 application = 7.0, originality = 7.0, feasibility = 6.5

The criterion correction is genuinely valuable but it is a *clarification*, and the audit flagged a
real possibility that the wave-vector/energy-vector misalignment is negligible near the fold, which
would make the paper vacuous (Gate at month 5). Speculative but well-founded: **7**.

---

## 4. Named awards

### BEST OVERALL — **P1**: *Bias-tuned negative reflection in a polarized piezoelectric plate*

The only problem that is simultaneously in the top two on originality-adjacent criteria (physical
significance 8.0), the strongest on alignment among the high-novelty options (9.0), and free of any
structural mathematical obstacle. It also inherits the largest amount of already-completed,
already-validated work: the coupled solver exists, the fold has already been located at
f = 1.130 MHz / k = 809 m⁻¹, the electrical-insensitivity of the mode has already been measured, and
the mode's antisymmetric potential has already been recovered. The project starts from month 4 of a
typical cold start.

**Its weakness, stated plainly:** originality 7.0 is not 9.0, and the bias-dependent constitutive data
are unverified. Both are addressed in `10` §4.

### BEST HIGH-RISK / HIGH-REWARD — **P3**: *Higher-order boundary conditions and the fate of forward↔backward conversion*

Top score on originality (9.0) and joint-top on physical significance (9.0), lowest on feasibility
(5.5). It is the only candidate whose answer would **bound the entire phenomenon** — either negative
reflection survives into the gradient regime with an explicit penalty law, or it does not exist there
at all. Either result is more valuable than another parameter study.

**Its risk is real and quantified:** a singular ℓ → 0 limit would invalidate the formulation, and the
audit places that at "medium–high". Gate 1 exists for exactly this, and P3 is auto-rejected if Gate 1
fails.

### BEST MATHEMATICALLY TRACTABLE — **P2**: *Grading as mode-shape control*

Passes 15/15 checks with no structural conditional, has an independent validation method at every
parameter point, reuses the existing solver with a purely algebraic extension, and has a pre-registered
go/no-go at month 3. **P4 is effectively tied** (also 9.0 on feasibility, with the strongest
limits-based validation of the five) and would be the better choice if an industrial partner for bond
testing is available; P2 edges ahead only on interior validation and on scope.

---

## 5. What the ranking does *not* say

- It does not say P3 is a bad project. It says P3 is a worse *bet* under a weighting that values
  feasibility as much as originality. A committee that weights originality more heavily — or that can
  assign a second investigator — should prefer P3. That trade-off is explicit, and the weights were
  fixed by the brief beforehand, which is the point of fixing them.
- It does not say P5 is uninteresting. It says P5 is speculative, expensive, and possibly vacuous.
- It does not rank the *quality* of the science. P3 is the best science here; P1 is the best project.
  Those are different statements and both belong on the record.
- **No score is a prediction of publication.** No journal has been assessed for suitability, scope or
  quality tier; per the brief, no Q-rank or acceptance probability is claimed anywhere in this work.

---

## 6. Sensitivity of the ranking to the weights

A ranking that flips under small weight changes is not a ranking. Check: re-weight to
originality 30 %, feasibility 15 %, significance 20 %, alignment 15 %, reproducibility 5 %,
application 10 %, scope 5 %:

| Problem | New total | Rank |
|---|---|---|
| P1 | 0.30(7)+0.15(7)+0.20(8)+0.15(9)+0.05(8.5)+0.10(8)+0.05(8) = 2.10+1.05+1.60+1.35+0.425+0.80+0.40 = **7.73** | 2 |
| P2 | 1.50+1.35+1.20+1.275+0.475+0.70+0.45 = **6.95** | 5 |
| P3 | 2.70+0.825+1.80+1.275+0.275+0.50+0.275 = **7.65** | 3 |
| P4 | 1.50+1.35+1.20+1.425+0.45+0.90+0.45 = **7.28** | 4 |
| P5 | 2.10+0.975+1.40+1.05+0.35+0.70+0.325 = **6.90** | 6 (of 5, i.e. last) |

P1 **stays first**; P3 rises from 4th to 3rd and closes to within 0.08 of P1. So:

- **P1's lead is robust** under a doubling of the originality weight. That is the single most
  important fact for the recommendation.
- **P3's position is weight-sensitive**, as expected for a high-risk problem, and it becomes a serious
  competitor under an originality-heavy rubric.
- **P2's and P4's positions are weight-sensitive in the other direction** — they are safe bets that a
  glamour-weighted committee would demote.

The recommendation in `10` accounts for this: P1 as primary, P4 as companion (they share machinery),
P3 as the named escalation path.

---

## 7. Rejections by hard rule

| Problem | Hard-rule status |
|---|---|
| P1 | none |
| P2 | none |
| P3 | **auto-reject if Gate 1 (ℓ → 0 regularity) fails** — condition now on record |
| P4 | none |
| P5 | **auto-reject if Gate (misalignment < 0.1°) fails** |

Two Phase-4 candidates were rejected before scoring, on mathematical grounds, and are recorded in
`08` §4: Eringen integral nonlocal edge scattering (no settled boundary conditions) and
exceptional-point reflection (a different mechanism requiring gain).
