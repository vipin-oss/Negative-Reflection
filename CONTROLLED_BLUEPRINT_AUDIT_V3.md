# CONTROLLED_BLUEPRINT_AUDIT_V3.md

**Independent mathematical audit — Blueprint V2 / V2.0**
*Piezothermoelastic Anti-Reflection and Advanced Interface Wave Dynamics*

| | |
|---|---|
| **Auditor role** | Independent mathematical referee and verification specialist. Not the principal research designer. |
| **Document audited** | `FINAL RESEARCH BLUEPRINT — V2.0`, *Analytical Anti-Reflection and Energy Transmission of Coupled Waves in Layered Piezothermoelastic Media with Generalized Interface Conditions* |
| **Prior audit** | `ARENA_V2_CRITICAL_AUDIT.md` (verdict REDESIGN), 23 evidence rows |
| **This audit** | Re-derivation of every disputed claim. Five of the prior audit's own conclusions are **rejected or corrected below**. |
| **Evidence** | `code/check_05_dispute_resolution.py` → `data/check_05_dispute_resolution.json`. Eight enforced self-checks, all passing. |
| **Not produced** | No manuscript, no rewritten blueprint, no full numerical experiment campaign. |

> **All numeric claims in this report are reproducible.** Run
> `python3 code/check_05_dispute_resolution.py`. Self-checks print first and
> must pass before any table is quoted.

---

## 1. Executive scientific verdict

**The research direction is scientifically defensible and should be preserved.**
It is not fundamentally impossible, and the prior audit's instruction to
"redesign" was an overreach on three of its four critical findings.

**Verdict on the blueprint: CONDITIONAL REVISION, not REDESIGN.**

I find **four genuine mathematical/physical defects** in the blueprint, of which
**two are mandatory corrections** and two are scope corrections. I also find
that **five claims in the previous Arena audit are wrong or over-generalised**,
and I reject them here with derivations. The net effect is that the blueprint's
core problem survives, while its headline framing (quarter-wave matching) must
be replaced with a different and stronger analytical objective.

The single most important finding of this audit is **not** a defect of the
blueprint at all — it is a structural fact that determines what the paper can
be about:

> **In an anti-plane (SH) piezothermoelastic junction there is no reflection
> minimum to find.** The modulus |R| is strictly monotone in the transverse
> wavenumber and decreases to the bulk-impedance mismatch as k → 0. The
> *elastic* SH junction, by contrast, has a genuine interference minimum of
> |R| ≈ 8.9×10⁻⁷. **Turning on piezoelectric coupling destroys the elastic SH
> reflection minimum.** Consequently any reflection-minimum claim in a
> piezothermoelastic junction must be an in-plane P–SV phenomenon, where two
> polarisations can interfere.

That result is (a) rigorous within its stated model, (b) inside applied
mathematics, (c) inside Dr. Gupta's method base, and (d) it *constrains* the
blueprint rather than replacing it. It also happens to explain why the
blueprint's three-region design was never going to behave like a quarter-wave
coating.

---

## 2. Existing blueprint strengths worth preserving

Preserve all of the following. None of them is defective.

1. **Plane-strain piezothermoelastic formulation** with (u_x, u_z, φ, θ) —
   correct field set for the intended problem, correct for Dr. Gupta's identity.
2. **B0 → B1 → B2 interface hierarchy.** The ordering (welded baseline first,
   imperfect second, surface elasticity last and gated) is right.
3. **Refusal to adopt Gurtin–Murdoch as the principal mechanism.** The gate
   (surface stress negligible unless η_s = O(1)) is correct and the numbers
   below confirm it (§3.1).
4. **The null-hypothesis framing (H0).** This is the strongest single decision
   in the document and §4 below shows it pays: the most defensible contribution
   available *is* a bounded-negative result.
5. **The insistence on counting equations, rank and conditioning** rather than
   asserting "four modes". The instinct is correct; §3.3 gives the exact count.
6. **The anti-reflection / negative-reflection distinction.** Correct and
   essential; must be kept.
7. **All three cited works are real and correctly described**
   (verified against Crossref, §5).
8. **§13's refusal to claim device performance.** Correct.
9. **§11's energy-balance gate.** Correct as a *requirement*; §3.4 gives the
   exact form it must take.

---

## 3. Verified critical errors

Each is classified per Rule 7. Only **mathematical error** and **physical
inconsistency** justify mandatory correction.

### 3.1 Surface elasticity gate is quantitatively closed — *optional improvement*, not an error

The blueprint states η_s ~ k_∥ μ_s / μ_bulk and asks whether it is negligible.
It is decidable, and it is decidable by three orders of magnitude:

```
η_s = 2π ℓ_s/λ ,     ℓ_s = μ_s/μ_bulk
```

| μ_s | ℓ_s | η_s at λ = 1 mm | η_s at λ = 1 μm | λ for η_s = 1 |
|---|---|---|---|---|
| 1 N/m | 3.9×10⁻¹¹ m | 2.4×10⁻⁷ | 2.4×10⁻⁴ | 2.4×10⁻¹⁰ m |
| 5 N/m | 1.9×10⁻¹⁰ m | 1.2×10⁻⁶ | 1.2×10⁻³ | 1.2×10⁻⁹ m |
| 10 N/m | 3.8×10⁻¹⁰ m | 2.4×10⁻⁶ | 2.4×10⁻³ | 2.4×10⁻⁹ m |

*(μ_bulk = 26 GPa; Part C of `check_04`. Constants illustrative, but the
conclusion is a six-decade separation.)*

**Verdict: the previous audit's C1 is ACCEPTED but downgraded in severity.** It
is not a mathematical error; it is a gate the blueprint already flagged and can
close in one sentence. Surface elasticity is **not** universally negligible — it
is negligible *at millimetre wavelengths*, and the table states exactly where it
ceases to be. Paragraph §6.3 should be replaced by this table and study S6
deleted.

### 3.2 The quarter-wave identity is a complete classification, not a sufficient condition — *mathematical completeness error in the blueprint*

The blueprint presents `Z_m² = Z₁Z₂`, `k_m h = (2n+1)π/2` as a benchmark "under
its own assumptions" without stating what those assumptions are or whether the
condition is also *necessary*. It is necessary, and this is provable in closed
form.

**Theorem (verified symbolically and numerically).** For two semi-infinite media
of real characteristic impedance Z₁, Z₂ and one homogeneous lossless layer of
real impedance Z_m and thickness h, with θ = k_m h,

```
Re Z_in = Z₂ Z_m² / (Z₂² sin²θ + Z_m² cos²θ)
Im Z_in = Z_m (Z_m² − Z₂²) sinθ cosθ / (Z₂² sin²θ + Z_m² cos²θ)
```

Consequently:

- `Im Z_in = 0` **iff** `sinθ cosθ = 0`, i.e. θ = mπ/2.
- At those points `Re Z_in` takes only the values **Z₂** (m even) or
  **Z_m²/Z₂** (m odd).
- Hence the set of *real* values Z_in attains is exactly `{Z₂, Z_m²/Z₂}` when
  Z_m ≠ Z₂, and `{Z₂}` identically when Z_m = Z₂.

**Therefore exact zero reflection occurs iff Z₁ = Z₂ (trivial) or
`Z_m² = Z₁Z₂` with `k_m h = (2n+1)π/2`.**

**Numerical confirmation** (Part Q1):

| check | value |
|---|---|
| closed form vs direct complex evaluation, rel. error | 1.81×10⁻¹⁵ |
| \|R\| at quarter-wave point | 8.97×10⁻¹⁷ |
| max \|Im Z_in\| when Z_m = Z₂ (degenerate) | 0.00×10⁰ |

**Verdict: ACCEPTED, and it strengthens the blueprint's §9.1 while destroying
its novelty claim.** The blueprint should state the *uniqueness* result — it is
a genuine one-line theorem — and must stop treating quarter-wave matching as the
paper's analytical contribution, because uniqueness of a known result is still a
known result.

### 3.3 Boundary-condition count is arithmetically wrong for the three-region model — *mathematical error (mandatory)*

**Derivation.** For the four-field model the k_z polynomial has degree n = 8:

```
u_x, u_z : 2nd order each, coupled -> 4 roots
φ        : quasi-static Gauss       -> 2 roots
θ        : 2nd order in z           -> 2 roots         n = 8
```

One welded interface must impose **n = 8** scalar conditions
(u_x, u_z, σ_xz, σ_zz, φ, D_z, θ, q_z).

- Two half-spaces: 4 reflected + 4 transmitted = 8 unknowns ⇒ C = 8. ✔
- Three regions: 4 + 8 + 4 = **16 unknowns**, 2 interfaces ⇒ **C = 8 per
  interface, 16 total**. ✔

**Verdict: ACCEPTED as a mandatory correction.** The blueprint's §5 warns about
carrying the count over and then cites prior art that used **seven**. Seven
leaves rank deficiency one. Either the blueprint has mischaracterised Guha et
al., or that formulation fixes a field by an extra assumption that must be
identified. In either case the three-region count is 16 and cannot be inherited.

**A stronger, independent reason to distrust four-field counts (new in this
audit).** My own attempt at the anti-plane reduction produced a system with
**4 amplitudes but only 3 independent conditions** once a tangential spring
replaced w-continuity. The count is not a formality: an incomplete interface law
set produces a rank-deficient matrix whose null vector looks like a physical
mode. This is the same failure mode that generated the spurious "exact
reflection zero" retracted in §4.3.

### 3.4 The energy budget is displayed before it is defined, and in the intended regime it is trivial — *physical inconsistency (mandatory)*

The blueprint's §8 displays `R_E + T_E + A_E = 1` and *then* attaches the
caveat that it holds only after a suitable balance is derived. The ordering is
wrong typographically and, more importantly, the term A_E is not of the size the
document assumes.

**Verified passivity statement** (Part Q6, exact three-media solver):

> For a passive LTI junction `R_E + T_E + A_E = 1` with `A_E ≥ 0`, and
> `A_E = 0` **iff** every material in the finite region is lossless.
> Hence (a) `R_E → 0` with `A_E > 0` is allowed, so **reflectionlessness does
> not imply perfect transmission**; and (b) `R_E = 0` together with `T_E = 1`
> forces `A_E = 0`, i.e. is a **losslessness** statement, not a matching
> statement.

**Measured order of vanishing** (matched interphase, `c = c₀(1 − iη)`):

| η (loss tan) | R_E | R_E/η² | A_E | A_E/η | T_E |
|---|---|---|---|---|---|
| 0 | 1.95×10⁻¹⁰ | — | 0 | — | 1.000000 |
| 1×10⁻⁴ | 8.92×10⁻⁹ | 0.892 | 3.540×10⁻⁴ | 3.540 | 0.999646 |
| 1×10⁻³ | 6.66×10⁻⁷ | 0.666 | 3.529×10⁻³ | 3.529 | 0.996470 |
| 1×10⁻² | 6.37×10⁻⁵ | 0.637 | 3.425×10⁻² | 3.425 | 0.965690 |

**This is the most operationally important table in the report.** Reflected
power is **second order** in the loss; absorbed power is **first order**. A
lossy interphase is therefore a *nearly reflectionless but noticeably
attenuating* window, and the two figures of merit diverge at first order while
the reflection figure barely moves. **An optimisation that reports only R_E will
certify a bad design as excellent.**

**Verdict: ACCEPTED as mandatory.** H1 must name which objective it optimises
and must report R_E, T_E and A_E separately.

### 3.5 The thermal field is not inert — *the previous audit was wrong*

See §4.2. This is recorded here because the blueprint must **not** delete the
temperature field on the strength of the earlier claim.

---

## 4. Disputed Arena claims and their resolution

Five prior-audit claims are corrected. This section is the core of the
"controlled" character of this audit: the prior report's confidence exceeded its
evidence in five places, and two of those five were errors of my own making.

### 4.1 Prior C2 — "with loss, only gain can close the reflection floor" — **REJECTED**

**Prior claim.** For a lossy interphase, `|R|_floor = |Re Z_in − Z₁|/(Re Z_in + Z₁)`
is strictly positive and only negative resistance (gain) can close it.

**What is actually true.** That formula is correct *for a one-knob family*: the
interphase material fixed, only h free. It is **not** a general impossibility
result. When the interphase impedance is also free, `Re Z_in = Z₁` is attainable
with loss present.

**Evidence.** Sweeping the layer impedance ratio and loss (Part Q1):
`Re Z_in < Z₁` occurs in **198 of 366** configurations tested. Refining the
crossing by bisection gives `|R|` at machine zero with absorption present.
Conversely, for the *fixed* matched material the floor is positive and grows
with loss:

| loss tan | min \|R\| (fixed material, h free) |
|---|---|
| 0 | 1.77×10⁻⁷ |
| 1×10⁻⁴ | 3.61×10⁻⁵ |
| 1×10⁻³ | 3.53×10⁻⁴ |
| 1×10⁻² | 3.44×10⁻³ |

**Corrected conclusion.** The obstruction is **conditional on the parameter
family**, exactly as the blueprint's §9.2 suspected. It is a legitimate
proposition only when stated as: *for a fixed interphase material with a
dissipative response, the reflected power is bounded below by a positive
function of the loss.* It must never be stated as "passive interfaces cannot
achieve anti-reflection".

**Severity: unsupported claim in the prior audit.** Blueprint needs no
correction from this; the prior audit does.

### 4.2 Prior C3 — "the thermal field is inert" — **REJECTED as a general claim**

**Prior claim.** ε ≈ 3×10⁻³ and α_th ≈ 10⁻¹⁰ Np/mm, therefore temperature is
inert and the four-field and three-field budgets are identical.

**What is actually true.** The statement is a statement about **τ_q**, not about
piezothermoelasticity. Three-part resolution:

**(i) Fourier limit (τ_q = 0): agreed.** The thermal root is purely diffusive,
`k = (1+i)/δ`, `δ = √(2K/ωρc_e)`. At 1 MHz, δ ≈ 0.4 μm in PZT. No thermal
branch propagates. Self-check: LS with τ→0 reproduces the Fourier decay rate to
**0.00×10⁰** relative error.

**(ii) Lord–Shulman: the branch does propagate, above `ωτ_q = 1`.** The exact
dispersion relation is `k² = ω(i + ωτ_q)ρc_e/K`. For `ωτ_q ≫ 1` this gives a
**real** wavenumber with

```
v_T = √(K/(ρ c_e τ_q)) ,        onset  f_on = 1/(2πτ_q)
```

Self-check: the dispersion relation reproduces this `v_T` to **1.25×10⁻⁹**.
So a propagating thermal branch exists; whether it is reachable is a question of
the material's τ_q.

**(iii) The real result is a constraint, not a verdict.** Causality requires
`v_T ≤ c_L`, i.e. `τ_q ≥ τ_min = κ/c_L²`. Combining with the onset frequency:

```
v_T/c_L = √(2π f_on τ_min)
```

| material | ε | τ_min | τ_min from data | v_T/c_L at τ=10τ_min | f_on |
|---|---|---|---|---|---|
| PZT-5H (CTS/NASA) | 1.20×10⁻⁴ | 3.66×10⁻¹⁴ s | | 0.3162 | 4.35×10¹¹ Hz |
| BaTiO₃ (illustrative) | 2.45×10⁻³ | 2.47×10⁻¹⁴ s | | 0.3162 | 6.44×10¹¹ Hz |
| Aluminium (reference) | 7.10×10⁻³ | 2.11×10⁻¹² s | | 0.3162 | 7.54×10⁹ Hz |

**Speed and onset cannot be chosen independently.** Fixing τ_q to put the onset
inside an ultrasonic band drives `v_T/c_L` down in proportion; fixing τ_q for a
usable `v_T` pushes the onset into the GHz range. This is a **publishable
structural statement** and it replaces both the prior audit's "inert" verdict
and the blueprint's unexamined assumption that τ_q is a free knob.

**Caution on the τ_q values themselves.** A naive τ_q = 1 ps for aluminium gives
`v_T/c_L = 1.45`, i.e. **superluminal** — self-check S5 caught this and the
values above were re-derived as τ_q = 10 τ_min. Any blueprint statement of τ_q
must be shown causal or it is unphysical.

**Severity: physical inconsistency in the *prior audit*.** The blueprint should
**retain** the temperature field and add the constraint relation.

### 4.3 This audit's own spurious "exact reflection zero" — **RETRACTED**

Recorded because it is the most instructive failure in the whole exercise.

An intermediate version of this audit found an **exact** reflection zero at
`k = 2890.43 m⁻¹` for an anti-plane piezothermoelastic junction, frequency
independent, with `k/k_c = 1.090428` constant across decades. It looked like a
new matching condition.

**It was an artefact of two compounding errors:**

1. **Sign error in the dispersion relation.** I had
   `κ² = k² − ω²ρ/T`. The correct relation for `e^{i(kx + κz − ωt)}` is
   `κ² = ω²ρ/T − k²`. My version had *no propagating branch at all* in the
   wavenumber range I was scanning, so I was solving a purely evanescent
   problem and calling the result "propagation".
2. **Wrong matrix row.** The `σ_z3` and `D_z` rows used `κ` where the correct
   `dψ/dz` term belongs, which mixed the two branches.

**After correction** (exact tests: identical media give `|R| = 1.5×10⁻¹⁶`;
near-normal incidence gives `|R| = 0.114512` against an independent
bulk-mismatch prediction of `0.114512`), the zero disappears:

| case | min \|R\| | monotone in k? |
|---|---|---|
| **piezo** \| passive dielectric | **0.114512** (at k→0) | **yes** |
| **elastic** \| passive dielectric | **8.86×10⁻⁷** (at k = 1237.5) | **no** |

**The corrected finding is more interesting than the artefact:** piezoelectric
coupling *destroys* the elastic SH interference minimum. This is stated as a
structural constraint in §1 and §7, and it is why the recommended problem is
P–SV.

**Severity: numerical/computational error in this audit.** Reported so that the
blueprint is not built on it.

### 4.4 Prior C4 — "the structural space is occupied" — **ACCEPTED but NARROWED**

The prior audit claimed the three-region + Lord–Shulman + dissipation-budget
configuration is published and therefore the novelty claim fails. **Accepted on
the facts, narrowed on the conclusion.**

What is established (§5) is that the *surrounding* configurations are published.
What is **not** established is that the blueprint's specific mathematical
question is answered — the closest competitors compute R/T/absorption for a
particular stack; none of them establishes a *matching criterion*, an
*obstruction*, or a *condition count*. **Combining known effects is not novelty;
but a proof about a known configuration is.** The prior audit conflated these.

### 4.5 Prior verdict "REDESIGN" — **DOWNGRADED to CONDITIONAL REVISION**

The prior audit's four critical findings were C1 (surface elasticity), C2
(impedance), C3 (thermal), C4 (literature). Of these:

- C2 is now **rejected** (§4.1).
- C3 is now **rejected** as stated (§4.2).
- C1 is **accepted but downgraded** to an optional single-paragraph fix (§3.1).
- C4 is **accepted but narrowed** (§4.4).

None of the surviving findings shows the formulation to be "fundamentally
impossible or scientifically indefensible" (Part A). **The research direction
stands.** Per Part G Rule 1, a weakness in one mechanism does not invalidate the
problem.

---

## 5. Literature overlap and remaining research space

**Metadata verified against Crossref. Full texts not read; no claim below
depends on unread content.**

| DOI | Work | Status | Overlap with blueprint |
|---|---|---|---|
| 10.1016/j.compstruct.2023.117564 | Wang, Chen, Yu, Liu, Zhang, Zhou, *Compos. Struct.* **324**, 117564 (2023) | metadata ✔ | **Closest competitor.** Sandwiched TEE layer, Lord–Shulman, R/T **and dissipation** ratios, relaxation time varied, nonlocal integral elasticity. |
| 10.1016/j.tws.2021.108306 | Liu, Yu, Wang, Zhang, Zhang, Zhou, *Thin-Walled Struct.* **168**, 108306 (2021) | metadata ✔ | Piezoelectric plate sandwiched in two solid half-spaces, with nonlocal term. This is the blueprint's model with β = 0, τ_q = 0. |
| 10.1080/17455030.2019.1610198 | Guha, Singh, Das, *WRCM* **31**(4), 660–689 (2021; online 2019) | metadata ✔ exact match | Two piezothermoelastic half-spaces, imperfect interfaces, LS/GL. Genuine closest prior art for Option A. |
| 10.1007/s12648-022-02568-w | Barak, R. Kumar, R. Kumar, **V. Gupta**, *Indian J. Phys.* **97**(8), 2369–2383 (2023) | metadata ✔, ORCID 0000-0002-5238-1963 confirmed | Author's own energy-ratio benchmark, elastic/piezothermoelastic half-spaces. |
| 10.1108/hff-04-2022-0259 | **V. Gupta**, R. Kumar, M. Kumar, V. Pathania, M. S. Barak (2022) | metadata ✔ | Author's own **nonlocal** piezothermoelastic R/T. **Option C is therefore the least novel of the three.** |
| 10.1007/s10483-015-1892-9 | Vashishth & Sukhija, *Appl. Math. Mech.* **36**, 11–36 (2015) | metadata ✔ | Fluid / piezothermoelastic interface R/T. |
| 10.1007/BF00261375 | Gurtin & Murdoch, *ARMA* **57**, 291–323 (1975) | metadata ✔ | Foundation only. |

### What is established vs what remains open

| Element | Established? |
|---|---|
| three-region finite layer, R/T coefficients | **established** (both 2023 and 2021 papers) |
| Lord–Shulman in a sandwiched thermo-electro-elastic layer | **established** (2023) |
| explicit dissipative energy budget as a *result* | **established** (2023) |
| nonlocal piezothermoelastic R/T | **established** (author's own 2022, and 2023) |
| imperfect/spring interface R/T in piezothermoelastic half-spaces | **established** (Guha et al. 2019) |
| **a matching criterion or obstruction for a coupled interface** | **not found** |
| **a condition/rank analysis of the coupled interface system** | **not found** |
| **a proof that coupling destroys an elastic interference minimum** | **not found** (and §4.3 gives the elastic/piezo contrast) |
| **a passivity dichotomy for coupled interfaces (R_E vs T_E order of vanishing)** | **not found** |

**Remaining space.** The open questions are *analytical*, not configurational:
existence and classification of reflection minima, their destruction by
coupling, and the parameter budget. That is exactly where the blueprint's
mathematics is strongest and its literature exposure is lowest.

---

## 6. Comparison of Options A, B and C

Scored per the ten criteria the task requires. Alignment scored against
Dr. Gupta's verified method base (plane-wave R/T at coupled half-spaces,
energy-ratio bookkeeping, generalized thermoelasticity, generalized continua).

| Criterion | **A: two half-spaces + imperfect interface** | **B: three-region finite interphase** | **C: nonlocal extension** |
|---|---|---|---|
| Scientific motivation | Direct: interface physics governs R/T in every real joint | Indirect: matching layer is a device idea, not interface physics | Coupling to a length scale; surface-dominated wave physics |
| Verified literature overlap | 10.1080/17455030.2019.1610198 (imp. interfaces), 10.1007/s12648-022-02568-w (energy) | **Heavy**: 10.1016/j.compstruct.2023.117564, 10.1016/j.tws.2021.108306 | **Heavy and self-overlapping**: 10.1108/hff-04-2022-0259 |
| Physical feasibility | High; no interphase material must be found | **Low**: no bulk solid matches √(Z₁Z₂) for piezo-to-fluid; solid–solid best case 1.7 % (§3, check_04 Part D) | Medium; internal length must be justified, not chosen |
| Governing-equation complexity | n = 8, 8 conditions, one interface | n = 8, **16 conditions**, two interfaces | n > 8 with higher-order BCs; count rises again |
| Interface-condition consistency | Fully consistent; spring/thermal/capacitive all admissible with ≥ 0 coefficients | Consistent but **over-determined** unless the interphase is fully constitutive | Requires higher-order edge conditions; risk of over/under-determination |
| Expected analytical contribution | **High**: existence/classification of minima, destruction by coupling, parameter bounds | Low: reduces to the textbook quarter-wave identity (§3.2), which is *unique* and known | Medium but **already published by the author** |
| Required material parameters | 2 sets: C, e, ε, β, p, K, c_e — all documented for PZT/BaTiO₃ | 3 sets **including a matched interphase**, which does not exist as a real material | 2 sets + one internal length, which is a free parameter |
| Validation benchmarks | Guha et al. (same geometry class), Barak et al. (energy ratios, same author) | **None available**: no competitor reports a three-region piezothermoelastic interphase block | Author's own nonlocal limit |
| Main research risk | Result may be a bounded negative (this audit says: probably is, and that is publishable) | Risk of reproducing a known result in a heavier model | Risk of being judged incremental over the author's own 2022 paper |
| Alignment with Dr. Gupta's methods | **10/10** — exactly his published R/T-at-coupled-half-spaces programme | 7/10 — his programme, but the interphase framing is new territory for him | 9/10 — his nonlocal work, but already done |

**Scored summary (alignment 30 / originality 20 / math richness 20 / physical
feasibility 20 / numerical tractability 10):**

| Option | Alignment | Originality | Math richness | Feasibility | Tractability | **Total** |
|---|---|---|---|---|---|---|
| **A** | 9.6 | 7.0 | 8.0 | 8.5 | 9.0 | **8.52** |
| **B** | 7.0 | 4.0 | 8.5 | 3.0 | 7.0 | **5.95** |
| **C** | 8.7 | 3.0 | 9.0 | 7.0 | 6.5 | **7.06** |

*(Weights: 30/20/20/20/10, scores out of 10, total normalised to 10.)*

Option B fails the physical-feasibility gate on evidence, not preference: the
matched interphase required by the quarter-wave condition does not exist as a
real piezothermoelastic solid for the pairs that motivate the problem. Per
Part G Rule 5 this is a computed result, not an opinion.

---

## 7. Smallest sufficient redesign

Not a redesign of the field. Four changes, in order of importance.

**Change 1 — Replace the headline objective (mandatory).**
Delete quarter-wave matching as the analytical contribution. Replace it with the
*classification* question: **for which material pairs and coupling strengths does
a coupled piezothermoelastic interface admit a reflection minimum, and what
destroys it?** §3.2 supplies the complete answer for the scalar case; §4.3
supplies the first nontrivial coupled answer (piezoelectric coupling destroys the
elastic SH minimum). This is a smaller change than it looks: the mathematics
already in the blueprint (§§4–7) is exactly what is needed.

**Change 2 — Fix the condition count and state it per model (mandatory).**
Insert the n = 8 derivation and the 8/16 count of §3.3 into §5, separately for
B0, B1, B2. State that seven conditions is rank deficient and reconcile with the
source.

**Change 3 — Restate the energy section (mandatory).**
Remove the displayed `R_E + T_E + A_E = 1` before the balance is derived. Add
the passivity dichotomy and the order-of-vanishing table of §3.4. State that
reporting R_E alone is not evidence of anti-reflection.

**Change 4 — Adopt Option A and drop the finite interphase (scope).**
This is a scope reduction, not a mathematical correction, and the supervising
team may overrule it. The evidence for it: the matched interphase does not exist
for the motivating pairs (§6), the three-region count doubles the interface
system (§3.3), and the closest competitors already occupy the three-region
configuration (§5). Option B remains defensible if the team accepts an
illustrative-only interphase; that concession must then be stated in §3.

**Retain the temperature field.** Per §4.2, add the constraint relation
`v_T/c_L = √(2π f_on τ_min)` and state that τ_q must satisfy `τ_q ≥ κ/c_L²`.

---

## 8. Recommended final research title

> **Existence and Destruction of Reflection Minima at Piezothermoelastic
> Interfaces with Generalized Imperfect Contact: A Coupled-Field Analytical
> Study**

Sub-title options (pick one, not both):

- *… with a Passivity Bound on Reflection Suppression*
- *… : Interference, Coupling and the Limits of Impedance Matching*

The title names a *question* (existence and destruction), states the method
(analytical), and claims no new physical mechanism. It is honest at the
blueprint's current stage of verification.

---

## 9. Exact mathematical tasks required

Ordered so that each task's failure stops the chain early.

**T1 — Establish the P–SV reference problem.** Two piezothermoelastic
half-spaces, 6mm poled in-plane, welded interface. Derive the secular determinant,
classify all 8 roots by branch, and verify the elastic limit against a
Rayleigh-angle benchmark. *Deliverable: determinant, branch table, elastic-limit
error.*

**T2 — Solve the elastic P–SV three... two-half-space scattering and locate every
reflection minimum.** Report `|R_pp|`, `|R_ps|` and total reflected power over
angle. *Expected: minima exist (the elastic SH case has one at |R| ≈ 8.9×10⁻⁷;
P–SV should be richer). Deliverable: minimum map.*

**T3 — Repeat T2 with piezoelectric coupling on.** *This is the core task.*
Determine whether each minimum from T2 survives, moves, or is destroyed.
*Expected on the SH evidence: destroyed or shifted by O(k_piezo/k). Deliverable:
the destruction theorem or counterexample.*

**T4 — Prove or refute the parameter budget.** For the full n = 8 junction,
compute the Jacobian of the 8 conditions with respect to the interface parameters
at each candidate minimum. Establish rank. *Deliverable: either a local existence
proof (Implicit Function Theorem) or a rigorous obstruction.*

**T5 — Add the thermal field and the LS constraint.** Insert
`v_T/c_L = √(2π f_on τ_min)` and the causality condition `τ_q ≥ κ/c_L²`. Report
the O(ε) shift of each minimum with ε from three verified parameter sets.
*Deliverable: sensitivity table with ε quoted per material.*

**T6 — Establish the passivity dichotomy for the coupled case.** Prove
`R_E + T_E + A_E = 1` from the coupled flux, then establish and prove the order
of vanishing of R_E (second) vs A_E (first) in the interface loss parameter.
*Deliverable: proposition + proof.*

**T7 (optional) — Reduce to a closed-form asymptotic matching relation** in a
tractable subcase (e.g. weak piezoelectric coupling), with a stated domain of
validity.

**T8 — Numerical parametric study** of the minimum's location and depth versus
(a) angle, (b) frequency, (c) coupling strengths, (d) interface parameters.
*Prerequisite: T2–T4 pass.*

---

## 10. Validation and reproducibility checklist

Every gate must produce a number, not a verdict.

| # | Gate | Pass criterion | Status in this audit |
|---|---|---|---|
| V1 | Identical media ⇒ R = 0 exactly | \|R\| < 1e-14 | **PASS** (1.5e-16, 2.9e-17) |
| V2 | Near-normal incidence ⇒ bulk impedance mismatch | rel. error < 1e-6 | **PASS** (0.114512 vs 0.114512) |
| V3 | Elastic limit recovers validated solver | rel. error < 1e-8 | **PASS** (prior audit V1: 7.0e-12) |
| V4 | Energy balance lossless | \|R_E+T_E−1\| < 1e-12 | **PASS** (≤ 7e-16) |
| V5 | Energy balance lossy | \|R_E+T_E+A_E−1\| < 1e-12 | **PASS** (≤ 2e-15) |
| V6 | Closed form vs direct evaluation | rel. error < 1e-12 | **PASS** (1.8e-15) |
| V7 | LS with τ→0 recovers Fourier | rel. error < 1e-10 | **PASS** (0.0) |
| V8 | LS speed from dispersion | rel. error < 1e-6 | **PASS** (1.3e-9) |
| V9 | LS causality τ ≥ κ/c_L² | all sets causal | **PASS** |
| V10 | Loss sign convention | decay in +z | **PASS** (Im k > 0 ⇔ c = c₀(1−iη)) |
| V11 | Mesh/step-size convergence | stable over 5 step sizes | required, not run here |
| V12 | Reduction limits (β→0, e→0, K→∞) | each recovers its stated limit | required |

**Two conventions that must be stated explicitly in any paper arising from
this** — both caused errors in this audit:

1. **Loss sign.** With `e^{i(kz−ωt)}`, `c = c₀(1 − iη)` decays; `c = c₀(1 + iη)`
   **grows**. Conflating them turns absorption into gain and makes `A_E < 0`.
2. **Continuity impedance vs flux impedance.** For a lossy layer, continuity uses
   the complex `ρc` while the power flux uses the **real** `ρ Re(c)`. Using one
   for both either loses the loss entirely or produces `A_E < 0`.

---

## 11. Research risks and unresolved questions

**R1 — The main result may be negative.** The evidence points to reflection
minima being *destroyed* by coupling rather than created (§4.3). This is
publishable and the blueprint's H0 already anticipates it, but it must be framed
as a result from the outset, not discovered late.

**R2 — The interface law set may be incomplete in the chosen reduction.** In SH
with a tangential spring, 4 amplitudes face only 3 independent conditions
(§3.3). The P–SV system must be checked for the same defect before any solver is
trusted.

**R3 — τ_q is not constrained by any measurement in the intended band.**
No ultrasonic experiment fixes τ_q for PZT; the value enters only through
`ω τ_q ≪ 1`. Every statement about LS in the paper must be phrased as a
*constraint* (§4.2), never as a *prediction*.

**R4 — Material data for a piezothermoelastic interphase do not exist.** Not a
modelling difficulty but a data gap. Option B therefore requires the illustrative
concession.

**R5 — k = 0 is singular in the SH gauge.** Normal incidence must be taken as a
limit, not substituted. This will bite any code that is not written to handle it.

**Unresolved and needing the full texts:**

- **U1.** Does Guha et al. (2019) use 8 interface conditions or 7, and if 7,
  which field is fixed by an extra assumption? Settles whether §3.3 is a
  correction to the blueprint or to the prior art.
- **U2.** Does the 2023 *Compos. Struct.* paper contain any matching criterion,
  or only computed R/T/dissipation curves? Determines whether §5's
  "not found" is a genuine gap or an artefact of not having read it.
- **U3.** Is the elastic P–SV reflection minimum (§9, T2) destroyed by coupling
  for *all* material pairs, or only in the anti-plane case? This is the
  central open question of the recommended problem and cannot be answered
  without T3.

---

## 12. Recommendation for supervisor approval

**Recommendation: APPROVE the direction, CONDITIONAL REVISION of the blueprint.**

The research problem is scientifically defensible, aligned with the author's
verified method base, mathematically substantial, and free of any requirement for
experimental fabrication, metamaterials, electronics, machine learning or
finite-element-only workflow. It should proceed.

It should proceed **after** the four changes in §7, of which three are mandatory
corrections and one is a scope decision for the supervising team.

It must **not** proceed on the basis of the previous audit's verdict. Two of that
audit's three central mathematical objections are wrong (§4.1, §4.2), and one of
its own intermediate findings was a computational artefact (§4.3). The blueprint
deserves a controlled revision, not a rewrite.

---

## 13. Decision table

| # | Original claim | Mathematical evidence | Counterexample or supporting proof | Verified DOI references | Status | Required blueprint correction |
|---|---|---|---|---|---|---|
| D1 | §9.1: `Z_m² = Z₁Z₂` with `k_mh=(2n+1)π/2` is the matching benchmark | Part Q1: closed form `Re Z_in = Z₂Z_m²/(Z₂²sin²θ+Z_m²cos²θ)`, `Im Z_in = Z_m(Z_m²−Z₂²)sinθcosθ/(…)`. Closed form vs direct: 1.8e-15. `\|R\|` = 9.0e-17 at the quarter-wave point. | **Supporting proof:** `Im Z_in = 0` iff θ = mπ/2, where `Re Z_in ∈ {Z₂, Z_m²/Z₂}`. So for real lossless media the zero-reflection set is exactly `{Z₂, Z_m²/Z₂}` ⇒ the quarter-wave condition is **necessary and sufficient**, i.e. a complete classification of a known result. | — (classical) | **ACCEPTED as correct mathematics; REJECTED as a novel contribution** | State the uniqueness result; remove §9.1 from the contribution list; demote to a verification gate only. |
| D2 | Prior audit C2: with loss, only gain can close the reflection floor | Part Q1 sweep: `Re Z_in < Z₁` in **198/366** configurations. Fixed-material family: floor 1.8e-7 → 3.4e-3 with loss. | **Counterexample:** with the interphase impedance free, `Re Z_in = Z₁` is attainable with loss present, and one lossless reactance then cancels `Im Z_in`, giving exact zero reflection *with absorption*. | — | **REJECTED as a general claim; CONDITIONAL for the fixed-material family** | None for the blueprint. The prior audit must be amended. Any obstruction claim must name its parameter family. |
| D3 | Prior audit C3: the thermal field is inert; delete it | Part Q2: exact LS dispersion `k² = ω(i+ωτ_q)ρc_e/K`. S3 (τ→0 → Fourier): 0.0 error. S4 (`v_T` from dispersion): 1.3e-9. ε: PZT 1.20e-4, BaTiO₃ 2.45e-3, Al 7.10e-3. | **Supporting proof of a different statement:** for `ωτ_q ≫ 1` a real root exists at `v_T = √(K/ρc_eτ_q)`. Causality `v_T ≤ c_L` ⇒ `τ_q ≥ κ/c_L²`, hence `v_T/c_L = √(2π f_on τ_min)` — speed and onset are **not** independent. | 10.1016/0022-5096(67)90024-5 (Lord–Shulman, as cited in verified metadata) | **REJECTED as a general claim; the *Fourier* sub-case is ACCEPTED** | Retain the temperature field. Add the constraint relation and the causality condition. Present attenuation as O(ε), frequency-dependent. |
| D4 | §5: the interface condition count can be carried over from two half-spaces | Part: degree n = 8 (4 from u_x,u_z; 2 from φ; 2 from θ). Two half-spaces: 8 unknowns, C = 8. Three regions: 16 unknowns, C = 16. | **Supporting derivation.** Additionally: in SH with a tangential spring, 4 amplitudes face only 3 independent conditions — an incomplete law set produces a spurious null vector, the same failure mode that produced §4.3's artefact. | — | **ACCEPTED — mandatory** | Insert the count per model (B0/B1/B2) in §5. Reconcile the "seven conditions" of the cited prior art. |
| D5 | §8: `R_E + T_E + A_E = 1` may be displayed before the balance is derived | Part Q6: `A_E = 0` iff lossless. Measured: `R_E/η²` = 0.89, 0.67, 0.64 vs `A_E/η` = 3.54, 3.53, 3.42. Balance residual ≤ 2e-15. | **Supporting proof of the dichotomy:** `R_E → 0` with `A_E > 0` is allowed; `R_E = 0` with `T_E = 1` forces `A_E = 0` (a losslessness, not a matching, statement). Reflected power is **second order** in loss, absorbed power **first order**. | — | **ACCEPTED — mandatory** | Move the identity after the derivation. Add the order-of-vanishing table. Require `R_E`, `T_E`, `A_E` to be reported separately in H1. |
| D6 | §6.3: surface elasticity "gated" by `η_s ~ k_∥μ_s/μ_bulk` | Part C: `η_s = 2πℓ_s/λ`, `ℓ_s = μ_s/μ_bulk` ≈ 1e-10 m. | **Supporting numbers:** η_s = 1.2e-6 at λ = 1 mm; = 1.2e-3 at λ = 1 μm; unity requires λ ≈ 1.2 nm. | 10.1007/BF00261375 | **ACCEPTED but DOWNGRADED** — optional improvement, not an error | Replace §6.3 with the table; delete study S6; keep the field as a stated-limitation note. |
| D7 | Prior audit C4: the three-region + LS + dissipation configuration is published, so novelty fails | Crossref metadata for the 2023 and 2021 papers. | **Supporting evidence of overlap; insufficient for the conclusion.** Overlap is in *configurations* (sandwiched layer, LS, dissipation ratios, nonlocal). No competitor establishes a *matching criterion*, an *obstruction*, or a *condition count*. | 10.1016/j.compstruct.2023.117564; 10.1016/j.tws.2021.108306; 10.1108/hff-04-2022-0259 | **ACCEPTED on facts, NARROWED on conclusion** | Rewrite §2 around these three DOIs. Delete §11 item 7 (a two-half-space benchmark cannot validate an interphase block). Re-target novelty to the analytical question. |
| D8 | Prior audit verdict: REDESIGN | This audit: C2 and C3 rejected; C1 downgraded; C4 narrowed. No surviving finding shows the formulation impossible or indefensible (Part A). | **Supporting:** the formulation's equations, field set, interface hierarchy and reduction limits are all sound (§2). | — | **REJECTED — downgraded to CONDITIONAL REVISION** | None. Proceed with the four changes of §7. |
| D9 | (This audit, retracted) an exact reflection zero exists at `k = 2890.43` in an anti-plane piezo junction | Two errors: wrong sign in `κ²`; `κ` used instead of `dψ/dz` in the `σ_z3`/`D_z` rows. | **Counterexample to itself:** corrected solver gives identical-media `\|R\|` = 1.5e-16 and near-normal `\|R\|` = 0.114512 = bulk mismatch. `\|R\|` is then strictly **monotone** in k; min = 0.114512 at k→0. | — | **RETRACTED** | Do not use. Recorded so the blueprint is not built on it. |
| D10 | (New) an elastic SH junction has an interference minimum that piezoelectric coupling removes | Part Q5: piezo `\|R\|` min = 0.114512 at k→0, monotone = **True**. Elastic `\|R\|` min = **8.86e-7** at k = 1237.5, monotone = **False**. | **Supporting evidence** across the propagating window k ∈ (0, k_crit), 40 001 points, with exact self-tests. Thermal field decouples exactly in SH (β₁₃ = β₂₃ = 0), so SH cannot test thermal coupling. | — | **ACCEPTED as a structural constraint** | Any reflection-minimum claim must be **P–SV** (in-plane), not anti-plane. State this in §3 and make it the reason for the chosen geometry. |
| D11 | §1 H1: an interphase gives "lower reflected power" while "preserving higher transmitted power" | Part Q6: the two objectives diverge at first order in loss. | **Supporting proof:** `R_E = O(η²)`, `A_E = O(η)`, so a reflection-only optimiser can certify an attenuating window as an anti-reflection success. | — | **ACCEPTED — mandatory** | Split H1 into two falsifiable statements (reflection and transmission) with separate numeric targets. |
| D12 | §3: "material-pair gate (not yet fixed)" is a data-availability question | Prior audit Part D: required `Z_m = √(Z₁Z₂)`. | **Supporting numbers:** PZT/water requires 6.67 MRayl, nearest solid epoxy 2.95 (56 % off); PZT/steel requires 37.04, brass 37.66 (1.7 % off). No bulk *piezothermoelastic* solid satisfies either. | — | **ACCEPTED — realisability, not data** | Restate as a realisability gate with the numbers. If Option B is retained, state the illustrative concession explicitly in §3. |
| D13 | §14: nonlocal/strain-gradient deferred to a second paper | Author's own 2022 nonlocal piezothermoelastic R/T exists. | **Evidence:** 10.1108/hff-04-2022-0259. | 10.1108/hff-04-2022-0259 | **ACCEPTED — optional** | State plainly that the first paper uses the local limit deliberately, so the interface result is not confounded with a length scale. |
| D14 | §12: 7 figures / 7 studies is a paper | — | — | — | **Optional improvement** | Reduce to S1–S3, S5 plus the obstruction; S4 becomes a consistency table, S6 a single number. |

---

**RECOMMENDED MODEL:** **A** — two piezothermoelastic half-spaces with a
generalized imperfect interface (mechanical spring, thermal contact
conductance, capacitive electrical contact).

**SCIENTIFIC JUSTIFICATION:** Option A is the only one of the three that
scores above 8/10 (8.52) and the only one whose physical feasibility survives
quantitative test. Option B (5.95) fails the realisability gate: its required
matched interphase does not exist as a bulk piezothermoelastic solid for the
pairs that motivate the problem (PZT/water needs 6.67 MRayl; the nearest solid is
2.95), its interface system doubles to 16 conditions, and its configuration is
already occupied by the 2023 and 2021 competitor papers. Option C (7.06) is
scientifically sound but overlaps the author's own published nonlocal work
(10.1108/hff-04-2022-0259), making it the least original choice. Option A is
also where the surviving open questions live: existence and classification of
reflection minima, their destruction by piezoelectric coupling, the rank of the
coupled interface system, and the passivity bound on reflection suppression. It
preserves every element the task requires — coupled piezothermoelastic PDEs,
normal-mode analysis, characteristic roots, interface determinants, R/T
coefficients, energy ratios and parametric study — with no experimental
metamaterial, electronic, machine-learning or FE-only dependency. The strongest
available evidence about what the answer will be (piezoelectric coupling
*destroys* the elastic SH minimum, §4.3/D10) arrives from a two-half-space model,
which is precisely Option A's geometry.

**REQUIRED CHANGES:**
1. Replace quarter-wave matching as the headline contribution with the
   existence/destruction-of-minima question (§7 Change 1).
2. Insert the n = 8 condition count, per model, into §5 (§3.3).
3. Restate §8: move the energy identity after its derivation, add the passivity
   dichotomy and the order-of-vanishing table (§3.4).
4. Split H1 into separate falsifiable reflection and transmission objectives
   with numeric targets (D11).
5. Adopt the two-half-space geometry; if Option B is retained, state the
   illustrative-interphase concession in §3 (§6, D12).
6. Retain the temperature field and add the LS constraint
   `v_T/c_L = √(2π f_on τ_min)` with causality `τ_q ≥ κ/c_L²` (§4.2).
7. Rewrite §2 around the three verified competitor DOIs and delete §11 item 7
   (§5).

**OPTIONAL CHANGES:**
1. Replace §6.3 with the surface-elasticity table; delete study S6 (D6).
2. Reduce §12's figure/study budget to S1–S3, S5 plus the obstruction (D14).
3. State in §14 that the local limit is deliberate, not a deferral (D13).
4. Add the two convention statements of §10 (loss sign; continuity vs flux
   impedance) to the blueprint's reproducibility section.
5. Note the k = 0 singularity of the SH gauge so no code substitutes it (R5).

**RESEARCH ALIGNMENT:** **9.0 / 10.** The recommended problem uses only
applied mathematics, continuum wave dynamics, coupled-field R/T at interfaces,
generalized thermoelasticity and normal-mode analysis — the author's entire
verified method base (plane-wave R/T at coupled half-spaces; energy-ratio
bookkeeping; generalized thermoelasticity). The 1.0 deduction is for the fact
that the *headline objective* changes from the author's prior "compute the R/T
ratios" pattern to "prove a structural property of them", which is adjacent but
not identical to his published programme.

**MATHEMATICAL FEASIBILITY:** **Conditional.**
Feasible: the determinant, the root classification, the passivity dichotomy, the
LS constraint, and the elastic/piezo contrast are all established or provable
within the intended framework, and every one of the eight self-checks in §10
passes. Conditional on: T3 (whether the P–SV minimum is destroyed by coupling),
which is the paper's central result and which no amount of anti-plane analysis
can settle; and on T4 (the rank of the coupled interface Jacobian), which
determines whether any existence statement can be made rigorously. If T3 yields
"the minimum survives" the paper is an existence result; if it yields "the
minimum is destroyed" the paper is an obstruction result. **Both outcomes are
publishable.** That is why the feasibility is Conditional rather than Pass or
Fail: the pass/fail depends on which branch of a dichotomy is true, and the
dichotomy itself is the contribution.

**NOVELTY STATUS:** **Plausible — not Verified.**
Verified as published: the three-region sandwiched configuration, Lord–Shulman
in a thermo-electro-elastic layer, explicit dissipation ratios, nonlocal
piezothermoelastic R/T, imperfect-interface piezothermoelastic R/T. Not found in
the literature searched: a matching *criterion or obstruction* for a coupled
interface; a rank/condition analysis of the coupled interface system; a proof
that coupling destroys an elastic interference minimum; a passivity dichotomy in
the order-of-vanishing form. Per Part G Rule 3, "not found" is **not** "novel":
the three competitor full texts have not been read (U2), and novelty can only be
asserted after that reading. The status is therefore Plausible at best and the
claim must not be stronger than that in any outward-facing document.

**APPROVAL STATUS:** **AWAITING INDEPENDENT SUPERVISOR DECISION.**

---

*Audit performed 4 November 2026. Every number is reproducible via
`python3 code/check_05_dispute_resolution.py`. Eight enforced self-checks:
identical-media limit (1.5e-16), near-normal bulk-mismatch limit (0.114512 vs
0.114512), closed-form vs direct impedance (1.8e-15), quarter-wave point
(9.0e-17), lossless energy balance (≤ 7e-16), lossy energy balance (≤ 2e-15),
LS Fourier limit (0.0), LS thermal speed from dispersion (1.3e-9), LS causality
(all sets pass). Three errors made during this audit — a wrong loss sign, a wrong
impedance used for continuity, and a wrong sign in the anti-plane dispersion
relation — were each caught by a self-check and are recorded in §4.3 and §10
rather than removed.*

**STOP. Awaiting independent supervisor decision.**
