# ARENA_V2_CRITICAL_AUDIT.md

**Document audited:** `FINAL RESEARCH BLUEPRINT — V2.0`,
*"Analytical Anti-Reflection and Energy Transmission of Coupled Waves in Layered
Piezothermoelastic Media with Generalized Interface Conditions"*, 4 November 2026.

**Referee role:** hostile mathematical referee. Not a co-author. No manuscript
text has been produced. No numerical result below is presented as a result of
the proposed paper; every number is either a property of the document's own
arithmetic, or the outcome of a *derivation completed in the audit* to test a
specific claim.

**Reproducibility:** all checkable numbers come from
`code/check_04_blueprint_audit.py` → `data/check_04_blueprint_audit.json`.
The script carries three mandatory self-checks (uniform-medium limit,
lossless energy balance, matched-layer quarter-wave limit) which are printed on
every run. Claim-by-claim evidence is in `ARENA_V2_CLAIM_EVIDENCE.csv`.

---

## 0. Verdict

**REJECT as written. REDESIGN — do not revise.**

The document is unusually disciplined. Its V2.0 retraction of Gurtin–Murdoch as
the principal mechanism, its separation of anti-reflection from negative
reflection, its null-hypothesis framing (H0) and its explicit refusal to let
parameter counting settle a feasibility question are all correct and were the
right calls. Those are genuine strengths and they are why the redesign below is
still possible.

It nevertheless fails, for four independent reasons, each of which alone is
sufficient:

1. **The declared headline contribution is a textbook identity**, and the audit
   reproduces it to 6×10⁻¹⁰ in the configuration's own elastic reduction. It
   cannot be a contribution.
2. **The matched condition survives only at normal incidence.** Off normal, the
   best achievable reflection is 1.9×10⁻², and the optimum leaves λ/4. The
   blueprint's own performance benchmark is unreachable for the reason the
   blueprint itself names but does not act on.
3. **The temperature field cannot carry the load the document assigns to it.**
   Its coupling parameter is 3×10⁻³ to 3×10⁻² and its effect on a travelling
   wave is an attenuation of 10⁻¹⁰ Np/mm at 1 MHz. It adds no branch, no
   absorption, and no measurable power signature. "Three-region
   *piezothermoelastic*" therefore advertises a fourth field that a referee will
   immediately test and find inert.
4. **The structural space is occupied.** A 2023 *Composite Structures* paper
   already computes reflection, transmission **and dissipation ratios** for a
   sandwiched, functionally graded, **thermo-electro-elastic** layer under
   **Lord–Shulman**. The single element the blueprint treats as its novelty —
   the third region plus generalized thermoelasticity plus an explicit
   dissipative energy budget — is a published configuration.

A point of process: the document's own literature list contains three items and
cites none of them in the body. The competitor above is absent. That is not a
citation-style complaint; it is the specific reason finding 4 was not caught
before the pivot.

---

## 1. What I could and could not verify

Stated explicitly, per the document's own standards.

**Verified against registries (Crossref).** The three cited works exist with the
titles, authors, journals, volumes and pages given:

| DOI | As cited in §2 | Verified |
|---|---|---|
| 10.1007/BF00261375 | Gurtin & Murdoch, *ARMA* **57**, 291–323 (1975) | ✅ title/cite match |
| 10.1080/17455030.2019.1610198 | Guha, Singh, Das, *WRCM* **31**, 660–689 (2021; online 2019) | ✅ exact; print vol 31(4) 2021, online 1 May 2019 |
| 10.1007/s12648-022-02568-w | Barak, R. Kumar, R. Kumar, Gupta, *Indian J. Phys.* (2023) | ✅ exact; **97**(8), 2369–2383, July 2023, ORCID 0000-0002-5238-1963 confirmed |

The Barak citation is a legitimate author-aligned benchmark. Its full bibliographic
data should be inserted now rather than "before submission".

**Located, not yet read in full, and material to the novelty claim:**

- `10.1016/j.compstruct.2023.117564` — Wang, Chen, Yu, Liu, Zhang, Zhou,
  *Composite Structures* **324**, 117564 (2023). *Reflection, transmission, and
  dissipation of plane waves in sandwiched functionally graded
  thermo-electro-elastic nanoplates via nonlocal integral elasticity theory.*
  LS theory, sandwiched TEE layer, R/T **and dissipation** energy ratios,
  relaxation time and nonlocal effect studied. **This is the paper the blueprint
  must be positioned against, and it is not in the document.**
- `10.1016/j.tws.2021.108306` — Liu, Yu, Wang, Zhang, Zhang, Zhou,
  *Thin-Walled Structures* **168**, 108306 (2021). Nonlocal piezoelectric plate
  **sandwiched in two solid half-spaces**. This is the β = 0, τ_q = 0 reduction
  of the proposed model, published with the nonlocal extension already in place.
- `10.1007/s10483-015-1892-9` — Vashishth & Sukhija, fluid/piezothermoelastic
  interface. Already cited *by* Barak et al.; relevant prior art.
- `10.1108/hff-04-2022-0259` — Gupta, Kumar, Kumar, Pathania, Barak, nonlocal
  piezothermoelastic R/T with an ideal fluid. **The author's own nearest work,
  and it already contains the nonlocal piezothermoelastic scattering machinery
  the blueprint defers to "a separate second paper".**

**Not verifiable in this environment.** Full texts of the three competitor
papers. I therefore do **not** assert that the proposed exact configuration is
published; I assert that the *neighbouring* configurations and, critically, the
*three-region + LS + dissipation-ratio* combination are published, and that the
blueprint contains no analysis of them. That is sufficient to reject the
novelty claim as unestablished, which is all a referee needs.

**Material constants** used in the audit are handbook-level and illustrative.
Every conclusion drawn from them is a separation of one decade or more, so none
of it is sensitive to ±50 % changes. Constants that would matter at the margin
are marked as such.

---

## 2. Findings

Severity: **C** = critical (blocks acceptance) · **M** = major · **m** = minor.

---

### C-1. The headline contribution is a textbook identity, and it only survives where the problem is trivial

**Location:** §9.1 (the `Z_m² = Z_1 Z_2` subsection, marked "compulsory"),
§9.2, §12 item 5.

**The defect.** §9.1 states the quarter-wave identity and correctly labels it a
benchmark. But it is then the *only* closed form the paper is designed to produce:
§9.2 offers three alternatives — (a) exact condition in a subcase, (b) rigorous
obstruction, (c) perturbative relation — and §12 lists "at least one analytic
reduction" as *mandatory*. Alternative (a) reduces to §9.1. So the mandatory
analytical content of the paper is a formula from a first-year course in
ultrasonics, and §9.1 itself concedes it "is not the solution to the full
coupled piezothermoelastic problem".

**Evidence.** Audit Part E. Lossless scalar layer, `Z_m = √(Z_1 Z_2)`,
`h = λ/4`:

```
min |R| = 1.9e-06     h/λ = 0.7500
```

the residual being numerical. Audit Part F5, taking a *mathematically ideal*
isotropic interphase with `Z_m = √(Z_1 Z_2)` and `μ_m = √(μ_1 μ_2)` inside the
full elastic **P–SV** three-region model:

```
θ = 0° : min R_E = 6.28e-10  at h/λ = 0.2500   T_E = 1.000000
```

Exact, to ten digits. The blueprint's own §12 item 6 requires exactly this
comparison and it passes — which is precisely the problem. A verification gate
that the trivial answer passes is not a research result.

**Required correction.** Do not build the analytical contribution on §9.1. Either
(i) derive a genuinely coupled obstruction (see C-3 and M-2), or (ii) delete §9.1
from the contribution statement and reclassify it as a numerical verification
gate only. The present document does neither and describes it as "compulsory".

---

### C-2. The matching condition fails off normal incidence, and the optimum is not λ/4

**Location:** §1 H1 ("at selected frequencies **and incidence angles**"),
§9.1 (`k_m h = (2n+1)π/2`), §10 S5 (angular tolerance).

**The defect.** §9.1 derives `Z_m² = Z_1Z_2` and `k_m h = (2n+1)π/2` and states
it holds "under its own assumptions" without ever stating that the assumption
that binds is **normal incidence**. H1 promises performance "at selected
frequencies and incidence angles". The identity is not an angular condition. It
degenerates the moment the incident wave is not normal, because
(a) `Z` becomes `Z/cos θ` (or `Z cos θ`) and is angle-dependent per branch;
(b) the two shear branches are reactivated, so `R` stops being one complex
number in one real variable.

**Evidence (Audit Part F5, same ideal matched interphase, full P–SV, 1 MHz):**

| θ | min R_E | optimal h/λ | \|R_pp\| | \|R_ps\| | T_E |
|---|---|---|---|---|---|
| 0° | 6.28e-10 | 0.2500 | 2.5e-05 | 0 | 1.000000 |
| 10° | 6.15e-03 | 0.2461 | 2.1e-02 | 5.3e-02 | 0.993853 |
| 20° | 1.88e-02 | 0.8255 | 7.6e-02 | 7.8e-02 | 0.981228 |
| 30° | 9.19e-03 | 0.8978 | 8.0e-02 | 3.5e-02 | 0.990813 |
| 40° | 2.69e-05 | 0.9608 | 4.9e-03 | 1.2e-03 | 0.999973 |

Energy residual at every row: ≤ 7×10⁻¹⁶ (self-check).

So even with a mathematically ideal interphase, the floor is **1–2 % reflected
power** across the intended angular range, and the optimal thickness jumps off
λ/4 by a factor of 3.3–3.9. H1's "selected frequencies and incidence angles" is
doing the work of hiding a factor-10⁷ degradation from normal incidence.

**Required correction.** §9.1 must state: *the identity is a normal-incidence
condition; off normal it is neither necessary nor sufficient, and the optimal
interphase thickness departs from λ/4 by a factor of order unity.* H1 must be
rewritten to name the angular window in which the claim is made, and must
specify whether the tolerance target is R_E < 10⁻², 10⁻³, or 10⁻⁶, because those
are different papers.

---

### C-3. The temperature field is inert at the intended scale; "piezothermoelastic" is a fourth field in name only

**Location:** §1, §4 ("Exactly four scalar primary unknowns"), §8, §10 S4
("relaxation time vs Fourier limit → absorption vs transmission"), §9's
`R_E + T_E + A_E = 1`.

**The defect.** Three separate claims in the document require the temperature
field to produce a *measurable power signature*:

- S4 asks for "absorption vs transmission" as a function of τ_q;
- §8 builds an absorption term A_E into the budget;
- §1 calls the medium "piezothermoelastic", i.e. four fields.

None of this is examined quantitatively anywhere in the document.

**Evidence.**

*Part A.* The thermal root has `|k_th| = √(ω/2κ)`, i.e. a diffusion branch. Its
ratio to the acoustic wavenumber:

| material | k_th/k_ac at 1 MHz | at 1 GHz | crossover |
|---|---|---|---|
| PZT-5H | 1.47×10³ | 46.6 | 2.2×10¹² Hz |
| BaTiO₃ | — | — | 2.8×10¹² Hz |
| Aluminium | 194 | 6.1 | 3.8×10¹⁰ Hz |

The ratio falls as ω^(−1/2) and reaches unity only where the Fourier model has
already failed. **The thermal root is never a propagating far-field channel in
any addressable band.** So there is no fourth *branch* to reflect or transmit.

*Part H.* Thermoelastic coupling parameter `ε = β²T₀/(ρ c_v M_L)` and the
resulting longitudinal attenuation:

| material | ε | α_th at 1 MHz | at 10 MHz |
|---|---|---|---|
| PZT-5H | 3.2×10⁻³ | 5.8×10⁻¹⁰ Np/mm | 5.8×10⁻⁸ Np/mm |
| BaTiO₃ | 8.2×10⁻³ | 8.8×10⁻¹⁰ Np/mm | 8.8×10⁻⁸ Np/mm |
| Steel | 9.9×10⁻³ | 1.2×10⁻⁸ Np/mm | 1.2×10⁻⁶ Np/mm |

*(Dimension check, recorded because it caught an error in this audit: `ε` must be
normalised by the longitudinal modulus `M_L = ρc_L²`, not by `c_L²`. The first
version of the script omitted a factor ρ and returned ε ≈ 24. The wrong value is
retained in the script's docstring so the check is reproducible.)*

So: **ε ~ 10⁻³–10⁻²** (temperature does shift the mechanical branch at the
few-tenths-of-a-percent level — real but small), and **α_th ~ 10⁻¹⁰–10⁻⁶ Np/mm**
(propagation over 10 wavelengths loses 10⁻⁹ of its amplitude). The document's
absorption term A_E is therefore zero to roughly nine decimal places. The
"absorption-dominated solution" against which H1 promises to compare does not
exist in the model.

*Part B.* `ωτ_q ≤ 6.3×10⁻⁴` for τ_q = 10⁻¹¹ s at 10 MHz. Lord–Shulman and
Fourier are numerically indistinguishable in the band the model can describe.
S4 will produce one curve, not a family.

**Required correction.** Three options, in order of preference:

1. **Demote the thermal field explicitly.** State in §1 that the fourth unknown
   is a *quasi-static potential* that couples the mechanical branch at O(ε) but
   contributes **no propagating mode and no O(1) dissipation** in the intended
   band, and cite Parts A/B/H of this audit as the screening calculation. Then
   delete A_E from the headline budget and delete S4 as a *results* study,
   retaining it as a consistency check. This is honest and costs one paragraph.
2. **Move to a theory where it is not inert.** A hyperbolic heat transport that
   keeps a propagating thermal mode requires either very high frequency or an
   explicitly non-Fourier state; state the frequency at which the thermal branch
   propagates and accept that this is a different physical regime.
3. **Delete "piezothermoelastic"** from the title and treat the model as
   piezoelectric with a thermal perturbation. This contradicts the author's
   identity and is not recommended.

Do **not** leave the current framing: a referee who computes A_E and finds
10⁻⁹ will read §8 as an inflation of the field count.

---

### C-4. The three-region + Lord–Shulman + dissipation-ratio configuration is published

**Location:** §2 (whole section), §11 item 7, §12.

**The defect.** §2 lists three works and states "Record exact DOI, title,
physical model, analytical result and overlap. **Do not claim a verified gap
until this comparison is completed.**" The comparison has not been done, and the
gap it was meant to establish does not survive it.

**Evidence.** `10.1016/j.compstruct.2023.117564` (Wang et al., *Compos. Struct.*
324, 2023): sandwiched layer, thermo-**electro**-elastic, **Lord–Shulman**,
reflection **and transmission and dissipation** ratios, relaxation time varied,
nonlocal integral elasticity included. The proposed paper's distinguishing
elements are therefore:

| Proposed element | Status against 2023 paper |
|---|---|
| three regions, one finite layer | **published** (there) |
| Lord–Shulman generalized thermoelasticity | **published** (there) |
| explicit dissipative energy budget A_E | **published** (there, and as a *result*) |
| nonlocal extension | **published** (there, as the headline) |
| piezoelectric coupling (not just electroelastic) | **open** (there: TEE, not piezo) |
| quarter-wave / impedance-matching **optimisation** of the interphase | **open** — this is the only surviving element |
| anti-reflection as an objective | **open** |

And `10.1016/j.tws.2021.108306` covers the β = 0, τ_q = 0 reduction (nonlocal
piezoelectric plate between two half-spaces) already with the nonlocal term.

So the surviving novelty is *one* element: **treating interphase matching as an
optimisation problem with a power metric**. That is a legitimate paper. It is not
the paper §1 describes, and it requires C-1 to be resolved first, because C-1
says the optimisation has a textbook closed form at normal incidence and a 1–2 %
floor off it.

**Required correction.** §2 must be rewritten around the two DOIs above, with an
explicit statement of which of the seven rows survive. §11 item 7 must be
deleted: benchmarking against Guha et al. (two half-spaces) and Barak et al.
(two half-spaces) does not test a three-region code, because neither has an
interphase. The correct benchmarks are the two sandwiched-layer papers.

---

### M-1. Internal arithmetic inconsistency: seven interface conditions cannot close a three-region problem

**Location:** §2 item 2 ("four modes and seven interface conditions already
addressed"), §5, §7 ("the count from two-half-space scattering cannot
automatically be carried over to a three-region problem").

**The defect.** The document is aware of the hazard and then reproduces it. The
count is not a matter of opinion; it follows from the polynomial degree.

**Derivation (Audit Part G).** For the four-field model:

```
u_x, u_z  : 2nd order each, coupled  -> 4 roots
phi       : quasi-static Gauss        -> 2 roots
theta     : 2nd order in z            -> 2 roots
                                   n = 8
```

One interface must impose `n = 8` continuity quantities:
`u_x, u_z, σ_xz, σ_zz, φ, D_z, θ, q_z`.

- Two half-spaces: 4 reflected + 4 transmitted = 8 unknowns ⇒ C = 8. ✅
- Three regions: 4 + 8 + 4 = **16 unknowns**, 2 interfaces ⇒ **C = 8 per
  interface, 16 total**. ✅

A seven-condition interface leaves rank deficiency one, which in a solver
appears as a spurious mode in `R` (this exact failure mode was hit repeatedly in
the author's own `check_03` work on the same class of problem: modes with
scalar-fraction 1.0000 and `cond(A_ss) = 6.2×10⁶` masquerading as physical
branches).

**Two readings, and only the source paper settles which.** Either the blueprint
has mischaracterised the prior art, or Guha et al. fix one field by an extra
assumption that must be identified and carried over. Either way the
three-region count is 16 and cannot be inherited.

**Required correction.** Replace §5's qualitative warning with the count above,
and state for **each** model (B0, B1, B2) the exact number of interface
conditions, the rank, and the number of admissible roots per region. This is the
one place where the blueprint could earn a genuine analytical contribution at
almost no cost — see the redesign in §4 below.

---

### M-2. Parameter counting is dismissed rhetorically where a theorem is available, and the theorem is the best result in the project

**Location:** §9.2 ("Do not infer a universal impossibility or possibility from
parameter counting alone: modal symmetries, constraints and passive admissibility
are decisive").

**The defect.** The statement is correct in spirit and is left as a warning. But
for the *passive, dissipative* case there is a clean, provable and useful
obstruction, and the document declines to state it.

**Derivation (Audit Part E2).** With one lossless reactive element `Z_e = iX_e`
in series with a layered input impedance `Z_in = R_ + iX_`:

```
|R|² = ((R_ − Z₁)² + (X_ + X_e)²) / ((R_ + Z₁)² + (X_ + X_e)²)
```

minimised over `X_e` at `X_e = −X_`, giving

```
|R|_floor = |Re Z_in − Z₁| / (Re Z_in + Z₁)
```

The reactive element removes the **reactive** mismatch and can never remove the
**resistive** one. Measured:

| case | min(Re Z_in − Z₁) | floor with ideal reactance | min\|R\| without element | needs gain? |
|---|---|---|---|---|
| lossless | +2.6×10⁻³ Pa | 4.4×10⁻¹¹ | 3.5×10⁻⁶ | no (already zero) |
| Q = 10⁴ | +1.05×10³ Pa | 1.77×10⁻⁵ | 1.80×10⁻⁵ | **yes** |
| Q = 10³ | +1.05×10⁴ Pa | 1.77×10⁻⁴ | 1.77×10⁻⁴ | **yes** |
| Q = 10² | +1.02×10⁵ Pa | 1.72×10⁻³ | 1.75×10⁻³ | **yes** |

Two consequences:

1. **In the lossless limit the element buys nothing** (floor already 4×10⁻¹¹).
2. **With loss the floor is strictly positive and only negative resistance —
   i.e. gain, which §1 explicitly excludes — can close it.**

This kills the intuitive hope that B1's "generalized imperfect contact" is a
matching resource. Passive interface elements redistribute mismatch; they cannot
absorb it.

**Integrity note on this audit.** An earlier draft of this finding claimed the
opposite — that one reactive element *restores* exact cancellation with loss. It
was wrong: `Re Z_in − Z₁` is positive for every `h` in the scan, so the floor
never closes. The wrong version has been replaced in the script, and the
correction is recorded here rather than silently dropped.

**Required correction.** Promote this from a warning to a stated proposition in
§9.2, as an obstruction for the declared interface family (option (b) of the
blueprint's own §9.2). It is rigorous, it is short, it is inside the author's
field, and it converts the project's main weakness into its main result.

---

### M-3. The "material-pair gate" is the real feasibility gate, and it fails for the only physically motivated pair

**Location:** §3 ("Material-pair gate (not yet fixed)").

**The defect.** The gate is stated as a data-availability problem ("documented
set of compatible parameters"). It is actually a **realisability** problem, and
it fails before any data question arises.

**Evidence (Audit Part D).** Exact quarter-wave matching needs a real
homogeneous solid whose longitudinal impedance is `√(Z₁Z₂)`:

| pair | Z₁ | Z₂ | required Z_m | nearest tabulated solid | error |
|---|---|---|---|---|---|
| PZT-5H \| water | 29.6 | 1.50 | **6.67** | epoxy 2.95 | 56 % |
| PZT-5H \| steel | 29.6 | 46.3 | 37.04 | brass 37.66 | 1.7 % |
| PZT-5H \| Al | 29.6 | 17.1 | 22.48 | Al 17.06 | 24 % |
| PZT-5H \| epoxy | 29.6 | 2.95 | 9.35 | glass 14.00 | 50 % |

(MRayl; c_L used for the impedance, consistent with normal incidence.)

The pair with the strongest physical motivation — a piezoceramic radiating into
a fluid load, which is the transducer-coupling story §13 invokes — has a required
impedance that **no bulk solid possesses**. Meanwhile the pair that matches
(PZT/steel, 1.7 %) is a solid–solid pair where a 1.7 % impedance error is already
a ~0.8 % power reflection, i.e. worse than the off-normal floor in C-2.

Worse, for a three-region model in which *all three* phases carry the full
piezothermoelastic field set, the interphase must itself be piezothermoelastic,
which introduces its own `e`, `ε`, `β`, `p` constants. Those are not free
parameters of a matching network; they are coupled to `C` and `ρ` by the same
thermodynamic potential. The realisability problem is therefore harder than the
table suggests, not easier.

**Required correction.** Replace the data-availability gate with a
**realisability gate**: state `Z_m^required = √(Z₁Z₂)` together with the c_L and ρ
of the candidate interphase, and require the pair to be *shown* realisable before
any optimisation is attempted. If no real pair exists, say so in §3 and make
that a scope decision, not a later discovery. Do not adopt "a controlled
theoretical parameter family clearly labelled illustrative" as the fallback
without stating that this concedes the realisability question.

---

### M-4. `R_E + T_E + A_E = 1` is asserted as a budget before being defined, and in the intended regime it degenerates

**Location:** §8.

**The defect.** §8 shows
`R_E + T_E + A_E = 1` in a display equation and then attaches the condition
"*only after* a suitable stationary, source-free integral balance is derived".
The ordering is correct in words and wrong typographically: the identity is
displayed first, and the caveat follows in italics. A reader, and a scoring
referee, will read the equation.

**Second problem.** Given C-3, `A_E ~ 10⁻⁹` for Fourier or LS at ultrasonic
frequencies. The budget therefore collapses to `R_E + T_E = 1`, which is the
lossless elastic budget. §8's entire apparatus — "thermal dissipation",
"stored-energy/entropy-flux accounting", "check the sign of dissipative work" —
is then describing a term that is numerically zero.

**Evidence.** Audit Parts E, F, F5: in the lossless elastic reduction the residual
`R_E + T_E − 1` is ≤ 7×10⁻¹⁶ at every one of nine configurations tested. There
is no dissipative term to account for.

**Required correction.**
1. Delete the display equation from §8 and state the budget in words only, with
   the identity appearing *after* the balance is derived — or, better, present
   the **modified** balance with explicit boundary/source terms as the primary
   object, since that is what the document says it will do.
2. Add: "in the Fourier/LS limit at ultrasonic frequencies `A_E < 10⁻⁸`; the
   budget is numerically the lossless one, and the absorbing-comparison arm of
   H1 is therefore vacuous." Then either delete that arm of H1 or introduce an
   explicit loss mechanism (a complex modulus, or an interphase conductance)
   that is *stated as a modelling choice*, not derived from thermoelasticity.
3. §8 must separate mechanical, electroelastic and thermal contributions to the
   flux *with their orders of magnitude*, in the manner of the author's own
   Barak et al. work. Otherwise "coupled energy flux" is a phrase, not a
   quantity.

---

### M-5. The verification-gate list contains a gate that cannot be passed and a gate that is already passed

**Location:** §11 items 1–8.

**Item 1** requires `h → 0` to recover the direct-welded interface "account for
any singular interphase scaling". With a spring/contact layer this is a
**singular perturbation**: the limit `K → ∞` is not attained uniformly, and
`h → 0` with fixed `K` drives `R → R_weld` only at rate `O(Kh/ω)`. As written the
gate is untestable ("appropriate" is doing the work). Specify the asymptotic
order required, e.g. `|R(h) − R_weld| = O(h)` with a stated constant.

**Item 6** requires comparison to "the scalar quarter-wave analytical formula".
Per C-1 this is already exact to 6×10⁻¹⁰. It is a *code check*, not a
*validation*, and listing it among the gates overstates the verification
apparatus.

**Item 7 is invalid** — see C-4. Guha et al. and Barak et al. are two-half-space
problems in which no interphase exists, so they cannot validate a three-region
solver's interphase block. They can validate the *bulk* mode solver, and §11
should say exactly that.

**Required correction.** Rewrite §11 as three lists: (a) code checks (item 6
belongs here); (b) reduction limits with stated asymptotic rates (items 1–5);
(c) external validation, which for the interphase requires the sandwiched-layer
papers and currently cannot be performed. Do not present a validation gate that
depends on a paper that does not test the code path.

---

### M-6. The novelty budget is spent in the wrong place: the physics-free parts are new, the physics is not

**Location:** whole document.

**The defect.** Reading the surviving novelty (C-4) against the inert field
(C-3):

- the **piezoelectric** coupling (as against electroelastic) is new relative to
  the 2023 competitor, and is the physically real part;
- the **impedance-matching optimisation framing** is new, and is the part with a
  textbook answer at normal incidence (C-1) and a 1–2 % floor off it (C-2);
- the **thermal field and the dissipation budget**, which occupy §4, §8 and
  study S4, are the parts the document argues hardest for, and are inert.

A referee will read §4 and §8 as complexity, not content. The document's own
§14 correctly warns against combining theories "merely to increase mathematical
complexity"; the warning is aimed at nonlocal/surface effects while the
thermoelastic machinery in §4 and §8 is exactly that kind of addition.

**Required correction.** Rebalance: keep piezothermoelasticity because it is the
author's identity and because the *piezoelectric* part of it is load-bearing; but
stop presenting the thermal field and A_E as contributors to the matching
problem. Instead make the **relative inertness of the thermal channel a stated
result**, with the numbers from C-3. "Which couplings can and cannot be used as
matching knobs, and why" is a stronger and more defensible contribution than
"we added temperature".

---

### m-1. The file does not compile

`\setlength{\headheight}{15pt}` appears as bare `setlength{...}` (missing
backslash). This is a hard error, not a warning:

```
! Undefined control sequence.
l.8 setlength{\headheight}{15pt}
```

Correct to `\setlength{\headheight}{15pt}`. Related: `\jump` and `\todo` are
defined in the preamble and never used — harmless, but remove them or use them.

I could not run a full compile (no TeX distribution in this environment), so this
is a source inspection, not a compiler transcript. The error is unambiguous.

---

### m-2. "Exactly four scalar primary unknowns" is right for the wrong reason, and the phrase will be attacked

**Location:** §4.

The count is correct (n = 8). But the sentence "corresponding to three physical
categories (mechanical, electrical, thermal)" invites the standard referee
objection that mechanical displacement components were counted as one category
while this is presented as a four-field theory. State the count as
**two displacement components + one electric potential + one temperature**, and
say explicitly that plane strain reduces the mechanical field from three
components to two. That is precise and pre-empts the objection.

---

### m-3. §14 defers exactly the work the author already owns

§14 defers nonlocal/strain-gradient theory to "a possible separate second paper".
But `10.1108/hff-04-2022-0259` (Gupta, Kumar, Kumar, Pathania, Barak) is the
author's own nonlocal piezothermoelastic R/T work. Deferring it reads as
deferring the author's own differentiator, and a referee who finds that paper
will ask why the first paper is the weaker one. Either bring the verified
nonlocal term into the *third region only* (a specific, contained, justified use)
or state plainly that the first paper deliberately uses the local limit so that
the interphase result is not confounded with a length scale. The second is
defensible; silence is not.

---

### m-4. Section count and figure budget still describe a programme

§12 targets "approximately 7 main figures" across seven sections covering
formulation, roots, scattering, analytical matching, numerics, sensitivity. For a
journal paper that is a stretch but possible; combined with §10's seven studies
(S1–S7) it is two papers. Note that S7 is explicitly "only after base
validation", so S1–S6 is the real budget — six studies, seven figures, three
tables, one analytical result that turns out to be a textbook identity. Reduce to
S1, S2, S3, S5 plus the obstruction; make S4 a consistency table and S6 a single
number (see m-5).

---

### m-5. The surface-effect gate can be closed with one line — close it

**Location:** §6.3, §10 S6.

The gate is stated qualitatively with a formula `η_s ~ k_∥ μ_s / μ_bulk` and an
instruction to check whether it is "not negligibly small". It is decidable:

```
η_s = 2π ℓ_s/λ,   ℓ_s = μ_s/μ_bulk
```

| μ_s | ℓ_s | η_s at λ = 1 mm | η_s at λ = 1 μm | λ for η_s = 1 |
|---|---|---|---|---|
| 1 N/m | 3.9×10⁻¹¹ m | 2.4×10⁻⁷ | 2.4×10⁻⁴ | 2.4×10⁻¹⁰ m |
| 5 N/m | 1.9×10⁻¹⁰ m | 1.2×10⁻⁶ | 1.2×10⁻³ | 1.2×10⁻⁹ m |
| 10 N/m | 3.8×10⁻¹⁰ m | 2.4×10⁻⁶ | 2.4×10⁻³ | 2.4×10⁻⁹ m |

O(1) requires λ ~ 1 nm. The gate is closed at mm wavelengths by three orders of
magnitude. Delete S6 as a study, replace with one sentence and this table, and
move the saved space to M-2's obstruction.

---

### m-6. Minor precision items

- §3 writes the state vector as `q^(j) = (u_x, u_z, φ, θ)ᵀ` and then §4.1 uses
  `q_h` for heat flux with the note "distinct from the primary state vector q".
  Two objects named `q` in a document that prides itself on precision; use `𝐯`
  for the state vector or `𝐪` for heat flux, not both as `q`.
- §5 assumes "conserved tangential `k_x`" without stating that this requires a
  time-invariant, translation-invariant interface — true here, but it is the
  assumption that fails first if any temporal modulation is ever added, and §1
  lists temporal modulation as out of scope, so say why.
- §9.1's `k_mh = (2n+1)π/2` mixes `k_m` (a wavenumber in the interphase) with the
  `h` used for the layer thickness elsewhere; consistent, but the document uses
  `h` both for thickness and (implicitly) for the heat-transfer coefficient
  `H_t`. Already disambiguated by capitalisation; keep it that way.
- §2 says the three sources' "title/DOI association cross-checked". All three
  check out (table in §1 above). State the check date and the registry used, as
  the author's own standards require.

---

## 3. What is right and should be kept

Stated so the redesign does not discard it:

- The **anti-reflection / negative-reflection separation** in §1 is correct and
  important, and was the right response to the earlier critique.
- **H0 as an acceptable outcome** is the single best sentence in the document,
  and M-2 is the proof that makes it pay.
- The **retraction of Gurtin–Murdoch as the principal knob** is correct and the
  order-of-magnitude argument in §6.3 is sound (m-5).
- **B0 before B1 before B2** is the right ordering.
- The **insistence on counts, rank and conditioning** in §5/§7 is exactly the
  right instinct; M-1 and M-2 are that instinct applied properly.
- **§13's refusal to claim device performance** is correct.
- The three cited works are **real and correctly described**. Cite them in the
  body.

---

## 4. Required correction, in order

**Step 1 (blocking).** Rewrite §2 around `10.1016/j.compstruct.2023.117564` and
`10.1016/j.tws.2021.108306`. Delete §11 item 7. Produce the seven-row overlap
table of C-4 and delete every novelty claim that does not survive it.

**Step 2 (blocking).** Replace §9.1's role. Add M-2 as a proposition in §9.2.
The paper's analytical claim becomes: *for the declared passive, dissipative
three-region family, the reflection floor is `|Re Z_in − Z₁|/(Re Z_in + Z₁)` and is
strictly positive; exact anti-reflection is attainable only in the lossless
normal-incidence limit, where it is the classical quarter-wave identity.* That is
a theorem, it is inside applied mathematics, it is inside the author's field, and
it is publishable.

**Step 3 (blocking).** Rewrite §1's H1 to name a frequency band, an angular
window, and a numeric `R_E` target, informed by C-2's table. Delete the
"absorption-dominated solution" comparison or add a stated loss mechanism (M-4).

**Step 4 (blocking).** Insert the DOF subsection of M-1 into §5, with per-model
counts for B0/B1/B2. This is the cheapest genuine contribution available, because
the author's own `check_03` work already establishes the machinery and the
failure modes.

**Step 5 (major).** Insert the screening table of C-3 into §4, and demote the
thermal field accordingly. Delete S4 as a results study.

**Step 6 (major).** Replace the material-pair gate with the realisability gate of
M-3. State the pair or state the concession.

**Step 7 (minor).** Fix `\setlength`. Fix `q`/`𝐯`. Insert full bibliographic
data for the three verified works. Rewrite §11 as three lists. Close S6 with the
m-5 table.

**Step 8.** After steps 1–7, the document is a materially different and stronger
proposal. Its honest scope is: *a rigorous passive-matching obstruction for
coupled-field layers, with the quarter-wave identity as the degenerate limit, the
material realisability bound for a piezoceramic interphase, and a verification
that the thermal channel is inert as a matching resource.* That is one paper with
one theorem, two bounds and one negative result. It is a good paper. It is not
the paper currently described.

---

## 5. Decision

**REDESIGN.**

Not a conditional revision, because three of the four critical findings
(C-1, C-3, C-4) attack the *identity* of the proposal rather than its execution,
and no amount of rewriting §11 or adding figures repairs an identity.

Not a flat reject, because M-2 is a real result waiting to be stated, the
author's method base and publication record support it, and the null-hypothesis
framing already anticipates it.

**Open questions the author must answer before the next version** (none of which
the audit can settle from outside the full texts):

1. Does Guha et al. (2019) use eight interface conditions or seven, and if seven,
   which field is fixed by an extra assumption? This determines whether M-1 is a
   correction to the blueprint's description or a correction to the prior art.
2. Is a piezothermoelastic interphase with `Z ≈ √(Z₁Z₂)` realisable at all, or is
   the interphase necessarily a third distinct material whose own constitutive
   constants move it off the matching value? M-3 says this is the binding
   constraint; it needs a table, not an assumption.
3. Is the intended angular window where C-2's floor is 10⁻³ or 10⁻²? This decides
   whether the paper claims "anti-reflection" or "reflection reduction".
4. Does the author's `10.1108/hff-04-2022-0259` already contain the mode-count and
   energy-flux machinery needed for M-1 and M-2? If so, the cost of the redesign
   is a fraction of the cost of the original plan.

---

*Audit performed 4 November 2026. Every quantitative claim is reproducible via
`python3 code/check_04_blueprint_audit.py`. Self-checks enforced: uniform-medium
limit (|R| = 3×10⁻¹⁶), lossless energy balance (residual ≤ 7×10⁻¹⁶), matched-layer
quarter-wave limit (R_E = 6.3×10⁻¹⁰). One error made in the course of this audit
was found and corrected by its own dimension check; it is recorded in M-2 rather
than removed.*
