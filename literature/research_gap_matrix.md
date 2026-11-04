# Research Gap Matrix — Blueprint V5.0, Phase 1

**Status labels.** `VERIFIED` = read in a source record (Crossref metadata and/or
abstract). `DERIVED` = derived in this Phase. `ASSUMED` = stated assumption, not
supported by a source. `UNRESOLVED` = could not be settled in Phase 1.

**Rule applied throughout.** No novelty claim is made from the *absence* of a
keyword in a title or abstract. Every row below records what was actually read,
and every "gap" is phrased as a question to be tested against full text
(gate G0), not as an established absence.

**Read status caveat, stated once and applying to every row.** For all ten
records in `verified_references.csv`, the *Crossref metadata* was verified
(title, authors, year, journal, volume, pages, DOI). Full text was **not**
accessed for any record. R03 and R05 additionally had their abstracts read.
Consequently:

- Every "relationship" judgement here is made from metadata and abstract, not
  from equations or figures.
- Gate G0 (full-text prior-art comparison) is **NOT passed**. It is
  `UNRESOLVED` and is carried forward to Phase 2.
- Nothing in this matrix may be used as a novelty claim.

---

## 1. Direct competitors

Direct competitors are studies that combine (i) two dissimilar *piezothermoelastic*
half-spaces, (ii) an *imperfect or generalized contact*, and (iii) reflection and
transmission coefficients with energy accounting.

| ID | Study | Coupled fields | Thermal theory | Contact | R/T + energy | Why it is a direct competitor | Evidence |
|----|-------|----------------|----------------|---------|--------------|-------------------------------|----------|
| R02 | Guha, Singh & Das (2019/2021) | u, φ, θ (piezothermoelastic) | LS **and** GL | mechanically, electrically **and** thermally imperfect | yes | Same three physical categories; same contact class | `VERIFIED` metadata |
| R03 | Karmakar, Sahu & Goyal (2023) | u, φ, θ, FG gradation | LS **and** GL | mechanical + electrical + thermal loose bonding | yes, with conservation checked graphically | Closest published analogue of the present bulk problem; explicitly retains four branches | `VERIFIED` metadata + abstract |
| R09 | Guo, Wei, Li & Tang (2015) | u, φ | isothermal | mechanical **and dielectric** imperfection | yes | Origin of the dielectric-imperfection contact law; no thermal field | `VERIFIED` metadata |

**Assessment for R02/R03/R09.** All three already compute exactly the quantity the
present work computes (R and T for a coupled half-space pair with an imperfect
contact). The present work therefore **cannot** claim novelty for:

- the four-field coupled formulation,
- the spring/dashpot contact law,
- the dielectric-imperfection (capacitive-gap) law,
- computing R and T vs angle for such a pair,
- energy ratios in a piezothermoelastic half-space pair.

These are recorded as **VERIFIED LITERATURE OVERLAP** and are treated as
established prior art, not as openings.

## 2. Contextual studies (adjacent, not competitors)

| ID | Study | Why contextual | Evidence |
|----|-------|----------------|----------|
| R04 | Barak, Kumar, Kumar & Gupta (2022) | Heterogeneous mode count (2 reflected, 4 transmitted) because medium 1 is purely elastic; nonlocal + memory + two-temperature, all outside Phase-1 scope | `VERIFIED` metadata + abstract |
| R05 | Gupta, Kumar, Kumar & Barak (2023) | Welded interface; higher-order two-temperature + dual-phase-lag + memory-dependent derivatives; nonlocal | `VERIFIED` metadata + abstract |
| R06 | Ghazaryan, Piliposyan & Piliposian (2025) | SH waves in a *periodically stratified* piezoelectric reflector; the mechanism is periodicity (Bragg), not a contact law | `VERIFIED` metadata |
| R07 | Wang et al. (2023) | Nanoplate sandwich, nonlocal integral elasticity; different structural class | `VERIFIED` metadata |
| R08 | Tuan, Dung & Nga (2025) | Weakly nonlocal piezoelectric half-spaces + imperfect interface **plus an inverse problem** | `VERIFIED` metadata |
| R10 | Goyal, Sahu & Mondal (2020) | FG piezoelectric half-spaces, QP wave only, isothermal | `VERIFIED` metadata |
| R11 | Pang & Liu (2011) | Piezoelectric/piezomagnetic pair, isothermal | `VERIFIED` metadata |

## 3. What each competitor does *not* appear to do

Each statement is a **question for gate G0**, not a finding.

| # | Question (candidate gap) | Nearest competitor | Why it might still be open | Evidence status |
|---|--------------------------|--------------------|-----------------------------|-----------------|
| G-1 | Is there an explicit, falsifiable *condition* for reflection **suppression** (not merely a computed R curve that happens to dip) with a declared transmission cost? | R02, R03 | Both compute R vs parameters and show dips. Neither abstract states a suppression *criterion* of the form "ΔR ≤ −ε_R with ΔT ≥ −δ_T". | `UNRESOLVED` — needs full text and figure inspection |
| G-2 | Is the reflected-power reduction ever separated from *dissipation*, i.e. is A reported alongside R and T for the contact itself? | R03, R07 | R03 checks conservation graphically; R07 reports a dissipation fraction but for a nonlocal nanoplate, not a contact. Whether the contact dissipation is separated in R02/R03 is unknown. | `UNRESOLVED` |
| G-3 | Is the scalar no-go statement (a passive zero-thickness contact cannot increase transmitted power over the welded baseline) stated anywhere? | none found | The statement is elementary transmission-line algebra and may well be folklore in the impedance-matching literature, which was **not** searched in Phase 1. | `UNRESOLVED` — an impedance-matching literature search is required before any novelty wording |
| G-4 | Does any competitor treat the finite interphase as the *primary* configuration with the two-half-space imperfect interface as a declared baseline? | R02, R03, R06 | R06 is layered but periodic and SH; R02/R03 are two-half-space. A finite piezothermoelastic *matching* interphase between two half-spaces does not appear in the ten records. | `UNRESOLVED` — absence in ten records is not absence in the literature |
| G-5 | Is the h → 0 reduction of the interphase-to-contact composition derived anywhere? | none found | Standard in asymptotic homogenisation of thin interphases; that literature was not searched. | `UNRESOLVED` |
| G-6 | Does any competitor report the *independent-equation count* and rank of the interface system rather than asserting a mode count? | none found | Likely treated as routine in the analytical mechanics literature. | `UNRESOLVED` |

## 4. Gaps that Phase 1 did **close** (DERIVED, not literature-dependent)

These are mathematical results obtained here; they do not depend on the
literature search and are not novelty claims about the literature.

| # | Result | Where | Status |
|---|--------|-------|--------|
| D-1 | The 6 mm piezothermoelastic characteristic determinant is a polynomial of degree 8 in ξ, **even** in ξ, of the form P₄(ξ²), with exactly five non-zero coefficients (ξ⁰, ξ², ξ⁴, ξ⁶, ξ⁸). | `theory/characteristic_matrix.tex`, check C1 | `DERIVED` |
| D-2 | Root count: 8 roots = 4 branches × 2 directions. The four branches are quasi-P, quasi-SV, an electrically evanescent branch with ξ ≈ i k √(ε₁₁/ε₃₃), and a thermal branch with ξ ≈ (1+i)√(ωρc_E/2κ). | check C2 | `DERIVED` + `TESTED` |
| D-3 | At ξ = 0 (propagation normal to the poling axis) the system block-decomposes into {u_x, θ} ⊕ {u_z, φ}: the shear-thermal block has **no** piezoelectric coupling and the longitudinal block has **no** thermal term. At k = 0 it decomposes into {u_x} ⊕ {u_z, φ, θ}. | check C1 | `DERIVED` |
| D-4 | The only gauge root is the constant-potential / zero-D_z state at (k, ξ) = (0, 0). It does not persist for k ≠ 0, so normal incidence must be taken as a limit. | check C1 | `DERIVED` |
| D-5 | The Lord–Shulman equation of Blueprint Eq. (5), **including** the placement of (1 + τ_q ∂_t) around the whole bracket, follows exactly from a single entropy balance plus one Cattaneo flux law. | `theory/governing_equations.tex`, check B2 | `DERIVED` + `TESTED` |
| D-6 | The constitutive set (σ, D, ρη) is generated by a single Helmholtz potential, with all Maxwell cross-relations closing. | `theory/governing_equations.tex`, check B1 | `DERIVED` + `TESTED` |
| D-7 | A causal floor τ_q ≥ κ/(ρ c_E c_L²) exists below which the LS thermal wave is superluminal. All four material sets are checked against it. | check H1 | `DERIVED` + `TESTED` |

## 5. Honest summary

- **Prior art is dense.** Three direct competitors already do coupled
  piezothermoelastic R/T with imperfect contacts. The formulation is not new.
- **The only defensible candidate contribution** is a *criterion* for
  reflection suppression with an explicit transmission-cost budget, together with
  a finite interphase as the primary configuration. Whether that is genuinely
  absent is **`UNRESOLVED`** until gate G0 (full-text) and an impedance-matching
  literature search are done. Phase 1 does **not** claim it is absent.
- **A large secondary contribution is the verification apparatus itself**: the
  exact determinant structure, the root/branch/parity inventory, the corrected LS
  reduction, the flux functional verified by two independent routes, and the
  identification of three numerical traps (see `reports/UNRESOLVED_ISSUES.md`).
  These are useful and honest, but they are *not* a novelty claim.
