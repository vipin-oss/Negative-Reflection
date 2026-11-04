# Validation and Benchmark Plan — Blueprint V5.0, Phase 1 (Task H)

**Purpose.** Blueprint V5.0 §7 makes validation a compulsory *result*: no
statement of agreement is permitted without an aligned quantity, a stated source
and a difference. This document specifies every benchmark to be run, what it
tests, what answer is expected, and what counts as acceptance — **before** any of
them has been run in the coupled case. Statuses below are the Phase-1 status.

**Global rule.** An energy residual is an internal consistency check, not
independent validation. Every entry states which it is.

**Status legend.** `CLOSED` = executed and passing in Phase 1.
`PARTIAL` = partly executed. `PLANNED` = specified, not executed.
`CONDITIONAL` = blocked pending a compatibility check or unavailable data.
`UNAVAILABLE` = cannot be done with accessible resources.

---

## A. Analytical limit benchmarks (route V2 — independent of any publication)

These require no external source: the reference value is a closed-form result
from classical elastodynamics or transmission-line theory. They are the strongest
class of check available and must all be closed before any coupled result is
quoted.

### B-01 — Identical welded half-spaces
| Field | Value |
|---|---|
| Assumptions | Both half-spaces the same material; welded contact; any incidence angle |
| Expected | The reflected field vanishes identically: every reflected amplitude is zero. With the flux definition in `theory/energy_flux_derivation.tex`, R → 0 at round-off |
| Source | Self-evident from translation invariance; no external source needed |
| Method | Assemble the interface operator for identical media, solve, evaluate R |
| Acceptance | R < 1e-14 (round-off only). Any larger value indicates a gauge or mode-selection error, not a physical effect |
| Type | **Internal consistency** (not independent validation) |
| Phase-1 status | `PARTIAL` — the exact roots and null vectors are verified for identical media (checks C2, E1), but the scattering solve is not built. The G1 scalar benchmark reproduces R = 0 exactly for Z₁ = Z₂ |

### B-02 — Uncoupled elastic normal incidence (Fresnel)
| Field | Value |
|---|---|
| Assumptions | e = β = p = 0 (as a *limit*, see warning below); normal incidence; lossless |
| Expected | r = (Z₂−Z₁)/(Z₂+Z₁), R = \|r\|², T = 1−R, with Z = ρc_L |
| Source | Classical elastodynamics; also reproduced analytically by Blueprint Eqs. (9)–(10) with the contact switched off |
| Method | Full mode-matching solve in the elastic limit, over an impedance-ratio sweep |
| Acceptance | \|ΔR\| < 1e-9 and \|R+T−1\| < 1e-12 |
| Type | **Independent validation** (external formula) |
| Phase-1 status | **`CLOSED` for the determinant and root structure** (check I1: exact roots match closed-form elastic branches to 1.7e-15 over 24 cases). **`UNAVAILABLE` for the full scattering solve** — see U-3 |

> **Critical procedural warning (found in Phase 1).** Setting e = β = p = 0
> *exactly* makes the 8×8 interface system **singular**: φ and θ then appear in
> no other equation, so four of the eight conditions become identically zero rows.
> The decoupled case must be reached as a limit with couplings scaled by a small
> λ → 0. Any benchmark that sets the couplings to zero will appear to fail for a
> reason that has nothing to do with the physics.

### B-03 — Classical elastic P–SV oblique incidence, welded
| Field | Value |
|---|---|
| Assumptions | Isotropic elastic half-spaces; welded; oblique incidence; critical angle and evanescent treatment exercised |
| Expected | The classical P–SV reflection/transmission coefficients (Aki & Richards / Achenbach closed forms) |
| Source | Standard seismology/elastic-wave texts; must be re-derived independently and the derivation recorded |
| Method | Assemble the 4×4 mechanical system, solve, compare both amplitudes r_PP, r_PS and the energy sum |
| Acceptance | \|Δr\| < 1e-10 per coefficient; \|R+T−1\| < 1e-12 below the critical angle |
| Type | **Independent validation** |
| Phase-1 status | `PLANNED`. Note the PZT-5H material is **not** sagittally isotropic (C₁₃ ≠ C₁₁−2C₄₄), so an isotropic closed form may not be used for it — this cost one failed check in Phase 1 before being spotted |

### B-04 — Scalar reactive / resistive contact propositions (Proposition 1)
| Field | Value |
|---|---|
| Assumptions | S1–S4 of Blueprint §5.2: one channel each side, normal incidence, real positive impedances, zero-thickness local passive contact, two topologies |
| Expected | (a) R ≥ R₀ for reactive contacts, equality iff B = 0; (b) exact R = 0 requires dissipation of the right topology, giving T = Z_min/Z_max and A = 1 − T; (c) T ≤ T_w for every passive contact, equality iff the contact vanishes |
| Source | Derived in Blueprint §5.3–5.4; **transmission-line algebra, not a new theorem** |
| Method | Analytic proof plus direct numerical evaluation; random sampling over passive contacts |
| Acceptance | Analytic identities to machine precision; random sweep must find no violation |
| Type | **Internal consistency** of a derived lemma |
| Phase-1 status | **`CLOSED`** (check G1). Measured: reactive sweep minimum R = R₀ exactly at \|B\| = 0; resistive matching gives R = 0, T = 0.5 = Z_min/Z_max, A = 0.5, R+T+A−1 = 0; Voigt mapping with K = 0, C = Z₁Z₂/(Z₂−Z₁) gives R = 0; 200 000 random passive contacts give max(T − T_w) = −2.5e-3 ≤ 0 |

### B-05 — Quarter-wave layer
| Field | Value |
|---|---|
| Assumptions | Lossless homogeneous layer of impedance Z_ℓ and phase φ = kd between the same two media |
| Expected | Z_in = Z_ℓ(Z₂ + iZ_ℓ tanφ)/(Z_ℓ + iZ₂ tanφ); R = 0 iff Z_ℓ² = Z₁Z₂ at φ = π/2, with T = 1 > T_w |
| Source | Standard transmission-line result; **must not be called novel** (Blueprint §5.4 says so explicitly) |
| Method | Direct evaluation of the input-impedance formula |
| Acceptance | R < 1e-18 at φ = π/2 − 1e-9 |
| Type | **Internal consistency** of a standard result |
| Phase-1 status | **`CLOSED`** (check G1): R = 1.25e-19, T = 1.000000000000 |

### B-06 — Lord–Shulman to Fourier reduction
| Field | Value |
|---|---|
| Assumptions | τ_q → 0 |
| Expected | The LS equation reduces to the classical Fourier conduction equation coupled to the same sources |
| Source | Definitional |
| Method | Symbolic substitution τ_q = 0 in the mode-form LS equation |
| Acceptance | The symbol Λ → 1 exactly and the thermal branch wavenumber becomes the parabolic-diffusion form |
| Type | **Internal consistency** |
| Phase-1 status | **`CLOSED` analytically** (check B2 shows Λ = 1 − iωτ_q; at τ_q = 0, Λ = 1). Note this is *not* the same as deleting the thermal field |

### B-07 — Coupling reduction β, p, e → 0
| Field | Value |
|---|---|
| Assumptions | Couplings scaled by λ → 0 with λ > 0 |
| Expected | The coupled system converges to the uncoupled one; the deviation is O(λ²) because the piezo contribution to the modulus is e²/ε |
| Source | Dimensional/asymptotic argument given in the theory notes |
| Method | Solve at λ = 1e-1, 1e-2, 1e-3 and extrapolate |
| Acceptance | Deviation ratio between successive λ decades ≈ 100 (i.e. O(λ²)) |
| Type | **Internal consistency** |
| Phase-1 status | **`UNAVAILABLE`** — see U-3. The λ-scaling was attempted but the amplitude normalisation of the mode-matching solve could not be established, so no reliable λ-trend was obtained |

### B-08 — Thermal / piezo uncoupling of anti-plane SH
| Field | Value |
|---|---|
| Assumptions | 6 mm symmetry, poling along z, anti-plane motion u_y(x,z) |
| Expected | The thermal coupling vanishes **identically** (β₁₃ = β₂₃ = 0 for 6 mm), so SH motion cannot exercise the thermal channel |
| Source | Derived in `theory/characteristic_matrix.tex`, Remark on SH |
| Method | Symbolic substitution of the SH ansatz into the LS equation |
| Acceptance | The thermal source term is exactly zero, symbolically |
| Phase-1 status | **`CLOSED`** — this is a structural fact, verified symbolically. A previous session had treated it as a solver failure; it is not |

### B-09 — Boundary-condition rank and residual
| Field | Value |
|---|---|
| Assumptions | Any configuration |
| Expected | The number of independent interface relations equals the number of unknown amplitudes after gauge removal; the assembled operator has full rank; the normalised residual ‖Ba − b‖ is small |
| Source | Derived per configuration in `theory/interface_conditions.tex` |
| Method | Rank by exact elimination on the equilibrated operator; residual of the solve |
| Acceptance | Rank = 8 (baseline) / 16 (interphase); residual < 1e-10 after equilibration; singular values reported |
| Type | **Internal consistency** |
| Phase-1 status | `PARTIAL` — **rank 8 confirmed** for the welded piezothermoelastic interface at three angles (check I2). The *residual* is `PLANNED` |

### B-10 — h → 0 interphase reduction
| Field | Value |
|---|---|
| Assumptions | Homogeneous interphase 0 < z < h with passive faces, h → 0 |
| Expected | Composite map s(h⁺) = P₂ e^{hA} P₁ s(0⁻) reduces to a jump law whose effective compliance is the sum of the two face compliances plus an O(h) layer contribution |
| Source | Derived to first order in `theory/interface_conditions.tex`, Proposition (h→0) |
| Method | Numerical composition at h = 1e-3, 1e-4, 1e-5 × λ and extrapolation to h = 0; compare against the direct two-half-space imperfect-contact solve |
| Acceptance | The h → 0 limit coincides with the **imperfect-contact** baseline (not the welded one) to within the O(h) term; the measured slope in h matches the predicted layer compliance |
| Type | **Internal consistency** plus an **independent consistency** check against the baseline solve |
| Phase-1 status | `PLANNED`. First-order reduction `CLOSED` analytically; the second-order term is `UNRESOLVED` (U-6) |

### B-11 — Perfect-versus-imperfect contact along stated paths
| Field | Value |
|---|---|
| Assumptions | K, C, H_T → ∞ (perfect) versus → 0 (free) |
| Expected | As K, C, H_T → ∞ the imperfect solution converges to the welded solution; the convergence path distinguishes stiffness (K) from damping (C) and from thermal conductance (H_T) |
| Source | Definitional limit of the contact laws |
| Method | Sweep each contact parameter over 6 decades and record the deviation from the welded solve |
| Acceptance | Monotone convergence with the correct order; the three parameters must not be conflated |
| Type | **Internal consistency** |
| Phase-1 status | `PLANNED` |

---

## B. Published-comparison benchmarks (route V1 — conditional)

Blueprint §7.1 is explicit: **no single published paper is mandated**, and a
published result may be used only after a documented compatibility check
establishing, from the *full accessible text*, the bulk constitutive model, the
thermal theory, the modes retained on each side, the boundary-condition type and
sign conventions, the material symmetry and poling direction, the incident wave
and angle definition, and which parameters can be switched off to reach the
present model.

**All three candidates below have `read_status = metadata only` in
`literature/verified_references.csv`. A compatibility check is therefore
impossible in Phase 1 and all three are `CONDITIONAL`.**

### B-12 — Guha, Singh & Das, `10.1080/17455030.2019.1610198`
| Field | Value |
|---|---|
| Assumptions / compatibility | Requires full text. Must confirm: LS *and* GL both used (so a switch-off path to LS exists), the exact loose-bonding contact law and its symbol convention, and the energy-ratio normalisation |
| Expected | If compatible, the LS-only limit of their model must reproduce the present R/T at matched parameters |
| Source | *Waves in Random and Complex Media* 31(4):660–689 |
| Method | Digitise their R/T curves; declare an extraction uncertainty; compare at matched (ω, θ) |
| Acceptance | Discrepancy within their reported precision and the extraction uncertainty; ≥3 convergence settings |
| Type | **Independent validation** |
| Phase-1 status | `CONDITIONAL` — full text not accessible; compatibility check not possible |

### B-13 — Karmakar, Sahu & Goyal, `10.1177/10812865221099502`
| Field | Value |
|---|---|
| Assumptions / compatibility | Abstract confirms 4 retained waves (quasi-P, quasi-SV, electroacoustic, thermal), LS **and** GL, all three loose-bonding channels, CdSe and PZT-5A. Still required from full text: the loose-bonding law's exact form and sign, the FG gradation profile (must be switched off to reach the homogeneous case), and the energy-ratio definition |
| Expected | The homogeneous-limit LS-only reduction must reproduce the present R/T |
| Source | *Mathematics and Mechanics of Solids* 28(4):891–919 |
| Method | As B-12; additionally verify that the gradation-off limit is well defined in their formulation |
| Acceptance | As B-12 |
| Type | **Independent validation** |
| Phase-1 status | `CONDITIONAL` — abstract read only; full text not accessible |

### B-14 — Barak et al., `10.1007/s12648-022-02568-w`
| Field | Value |
|---|---|
| Assumptions / compatibility | Abstract states two reflected and four transmitted waves, with medium 1 **elastic**. This is a *heterogeneous* mode count, not the symmetric 4/4 case. To use it, nonlocality, memory and the two-temperature option must all be switched off, the thermal law must be shown to reduce to LS, and the heterogeneous mode count must be derived independently |
| Expected | Only after all of the above: agreement in the stated limit |
| Source | *Indian Journal of Physics* |
| Method | As B-12, plus an independent derivation of the 2/4 mode count |
| Acceptance | As B-12, plus the mode count must be derived and matched, not assumed |
| Type | **Independent validation** |
| Phase-1 status | `CONDITIONAL`. Blueprint §7.1 warns that if the limit cannot be demonstrated the candidate must be marked **NOT USABLE AS BENCHMARK**. That is the current status |

### B-15 — Self-benchmark against the author group's own papers
| Field | Value |
|---|---|
| Assumptions | — |
| Expected | — |
| Source | — |
| Method | — |
| Acceptance | Blueprint §7.1 prefers a *compatible independent* benchmark over an author-group self-benchmark |
| Type | — |
| Phase-1 status | `DEPRIORITISED` by Blueprint §7.1. Recorded here so the choice is explicit rather than accidental |

---

## C. Numerical-correctness benchmarks (route V3)

### B-16 — Eigenmode residual against an independent root formulation
| Field | Value |
|---|---|
| Expected | ‖Ma‖/(‖M‖‖a‖) below round-off-scaled tolerance, where the roots come from an *independent* formulation (closed-form branches or a state-space companion matrix) rather than from the same polynomial |
| Method | Compare polynomial roots against closed-form elastic branches; compute residuals |
| Acceptance | Residual < 1e-20 at 40-digit precision |
| Phase-1 status | **`CLOSED`** — check C2 measures a maximum residual of 2.2e-27 over all 8 roots; check I1 matches roots to closed-form elastic branches over 24 cases to 1.7e-15 |

### B-17 — Convergence in tolerance and resolution
| Field | Value |
|---|---|
| Expected | R, T, A and the location of any minimum are stable under at least three tightened tolerances and two numerical resolutions |
| Acceptance | Changes below the declared reporting tolerance; the *location* of a minimum must be reported with its convergence history |
| Phase-1 status | `PLANNED` (no scattering solve yet) |

### B-18 — Robustness against artificial minima
| Field | Value |
|---|---|
| Expected | Any reflection minimum survives refinement and perturbation; branch swaps, non-physical parameters, wrong attenuation signs and absorption mislabelled as transmission are excluded |
| Acceptance | Minimum persists; A reported alongside ΔR and ΔT |
| Phase-1 status | `PLANNED` |

### B-19 — Energy gates E1/E2/E3
| Gate | What is required | Phase-1 status |
|---|---|---|
| E1 | LS energy identity and normal total-field flux derived, including the thermal term and the τ_q-dependent stored quantities; bulk dissipation positivity verified | **`NOT PASSED`** — Phase 1 stated an invalid linear thermal term (`+½Re[q_z]`) and is withdrawn. The corrected bilinear form `+½Re[q_zΘ*]/T0` is now derived and verified: elastic limit to 1.95e-16, antisymmetry exact, z-conservation to 2e-41, Second-Law sign reproduced to 1.2e-41. The full coupled local-balance identity closes to 7.1e-17 in an independent dps=60 evaluation but only to 1.9e-5 through the repository null-vector routine (open item U-3). Check D2 is reported **UNRESOLVED**. Until that closes, E1 is not passed. |
| E2 | Mechanical, thermal and electrical contact dissipations derived as non-negative forms with units checked and signs fixed | **`NOT PASSED`** — derived and unit-checked (checks F1, G1), but the ½-bookkeeping between the mechanical and thermal forms is not yet unified for the Phase-2 flux normalisation (U-5). Phase 1 wrongly described this as a difference of physical origin; it is a convention question only. |
| E3 | Numerical energy balance evaluated on planes in the half-spaces with convergent residuals | **`UNAVAILABLE`** — requires the coupled solver. **No energy ratio may be quoted until this closes** |

---

## D. Acceptance summary and honest status

| Benchmark | Type | Phase-1 status |
|---|---|---|
| B-01 identical media | internal | PARTIAL |
| B-02 Fresnel elastic | **independent** | **CLOSED in the reduced system** (check J2: 1.67e-16 at normal incidence; energy closure ≤1e-15 over 7 angles). UNAVAILABLE in the coupled system (U-3) |
| B-03 elastic P–SV oblique | **independent** | PLANNED |
| B-04 Proposition 1 | internal | **CLOSED** |
| B-05 quarter-wave | internal | **CLOSED** |
| B-06 LS → Fourier | internal | **CLOSED** |
| B-07 coupling → 0 | internal | **SUPERSEDED** — Phase 1B shows the e=β=p=0 limit is a *reduction* (two gauge variables and their 0=0 conditions), not a λ→0 continuation. Established symbolically by check J1. |
| B-08 SH thermal decoupling | structural | **CLOSED** |
| B-09 BC rank / residual | internal | PARTIAL (rank CLOSED) |
| B-10 h → 0 | internal + independent | PLANNED |
| B-11 perfect/imperfect paths | internal | PLANNED |
| B-12 Guha et al. | **independent** | CONDITIONAL |
| B-13 Karmakar et al. | **independent** | CONDITIONAL |
| B-14 Barak et al. | **independent** | CONDITIONAL / NOT USABLE |
| B-15 self-benchmark | — | DEPRIORITISED |
| B-16 eigen residual | internal | **CLOSED** |
| B-17 convergence | internal | PLANNED |
| B-18 artificial minima | internal | PLANNED |
| E1 flux | gate | **CLOSED** |
| E2 interface dissipation | gate | CLOSED except U-5 |
| E3 numerical balance | gate | **UNAVAILABLE** |

**What this means.** The *analytical* foundation is validated: the determinant
structure, the root inventory, the flux functional and the scalar propositions are
all closed and independently cross-checked. The *scattering* layer is not: no
coupled reflection or transmission coefficient has been computed, and gate G5/E3
remains unavailable. Under Blueprint §7.4, **failure of any energy gate blocks the
reporting of R, T or A** — so Phase 1 reports none.

**Stop rule.** No submission-ready numerical result is generated in Phase 1.
