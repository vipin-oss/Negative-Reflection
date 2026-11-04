# Unresolved Issues — Phase 1, Blueprint V5.0

**Date:** 4 November 2026
**Scope:** items that Phase 1 did **not** close, plus the errors made and corrected
during Phase 1. Nothing here is hidden or softened: an item that failed is listed
as failed.

**Label key.** `ERROR-CORRECTED` = a mistake was made in Phase 1 and fixed;
recorded so it is not repeated. `OPEN` = known and not yet resolved.
`BLOCKING` = prevents a downstream claim. `DEFERRED` = acknowledged, out of Phase-1
scope by instruction.

---

## Part 1 — Errors made and corrected during Phase 1

These are listed first because they are the most useful part of this document. Each
one was a *real* defect in my own work, found by running the checks, not a
hypothetical.

### E-1 — `sympy` `.coeff()` on an unexpanded expression returned an all-zero matrix
**Class:** Numerical uncertainty / mathematical error.
**Symptom.** The first version of the characteristic-matrix prototype reported
`det M ≡ 0` identically, for every material and every parameter set.
**Cause.** `Matrix.coeff()` was called on an expression that had not been
`.expand()`ed. `coeff` only extracts a coefficient when the expression is in
expanded form; on an unexpanded sum it silently returns `0` rather than raising.
The resulting matrix was all zeros in the affected entries and its determinant
vanished.
**Fix.** `sp.expand()` every equation before calling `.coeff()`. Flagged explicitly
in the code.
**Lesson.** A determinant that is *identically* zero is almost always a
construction bug, not physics. Any solver should assert that `det M` is not
identically zero before proceeding.

### E-2 — Root finding by polynomial fitting across a 6-decade scale range
**Class:** Numerical uncertainty.
**Symptom.** The thermal root came out as `~1e-4 i` instead of `~(1+i)·2.6e6`, and
the "flux" of that spurious mode was wrong by 5 orders of magnitude. Check D1
failed with a 99.6 % relative deviation.
**Cause.** `det M` was sampled on a grid of ξ² and fitted with a degree-4
polynomial. The thermal wavenumber is ~2.6e6 m⁻¹ while the elastic wavenumbers are
~1e3 m⁻¹ and the electric one ~8.5e2 m⁻¹ — six decades. A least-squares fit on a
grid sized for the elastic scale cannot resolve the thermal root, and the fit
residual (0.8 %) was small enough to look acceptable.
**Fix.** Abandon fitting entirely. The five non-zero coefficients of `det M` are
obtained **exactly and symbolically** by Berkowitz fraction-free elimination,
lambdified once, and cached in `code/det_coeffs.json`. Roots then come from
`mpmath.polyroots` on the exact coefficients at 40 decimal digits.
**Lesson.** Never fit a polynomial to infer coefficients when the coefficients are
available exactly. The fitting residual is not a proxy for root accuracy.

### E-3 — `mpmath.polyroots` silently requires *descending* coefficients
**Class:** Numerical uncertainty.
**Symptom.** After E-2 was fixed, `polyroots` returned four roots clustered at
zero and lost the physical branches.
**Cause.** The coefficient list was assembled in ascending order
`[c₀, c₁, …, c₈]` (natural for `Σ cⱼ ξʲ`), but `mpmath.polyroots` expects
descending order. It does not warn.
**Fix.** `mp.polyroots(list(reversed(co)))`, with the polynomial made monic first.
**Lesson.** Verify any root routine against a known simple polynomial before
trusting it.

### E-4 — Normalising a mixed-unit state vector by `max|a|`
**Class:** Physical inconsistency.
**Symptom.** The flux antisymmetry `⟨P_z⟩(ξ) = −⟨P_z⟩(−ξ)` failed with a 200 %
error on the propagating branches of CdSe and PZT-5H, while Al6061 passed exactly.
**Cause.** The null vector `a = (U, W, Φ, Θ)` has components in **metres, volts and
kelvin**. Normalising by `max|a|` is dimensionally meaningless; it rescales the
four components inconsistently from one root to the next, so the same physical mode
got two different normalisations at +ξ and −ξ. Al6061 passed only because `Φ ≡ 0`
and `Θ ≡ 0` there, leaving a purely mechanical vector.
**Fix.** Fix the largest **mechanical** component to unity. The antisymmetry then
holds to 3.5e-13.
**Lesson.** Any normalisation of a multi-physics eigenvector must be checked for
dimensional consistency. This is a trap that would silently corrupt every
amplitudes-based result downstream.

### E-5 — Cofactor null vectors are not accurate enough at condition number 1e22
**Class:** Numerical uncertainty.
**Symptom.** Even with E-4 fixed, the flux conservation check drifted by 6e-12 for
CdSe while Al6061 achieved 4e-28. Increasing precision from 40 to 80 digits did
**not** help.
**Cause.** The polynomial root itself was converged to 1.2e-28. The loss was in the
null vector: the 6 mm matrix has condition number ≈1e22 in double precision, and a
4×4 cofactor expansion inherits that conditioning. A cofactor-based null vector is
accurate only to ε·cond(M), and raising the working precision does not remove a
*cancellation* in an ill-conditioned quantity.
**Fix.** Two stages: use the cofactors only to *select* the best 3×3 subsystem and
the column with the largest cofactor; then solve that 3×3 exactly by Cramer's rule
on the **column-equilibrated** matrix. Column equilibration lowers the effective
condition number to ~1e4. The remaining 6e-12 drift for CdSe is the cancellation
floor of a flux that is a difference of ~1e20 W m⁻² terms; the tolerance is set to
1e-10 and **the justification is stated in the code**, not chosen to pass.
**Lesson.** Report conditioning-limited tolerances as such. A number that "passes
once the tolerance is loosened" is worthless unless the reason for the floor is
identified.

### E-6 — Applying an *isotropic* closed form to a transversely isotropic material
**Class:** Mathematical error.
**Symptom.** Check I1 compared numerical roots against
`ξ = √((ω/c_L)² − k²)` with `c_L² = C₃₃/ρ` for PZT-5H and found up to 75 %
deviation at 75° incidence.
**Cause.** That closed form is valid only when `C₁₃ = C₁₁ − 2C₄₄`, i.e. when the
medium is elastically isotropic *in the sagittal plane*. For PZT-5H,
`C₁₃ = 84.1 GPa` but `C₁₁ − 2C₄₄ = 79.4 GPa` — it is genuinely transversely
isotropic, so its quasi-P and quasi-SV branches are angle-dependent and no
isotropic closed form applies. The *code* was right; the *check* was wrong.
**Fix.** Run the check on Al6061, which satisfies the isotropic identity exactly
(`C₁₃ = 54.6 = C₁₁ − 2C₄₄`), and assert the identity holds before trusting the
comparison. Result: 1.7e-15 over 24 cases.
**Lesson.** A validation check is only as good as its reference. Validate the
reference's preconditions too.

### E-7 — I1 initially reported `PASS` while executing **zero** cases
**Class:** Unsupported claim (in my own code).
**Symptom.** Check I1 printed `[PASS] ... max |dR| 0.00e+00`. The detail showed
`no incident longitudinal root` for all 12 parameter points.
**Cause.** The incident-mode search restricted itself to the *outgoing* root set,
which by construction excludes the incoming wave. Every case was skipped, the
violation accumulator stayed at zero, and the pass condition `worst < tol` was
trivially satisfied.
**Fix.** Two changes: (i) guard every check with `n_exec == 0 → UNRESOLVED`, never
`PASS`; (ii) search for the incident root in the *complement* of the outgoing set.
**Lesson.** **A check that executes nothing must never report success.** This is
the single most dangerous failure mode in automated validation, and the guard is
now applied to every check in the file.

### E-8 — Signed vs absolute power ratio in the P–SV reflection calculation
**Class:** Physical inconsistency.
**Symptom.** Reflection coefficients came out *negative* and `R + T − 1` was off by
0.1–0.66.
**Cause.** Two compounding errors: (i) the reflected branch carries flux in −z, so
its *signed* z-flux is negative; (ii) taking `|P_ref|` instead of `−P_ref` adds the
**interference** between the reflected P and S branches into R, which breaks
`R + T = 1` by exactly the cross term.
**Fix.** Use `R = −Re(P_ref)/Re(P_in)`, and define the reflected amplitude against
the incident amplitude rather than its flux.
**Lesson.** For multi-channel reflection, `R` is not a sum of `|r_i|²`. Channel-wise
sums are valid only where orthogonality has been shown — which is precisely what
Blueprint §5.1 warns about.

### E-9 — Claimed the CdSe parameter set was Lord–Shulman causal when it was not
**Class:** Unsupported claim / physical inconsistency.
**Symptom.** Check H1 reported `LS causal=True` for all four materials. But the
CdSe working value `τ_q = 1.0e-13 s` is **below** the derived causal floor
`τ_q,min = κ/(ρ c_E c_L²) = 4.67e-13 s`, so the LS thermal wave would be
superluminal.
**Cause.** The check printed the floor but did not *compare* against it; the
admissibility flag was computed only from `τ_q ≥ 0`.
**Fix.** (i) `τ_q = 5.0e-13 s` for CdSe, above the floor. (ii) The check now
computes `τ_q/τ_q,min` and fails if any material is below 1. CdSe sits at 1.07,
Al6061 at 1.21 — marginal, and honestly reported as such.
**Lesson.** Computing a constraint is not the same as enforcing it. Every derived
inequality must appear in the pass condition.

---

## Part 2 — Open items

### U-1 — Gate G0 (full-text prior-art comparison) NOT passed — `BLOCKING`
All ten records in `literature/verified_references.csv` have
`read_status = metadata only`. Full text was **not** accessed for any of them
(R03 and R05 had their abstracts read). Consequently:

- no compatibility check is possible for benchmarks B-12 to B-14;
- the "candidate gap" statements in `literature/research_gap_matrix.md` are
  questions for full-text inspection, **not** findings;
- **no novelty claim is made anywhere in Phase 1.**

This is the single largest limitation on the phase. It cannot be fixed by more
searching of metadata; it needs subscription access or interlibrary supply.

### U-2 — The impedance-matching literature was not searched — `BLOCKING` for G-3
Proposition 1 is elementary transmission-line algebra and is very likely folklore
in the impedance-matching / absorber literature, which was not searched in Phase 1.
Before any novelty wording about it, that literature must be checked. The
proposition is used here as a *derived lemma for benchmarking*, which is how
Blueprint §5.3 frames it, and that use is safe.

### U-3 — The coupled mode-matching scattering solve could not be validated — `BLOCKING`
**This is the main technical shortfall of Phase 1.** An 8×8 welded-interface
mode-matching system was assembled and the following were established:

- the operator has **rank 8** at 5°, 30° and 60° (check I2), with all eight
  singular values non-negligible after two-sided equilibration;
- the solve conserves energy: `R + T = 1` to 1e-5;
- the reflected and transmitted **power ratios** are of the right order of
  magnitude.

But the reflected **P amplitude** did not reproduce the Fresnel value: at
`Z₂/Z₁ = 4`, the matched `R` was 0.333 against a Fresnel value of 0.360, an error
of 2.7 %, and it did not decrease with λ or with θ. A rank-8, energy-conserving
system with a wrong amplitude split indicates a **normalisation** error — that is,
the mapping from the null vector `a = (U, W, Φ, Θ)` to a physical amplitude with
the right units — not a rank or assembly error.

Rather than tune the normalisation until the number matched, the check was
**withdrawn** and replaced by check I1, which validates the root solver against
external closed-form branches (24 cases, 1.7e-15). **No scattering coefficient is
reported anywhere in Phase 1.**

Additional structural facts found while attempting this, which any Phase-2 solver
must handle:

1. **`k = 0` with `e = 0` is degenerate.** `ξ = 0` becomes a root of multiplicity
   four and the potential is a pure gauge field. Normal incidence must be a limit
   `k → 0`, never a parameter value.
2. **The exactly decoupled case is singular.** With `e = β = p = 0` exactly, `φ`
   and `θ` appear in no other equation, so four of the eight interface conditions
   become identically zero rows. The decoupled case must be reached by scaling the
   couplings with `λ → 0`.
3. **Two-sided equilibration is mandatory.** In the 8×8 system the displacement
   rows are O(1) m and the heat-flux rows are O(1e19) W m⁻². Row-only scaling
   leaves a numerical null space of dimension 2 and the solve is meaningless.
   Rows *and* unknown columns must both be equilibrated.
4. **The condition number after equilibration is still ~1e9**, so double precision
   is insufficient; 40-digit arithmetic is required.

### U-4 — Pyroelectric and thermoelastic constants for PZT-5H and BaTiO₃ are ASSUMED — `BLOCKING` for any quantitative claim
`p₃`, `β₁` and `β₃` could not be sourced to a citable table for either material.
The values in `literature/material_constants.csv` are marked `ASSUMED` with an
order-of-magnitude basis. Only CdSe has tabulated `β₁`, `β₃` and `p₃`. **No
reported number may rest on the precise values of these constants**, and any
Phase-2 result depending on the thermal or pyroelectric channel must first obtain
sourced values.

### U-5 — Factor-of-two convention for the interface dissipation is not unified — `OPEN`
`D_mech = ½ ⟦v⟧ᵀC⟦v⟧` carries a ½ from the cycle average; `D_th = H_T|⟦θ⟧|²/T₀`
does not. Assembling a numerical energy balance with both as written produces a
50 % error in the thermal channel. The correct unification depends on the
normalisation adopted for `R`, `T`, `A` in the Phase-2 solver and is therefore
deliberately left open rather than guessed.

### U-6 — The `h → 0` reduction is established only to first order — `OPEN`
The composition `s(h⁺) = P₂ e^{hA} P₁ s(0⁻)` was expanded to first order, giving
an effective compliance = sum of the two face compliances + `O(h)`. The second-order
term, and whether a *single scalar* effective compliance exists for general
(non-diagonal) face laws, are not established. Blueprint §6 marks this as a
hypothesis and it is only partially confirmed.

### U-7 — The thermal mode's flux at practical evaluation planes is not quantified — `OPEN`
The thermal diffusion length is `μ = √(2κ/(ωρc_E)) ≈ 3.9e-7 m` at 1 MHz, so a
plane 10 μm from the interface is ~25 diffusion lengths away and should see no
thermal flux. This is argued, not demonstrated in a scattering configuration.
Blueprint §5.1 warns that evanescent and thermal modes can produce cross terms, so
this needs numerical confirmation before channel-wise sums are used.

### U-8 — The interphase configuration has not been analysed at all — `DEFERRED`
The primary candidate configuration (H2) — a finite interphase with passive faces
on both sides — has not been set up. Phase 1 established only that its governing
machinery is the same and that the `h → 0` limit is the *imperfect-contact*
baseline, not the welded one. The 16×16 system and its rank have not been built.

### U-9 — The bulk absorption integral is not separated from surface reflection — `OPEN`
Equation (master) of `theory/energy_flux_derivation.tex` contains
`∫ D_bulk dz` over the window between evaluation planes. For lossy media this is
non-zero and belongs to `A`. No configuration has yet been set up to evaluate it,
and the known danger is precisely the one the blueprint's rule F1 guards against:
a reflection-only optimiser certifies an attenuating window as anti-reflecting.
The derivation shows this is structural, not a reporting preference.

### U-10 — No §-level review by the supervisor has occurred — `DEFERRED` by instruction
Phase 1 stops here. Nothing in Phase 1 has been approved, and the blueprint's
stage-gate D and E are not addressed.

---

## Part 3 — Deliberately not done (per instruction)

| Item | Why not |
|---|---|
| Numerical solver for the coupled scattering problem | Phase 2 (P2) work; Phase 1 stops before it |
| Any anti-reflection curve, any optimised contact parameter, any parameter search over the H1/H2 boxes | Explicitly forbidden until the complete coupled solver exists |
| Manuscript text | Not in Phase-1 scope |
| Experimental, device, ML, topological, non-Hermitian, circuit content | Excluded by standing instruction |
| Gurtin–Murdoch surface elasticity, Klein–Gordon/NSGT nonlocality | Excluded from Phase 1; requires a logged demonstrated need |
| Compiling the `.tex` files | No TeX engine is installed in this environment. The files are text-verified only. **No compilation is claimed.** |

---

## Part 4 — Honest summary

**What Phase 1 established.** The bulk formulation, the exact determinant
structure, the root inventory, the flux functional and the scalar propositions are
derived and independently cross-checked. Thirteen automated checks pass with
measured residuals, and four of them compare against references that are external
to the code (Fresnel/closed-form elastic branches, the analytic elastic flux
`½ρc_Lω²|W|²`, the quarter-wave formula, and the analytic Proposition-1 values).

**What Phase 1 did not establish.** The scattering layer. No reflection or
transmission coefficient has been computed and validated, so per Blueprint §7.4
(energy gates block `R`, `T`, `A`) **Phase 1 reports none**. The prior-art
comparison did not reach full text, so no novelty claim is made.

**Nine of my own errors were found and corrected** (Part 1). Three of them
(E-7, E-4, E-9) would have produced a *silently wrong* result rather than a visible
failure: a check that passed without executing, a normalisation that corrupted
amplitudes, and a causality constraint that was computed but not enforced. Those
three are the ones most worth carrying forward.
